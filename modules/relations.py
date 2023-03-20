import networkx as nx
import algebra as va
from .data import cache_get, cache_set
from .calling import sort_types

# TODO use consistent datastructure with OOP

def find_context(nodes, edges, as_edges=False):
    """Find the context (connected nodes) for a given set of nodes based on an edge list."""
    # TODO do this based on a networkx graph
    context = set()
    context_edges = list()
    for node in nodes:
        for s, t, d in edges:
            if s == node or t == node:
                if as_edges:
                    context_edges.append((s, t, d))
                else:
                    context.add(s)
                    context.add(t)
    if as_edges:
        return context_edges
    return context

def has_common_ancestor(graph, node1, node2):
    """Check if two nodes have a common ancestor in a directed graph."""
    if node1 not in graph.nodes() or node2 not in graph.nodes():
        return False
    # Do a parallel BFS on the directed graph and check if both have a node as a child
    # TODO possible to not do a complete BFS the first time?
    # TODO use networkx ancestors
    visited = {node: [False, False] for node in graph.nodes()}
    for i, start_node in enumerate((node1, node2)):
        queue = [start_node]
        while len(queue) > 0:
            current = queue.pop(0) # next node
            if visited[current][i]: # Already visited from this start node
                continue
            visited[current][i] = True
            if all(visited[current]): # Visited from both start nodes
                return True
            for neighbour in graph[current]:
                queue.append(neighbour)
    return False

def redundant_reflexive(graph):
    """Returns redundant reflexive relations"""
    return nx.selfloop_edges(graph)

def redundant_common_ancestor(subgraph_contains, subgraph_overlap, full=False):
    """Returns redundant overlap relations due to common ancestor
    
    If the graph is full than the ancestor check is not needed since all relationships would exist and can be checked directly.
    This function can be used more generally and be used faster for full graphs.
    """
    # TODO use ancestors function of nx instead and do all at once
    # TODO can do this faster with lowest common ancestors function (but cannot use nx implementation)
    to_remove = []
    for s, t in subgraph_overlap.edges():
        if full:
            if subgraph_contains.has_node(s) and subgraph_contains.has_node(t) and \
                    set(subgraph_contains[s]) & set(subgraph_contains[t]):
                to_remove.append((s, t))
            continue
        if has_common_ancestor(subgraph_contains, s, t):
            to_remove.append((s, t))
    return to_remove

def redundant_transitive(graph):
    """Return edges redundant due to transitivity."""
    to_remove = []
    for relation in (va.Relation.IS_CONTAINED, ):
        subgraph = graph.edge_subgraph([edge for edge in graph.edges() if edge not in redundant_symmetric(graph)])
        subgraph = subgraph.edge_subgraph([(s, t) for s, t, d in graph.edges(data=True) if d["relation"] == relation])
        subgraph_reduced = nx.transitive_reduction(subgraph)
        to_remove += [
            (s, t, d) 
            for s, t, d in subgraph.edges(data=True) 
            if not subgraph_reduced.has_edge(s, t)
        ]
        to_remove += [
            (t, s, d) 
            for s, t, d in subgraph.edges(data=True) 
            if not subgraph_reduced.has_edge(s, t)
        ]
    return to_remove

def dfs(subgraph_contained, subgraph_overlap, current, overlapping):
    """Find if multiple edges in a path overlap with the same node and return the less specific."""
    edges = []
    if current not in subgraph_overlap.nodes(): # Skip
        return edges
    # Check if overlapping already in path
    for connected in subgraph_overlap[current]: # Track overlapping at current depth
        if any((connected in overlap for overlap in overlapping.values())): # Already found in this path, is redundant
            edges.append((current, connected))
            edges.append((connected, current))
    overlapping[current] = set(subgraph_overlap[current])
    # Go to next depth in containment relation
    for overlap in subgraph_contained[current]:
        edges += dfs(subgraph_contained, subgraph_overlap, overlap, overlapping)
    # forget when taking other path
    del overlapping[current]
    return edges

def redundant_most_specific(subgraph_contained, subgraph_overlap):
    """Returns redundant overlap relations due to most specific (smallest contained) overlap"""
    to_remove = []
    for start_node in nx.topological_sort(subgraph_contained):
        # Do DFS from start of topological ordering
        # store overlapping, if already found the new one is less specific and can be removed
        # TODO don't repeat already visited nodes (triangle approach)
        #       only have to check overlapping nodes where one has a containment and contains relation
        overlapping = dict()
        to_remove += dfs(subgraph_contained, subgraph_overlap, start_node, overlapping)
    return to_remove

def redundant_symmetric(graph):
    """Return one direction of symmetric relations."""
    to_remove = set()
    for s, t, d in graph.edges(data=True):
        if d["relation"].name not in ("OVERLAP", "EQUIVALENT"): # Only symmetric relations
            continue
        if (s, t) in to_remove: # Already removing one direction
            continue
        to_remove.add((t, s)) # Remove other direction
    return to_remove

def redundant_equivalence(subgraph_contained, subgraph_contains, subgraph_overlap, subgraph_equivalence):
    """Remove redundant relations due to equivalence"""
    # TODO can use contracted nodes method of nx 
    to_remove = []
    for component in nx.weakly_connected_components(subgraph_equivalence):
        if len(component) == 1: 
            continue
        # Favour relations with core allele, then sub, etc. (is arbitrary)
        center = sorted(list(component), key=sort_types)[0] 
        for node in component:
            if node == center: # Keep relations from here
                continue
            # Remove equivalence not connected to center
            for other in subgraph_equivalence[node]: 
                if other == center:
                    continue
                to_remove.append((node, other))
                to_remove.append((other, node))
            # Remove all other relations not with center
            # TODO traverse contained in reverse?
            for subgraph in (subgraph_contained, subgraph_overlap, subgraph_contains):
                if not subgraph.has_node(node):
                    continue
                for other in subgraph[node]:
                    to_remove.append((node, other))
                    to_remove.append((other, node))		
    return to_remove 

def prune_relations(relations, cache_name=None):
    """Prune relations which are redundant.

    Symmetric relations should not be displayed twice (disjoint, equivalent, overlap).
    Reflexive relations should not be displayed (equivalence to self).
    Transitive relations can be reduced (equivalence, containment).
    A combination of containment and overlap can be used to filter out some relations:
    - Most specific: if A --> B and A -- C; B -- C then B -- C is redundant
    - Common ancestor: if A --> B; A --> C and B -- C then B -- C is redundant
    Disjoint relation can be left out since they are understood as 'no relation'.
    One direction of the containment relation can be left out since the inverse follows.
    For equivalence one relation can be left out: A == B and A - C; B - C one is redundant. The choice of which is left out is partially arbitrary dependent on what is clearest in the context.

    returns list of edges and nodes.
    """
    # TODO redo more efficiently by using un/directed graphs
    try:
        if cache_name: return cache_get(cache_name)
    except:
        pass
    # Create edge list in proper format for networkx
    nodes = set()
    edges = []
    for left, right, relation in relations:
        nodes.add(left)
        nodes.add(right)
        if relation.name == "DISJOINT": # Same as no relation
            continue
        if relation.name == "CONTAINS": # Only one direction of containment is needed
            continue
        edges.append((left, right, {"relation": relation}))
    # Construct networkx graph object from edge list
    graph = nx.DiGraph()
    graph.add_nodes_from(nodes)
    graph.add_edges_from(edges)
    # Create subgraphs for some reductions (are linked to graph)
    subgraph_equivalence = graph.edge_subgraph([(s, t) for s, t, d in graph.edges(data=True) if d["relation"].name == "EQUIVALENT"])
    subgraph_overlap = graph.edge_subgraph([(s, t) for s, t, d in graph.edges(data=True) if d["relation"].name == "OVERLAP"])
    subgraph_contained = graph.edge_subgraph([(s, t) for s, t, d in graph.edges(data=True) if d["relation"].name == "IS_CONTAINED"])
    subgraph_contains = nx.DiGraph([(t, s, {"relation": va.Relation.CONTAINS}) for s, t, d in graph.edges(data=True) if d["relation"].name == "IS_CONTAINED"]) # TODO use incoming nodes for is contained
    # Remove reflexive self-loops
    graph.remove_edges_from(redundant_reflexive(graph))
    # Remove common ancestor redundancy
    graph.remove_edges_from(redundant_common_ancestor(subgraph_contains, subgraph_overlap)) 
    # Remove transitive redundancies for single relations 
    graph.remove_edges_from(redundant_transitive(graph)) 
    # Remove most specific redundancy
    graph.remove_edges_from(redundant_most_specific(subgraph_contained, subgraph_overlap))
    # Remove redundant relations due to equivalence
    graph.remove_edges_from(redundant_equivalence(subgraph_contained, subgraph_contains, subgraph_overlap, subgraph_equivalence))
    # Remove redundant symmetric relations 
    graph.remove_edges_from(redundant_symmetric(graph))

    # Convert back to edge list
    edges = [(s, t, d["relation"]) for s, t, d in graph.edges(data=True)]
    if cache_name: cache_set((nodes, edges), cache_name)
    return nodes, edges    