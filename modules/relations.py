import networkx as nx
import algebra as va
from .data import cache_get, cache_set
from .calling import sort_types
import warnings
from itertools import combinations

# TODO use consistent datastructure with OOP

def find_context(nodes, edges, directional=False, extend=False, extended=None, depth=0):
    """Find the context (connected nodes) for a given set of nodes based on an edge list."""
    # TODO show some information of indirect matches?
    # TODO make graph the input
    eq_graph = nx.Graph([(left, right) for left, right, relation in edges if relation == va.Relation.EQUIVALENT])
    cont_graph = nx.DiGraph([(left, right) for left, right, relation in edges if relation == va.Relation.IS_CONTAINED])
    overlap_graph = nx.Graph([(left, right) for left, right, relation in edges if relation == va.Relation.OVERLAP])
    graphs = {
        va.Relation.EQUIVALENT: eq_graph,
        va.Relation.IS_CONTAINED: cont_graph,
        va.Relation.OVERLAP: overlap_graph
    }
    context_nodes = set()
    context_nodes |= nodes
    context_edges = set()

    for node in nodes:	
        neighbour_nodes = set()
        for rel, graph in graphs.items():
            if node not in graph.nodes():
                continue
            if rel == va.Relation.IS_CONTAINED:
                for neighbour, _ in graph.in_edges(node):
                    if extend and neighbour in extended or sort_types(neighbour) == 4:
                        continue
                    neighbour_nodes.add(neighbour)
                    context_edges.add((neighbour, node, rel))
                if directional: # Single direction of containment
                    continue
            for neighbour in graph[node]:
                if extend and neighbour in extended or sort_types(neighbour) == 4:
                    continue
                neighbour_nodes.add(neighbour)
                context_edges.add((node, neighbour, rel))
        context_nodes |= neighbour_nodes
        # Extend suballeles
        if extend:
            for neighbour in neighbour_nodes:
                if sort_types(neighbour) not in (1, 2):
                    continue
                extended.add(neighbour)
                extended_nodes, extended_edges = find_context({neighbour,}, edges, directional, sort_types(neighbour) == 2, extended, depth+1)
                context_nodes |= extended_nodes
                context_edges |= extended_edges
    return context_nodes, context_edges


def redundant_reflexive(graph):
    """Returns redundant reflexive relations"""
    return nx.selfloop_edges(graph)

def redundant_common_ancestor(subgraph_contained, subgraph_overlap, full=False):
    """Returns redundant overlap relations due to common ancestor
    
    If the graph is full than the ancestor check is not needed since all relationships would exist and can be checked directly.
    However this function works more generally.
    """
    # TODO use ancestors function of nx instead and do all at once
    # TODO can do this faster with lowest common ancestors function (but cannot use nx implementation)
    to_remove = []
    ancestors = {}
    for node in subgraph_overlap.nodes():
        if node not in subgraph_contained.nodes():
            continue
        ancestors[node] = nx.ancestors(subgraph_contained, node)
    for s, t in subgraph_overlap.edges():
        if s not in ancestors or t not in ancestors:
            continue
        if ancestors[s] & ancestors[t] != set(): # Overlapping nodes have common ancestor
            to_remove.append((s, t))
    return to_remove

def redundant_transitive(graph):
    """Return edges redundant due to transitivity."""
    raise DeprecationWarning("Now using transitive reduction directly.")
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

def dfs_less_specific(current, overlap, subgraph_contained, subgraph_overlap):
    """Find if multiple edges in a path overlap with the same node and return the less specific."""
    edges = []
    at_depth = set()
    if current in subgraph_overlap.nodes():
        # Check if overlapping already in path
        for overlapping in subgraph_overlap[current]: # All overlapping at current depth
            if overlapping in overlap: # Already found in this path, is redundant
                edges.append((current, overlapping))
            else:
                at_depth.add(overlapping) # Track which added at this level
                overlap.add(overlapping) # Track overlapping at all depths
    # Go to next depth in containment relation
    for next in subgraph_contained[current]:
        edges += dfs_less_specific(next, overlap, subgraph_contained, subgraph_overlap)
    # forget at these depth when taking other path
    for overlapping in at_depth:
        overlap.remove(overlapping)
    return edges

def redundant_most_specific(subgraph_contained, subgraph_overlap):
    """Returns redundant overlap relations due to most specific (smallest contained) overlap
    
    Start at nodes that contain no nodes and do a DFS while storing the most specific overlap.
    Alternatively you can do the triangle approach, which is more efficient.
    You find all nodes with two overlaps and remove the overlap if one is contained in the other.
    """
    # TODO use triangle approach?
    to_remove = []
    for start_node, degree in subgraph_contained.in_degree(): # For all nodes that contain no nodes (start of path)
        if degree != 0:
            continue
        to_remove += dfs_less_specific(start_node, set(), subgraph_contained, subgraph_overlap)
        # TODO use non dfs approach?
        # for node in nx.dfs_successors(subgraph_contained, start_node): # Follow path of containment
        #     print(node)
        #     if node not in subgraph_overlap.nodes():
        #         continue
        #     for adjacent in subgraph_overlap[node]:
        #         if adjacent not in overlapping: # Most specific overlap
        #             overlapping.add(adjacent)
        #         else: # Already found, is redundant
        #             to_remove.append((node, adjacent))
    return to_remove

def redundant_symmetric(graph):
    """Return one direction of symmetric relations."""
    raise DeprecationWarning("Redundant because of undirected graph")
    to_remove = set()
    for s, t, d in graph.edges(data=True):
        if d["relation"].name not in ("OVERLAP", "EQUIVALENT"): # Only symmetric relations
            continue
        if (s, t) in to_remove: # Already removing one direction
            continue
        to_remove.add((t, s)) # Remove other direction
    return to_remove

def redundant_equivalence(graphs):
    """Remove redundant relations due to equivalence"""
    for component in nx.connected_components(graphs[va.Relation.EQUIVALENT]):
        if len(component) == 1: 
            continue
        # Remove duplicate relations
        # Prefer corealleles over suballeles, etc.
        component = sorted(list(component), key=sort_types)
        for rel in (va.Relation.EQUIVALENT, va.Relation.IS_CONTAINED, va.Relation.OVERLAP):
            if component[0] not in graphs[rel].nodes():
                continue
            for other in component[1:]:
                if other not in graphs[rel].nodes():
                    continue
                graphs[rel] = nx.contracted_nodes(graphs[rel], component[0], other, self_loops=False)
        # Add removed relations back
        graphs[va.Relation.EQUIVALENT].add_edges_from([(component[0], other) for other in component[1:]])
    return graphs

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
    # TODO allow for simplifying already simplified relations
    try:
        if cache_name: return cache_get(cache_name)
    except:
        pass
    # Construct networkx graph objects from edge list
    # Using undirected graphs reduces symmetric redundancy
    subgraph_equivalence = nx.Graph([(s, t) for s, t, d in relations if d.name == "EQUIVALENT"])
    subgraph_overlap = nx.Graph([(s, t) for s, t, d in relations if d.name == "OVERLAP"])
    subgraph_contained = nx.DiGraph([(s, t) for s, t, d in relations if d.name == "IS_CONTAINED"])
    # Remove reflexive self-loops
    subgraph_equivalence.remove_edges_from(redundant_reflexive(subgraph_equivalence))
    # Remove transitive redundancies 
    subgraph_contained = nx.transitive_reduction(subgraph_contained)
    # Remove common ancestor redundancy
    subgraph_overlap.remove_edges_from(redundant_common_ancestor(subgraph_contained, subgraph_overlap)) 
    # Remove most specific redundancy
    subgraph_overlap.remove_edges_from(redundant_most_specific(subgraph_contained, subgraph_overlap))
    # Remove redundant relations due to equivalence (in place)
    graphs = {
        va.Relation.EQUIVALENT: subgraph_equivalence, 
        va.Relation.OVERLAP: subgraph_overlap,
        va.Relation.IS_CONTAINED: subgraph_contained
        # Only store is contained side of relation
    }
    graphs = redundant_equivalence(graphs)
    # Convert back to edge list
    nodes = set()
    edges = []
    for rel, graph in graphs.items():
        nodes |= set(graph.nodes())
        for s, t in graph.edges():
            edges.append((s, t, rel))
    if cache_name: cache_set((nodes, edges), cache_name)
    return nodes, edges    