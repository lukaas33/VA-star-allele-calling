import networkx as nx
import algebra as va
from .data import cache_get, cache_set
from .calling import find_type, Type, distance
import warnings
from itertools import combinations

# TODO use consistent datastructure with OOP

def find_context(nodes, edges, directional=False, extend=False, extended=None, overlap=True, depth=0):
    """Find the context (connected nodes) for a given set of nodes based on an edge list.
    
    Direct context is the nodes directly connected to the input.
    Extended context also includes edges while the relation becomes stronger so you would show equivalent nodes of a contained node but not overlaps.

    Doesn't include samples.
    """
    # TODO make graph the input
    # TODO make more efficient
    eq_graph = nx.Graph([(left, right) for left, right, relation in edges if relation == va.Relation.EQUIVALENT])
    cont_graph = nx.DiGraph([(left, right) for left, right, relation in edges if relation == va.Relation.IS_CONTAINED])
    overlap_graph = nx.Graph([(left, right) for left, right, relation in edges if relation == va.Relation.OVERLAP])
    graphs = {
        va.Relation.EQUIVALENT: eq_graph,
        va.Relation.IS_CONTAINED: cont_graph,
        va.Relation.OVERLAP: overlap_graph
    }
    context_nodes = set()
    context_edges = set()
    context_nodes |= nodes
    if extend: extended |= nodes
    # Depth limit for debugging
    if depth == 99:
        return context_nodes, context_edges
    # Find context for each node
    for node in nodes:
        for rel, graph in graphs.items(): 
            if node not in graph.nodes():
                continue
            if not overlap and rel == va.Relation.OVERLAP:
                continue
            direct = set()
            # Find direct context of node
            # Nodes
            if rel == va.Relation.IS_CONTAINED:
                direct |= set([n for n in graph.predecessors(node) if find_type(n) != Type.SAMPLE])
                if not directional:
                    direct |= set([n for n in graph.successors(node) if find_type(n) != Type.SAMPLE])
            else:
                direct |= set([n for n in graph[node] if find_type(n) != Type.SAMPLE])
            # Edges
            for n in direct:
                if find_type(n) == Type.SAMPLE:
                    continue
                if extend and n in extended:
                    continue
                if rel == va.Relation.IS_CONTAINED:
                    context_edges.add((n, node, rel))
                    if directional: # Single direction of containment
                        continue
                    rel = va.Relation.CONTAINS
                context_edges.add((node, n, rel))
            if not extend:
                context_nodes |= direct
                continue
            # Find extended context of node
            # Only allow overlap when relation is not stronger.
            # Equivalence will always be a terminal node
            novel = direct - context_nodes - extended
            context_nodes |= direct
            for n in novel:
                deeper_nodes, deeper_edges = find_context({n,}, edges, directional, extend, extended | direct, rel == va.Relation.OVERLAP, depth+1)
                context_nodes |= deeper_nodes
                context_edges |= deeper_edges
            
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
    # to_remove = [(s, t) for s, t in to_remove if find_type(s) != 4 and find_type(t) != 4] # Keep relations with samples
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
    # to_remove = [(s, t) for s, t in to_remove if find_type(s) != Type.SAMPLE and find_type(t) != Type.SAMPLE] # Keep relations with samples
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
        component = sorted(list(component), key=find_type)
        for rel in (va.Relation.EQUIVALENT, va.Relation.IS_CONTAINED, va.Relation.OVERLAP):
            if component[0] not in graphs[rel].nodes():
                continue
            for other in component[1:]:
                # TODO keep relations with samples
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

def find_path(s, t, cont_graph, eq_graph, overlap_graph, path=None, visited=None):
    """Find a path from s to t in the graphs"""
    # TODO fix
    raise NotImplementedError("Implementation contains bugs")
    # TODO integrate in web interface
    if path is None: path = [s]
    if visited is None: visited = set()
    if s in visited: return None
    visited.add(s)
    if s == t:
        return path
    for g in (cont_graph, eq_graph, overlap_graph):
        if s in g.nodes():
            for n in g[s]:
                if sort_types(n) in (4, 5): # Don't use samples for iteration
                    continue
                path.append(n)
                result = find_path(n, t, cont_graph, eq_graph, overlap_graph, path, visited)
                if result is not None:
                    return result
                path.pop()
    return None

def find_all_distances(relations_samples_extended, supremal_samples, supremal_extended, samples, reference_sequence, cache_name=None):
    def patch(sup1, sup2):
        if sup2 is None: # CYP2D6*1
            sup2 = va.Variant(sup1.start, sup1.end, reference_sequence[sup1.start:sup1.end])
        # Patch changes in area size of sup1
        # Needed for consistency of scores, we want to compare the score for the same observed allele
        interval = (sup1.start, sup1.end)
        seq = reference_sequence[interval[0]:interval[1]]
        seq1 = seq[:sup1.start-interval[0]] + sup1.sequence + seq[sup1.end-interval[1]+1:]
        seq2 = seq[:sup2.start-interval[0]] + sup2.sequence + seq[sup2.end-interval[1]+1:]
        return seq1, seq2
    try:
        if cache_name: return cache_get(cache_name)
    except:
        pass
    distances = {}
    for sample in samples:
        if sample.split('_')[1] != "all":
            continue
        distances[sample] = {}
        # For all contained and equivalent alleles to observed allele
        for left, right, rel in relations_samples_extended:
            if left != sample:
                continue
            if find_type(right) not in (Type.CORE, Type.SUB):
                continue
            if rel != va.Relation.CONTAINS and rel != va.Relation.EQUIVALENT:
                continue
            # Pairwise alignment and scoring
            distances[sample][right] = distance(*patch(supremal_samples[sample], supremal_extended[right]))
            print(left, right, distances[sample][right])
        distances[sample]["CYP2D6*1"] = distance(*patch(supremal_samples[sample], None))
    if cache_name: cache_set(distances, cache_name)
    return distances