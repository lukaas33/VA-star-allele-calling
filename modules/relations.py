import algebra as va
import networkx as nx
from .data import cache_get, cache_set
from .parse import parse_multi_hgvs

def find_relations(corealleles, reference_sequence):
    """Find the relation between all corealleles.

    Relations are cached since they take a long time to generate.

    Returns edge list.
    """
    # Find relation for each pair of corealleles directionally
    coreallele_names = list(corealleles.keys())
    name = "relations"
    try:
        edges = cache_get(name)
    except:
        relations = {
            coreallele: {
                coreallele2: None 
                for coreallele2 in coreallele_names
            } 
            for coreallele in coreallele_names
        }
        for i, left_coreallele in enumerate(coreallele_names):
            for right_coreallele in coreallele_names[i:]: # Only check each pair once
                # Parse allele
                left_hgvs = [variant["hgvs"] for variant in corealleles[left_coreallele]["variants"]]
                left_variants = parse_multi_hgvs(left_hgvs, reference_sequence, left_coreallele)
                right_hgvs = [variant["hgvs"] for variant in corealleles[right_coreallele]["variants"]]
                right_variants = parse_multi_hgvs(right_hgvs, reference_sequence, right_coreallele)
                # Store relation
                relation = va.compare(reference_sequence, left_variants, right_variants)
                relations[left_coreallele][right_coreallele] = relation
                # Store inverse relation (the same for most relations)
                if relation == va.Relation.CONTAINS:
                    inv_relation = va.Relation.IS_CONTAINED
                elif relation == va.Relation.IS_CONTAINED:
                    inv_relation = va.Relation.CONTAINS
                else:
                    inv_relation = relation
                relations[right_coreallele][left_coreallele] = inv_relation
        edges = []
        for l_node in corealleles:
            for r_node in corealleles:
                edges.append((l_node, r_node, relations[l_node][r_node]))
        cache_set(edges, name) # Cache
    return edges

def has_common_ancestor(graph, node1, node2):
    """Check if two nodes have a common ancestor in a directed graph."""
    if node1 not in graph.nodes() or node2 not in graph.nodes():
        return False
    # Do a parallel BFS on the directed graph and check if both have a node as a child
    # TODO possible to not do a complete BFS the first time?
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

def redundant_common_ancestor(subgraph_contains, subgraph_overlap):
    """Returns redundant overlap relations due to common ancestor"""
    to_remove = []
    for s, t in subgraph_overlap.edges():
        if has_common_ancestor(subgraph_contains, s, t):
            to_remove.append((s, t))
    return to_remove

def redundant_transitive(graph):
    """Return edges redundant due to transitivity."""
    to_remove = []
    for relation in (va.Relation.IS_CONTAINED, va.Relation.EQUIVALENT):
        subgraph = graph.edge_subgraph([(s, t) for s, t, d in graph.edges(data=True) if d["relation"] == relation])
        subgraph_reduced = nx.transitive_reduction(subgraph)
        to_remove += [
            (s, t, d) 
            for s, t, d in subgraph.edges(data=True) 
            if not subgraph_reduced.has_edge(s, t)
        ]
    return to_remove

def redundant_most_specific(subgraph_contained, subgraph_overlap):
    """Returns redundant overlap relations due to most specific (smallest contained) overlap"""
    to_remove = []
    for start_node in nx.topological_sort(subgraph_contained):
        # Do BFS from start of topological ordering
        # store overlapping, if already found the new one is less specific and can be removed
        # TODO don't repeat already visited nodes
        overlapping = set()
        queue = [start_node] 
        while len(queue) > 0:
            current = queue.pop(0)
            if current in subgraph_overlap.nodes():
                overlaps = [edge for edge in subgraph_overlap[current]] # Find overlap of current node 
                for target in overlaps:
                    if target in overlapping: # More specific overlap was known
                        to_remove.append((current, target))
                        to_remove.append((target, current))
                    else: # Add since this is a new overlap
                        overlapping.add(target)
            for neighbour in subgraph_contained[current]:
                queue.append(neighbour)
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

def prune_relations(nodes, relations):
    """Prune relations which are redundant.

    Symmetric relations should not be displayed twice (disjoint, equivalent, overlap).
    Reflexive relations should not be displayed (equivalence to self).
    Transitive relations can be reduced (equivalence, containment).
    A combination of containment and overlap can be used to filter out some relations:
    - Most specific: if A --> B and A -- C; B -- C then B -- C is redundant
    - Common ancestor: if A --> B; A --> C and B -- C then B -- C is redundant
    Disjoint relation can be left out since they are understood as 'no relation'.
    One direction of the containment relation can be left out since the inverse follows.

    returns list of edges.
    """
    # Create edge list in proper format for networkx
    # Filter out disjoint relations
    edges = [
        (edge[0], edge[1], {"relation": edge[2]}) 
        for edge in relations 
        if edge[2].name not in ("DISJOINT",)
    ]
    # Construct networkx graph object from edge list
    graph = nx.DiGraph()
    graph.add_nodes_from(nodes)
    graph.add_edges_from(edges)
    # Create subgraphs for some reductions (are linked to graph)
    subgraph_contains = graph.edge_subgraph(
        [(s, t) for s, t, d in graph.edges(data=True) if d["relation"].name == "CONTAINS"]
    )
    subgraph_overlap = graph.edge_subgraph(
        [(s, t) for s, t, d in graph.edges(data=True) if d["relation"].name == "OVERLAP"]
    )
    subgraph_contained = graph.edge_subgraph(
        [(s, t) for s, t, d in graph.edges(data=True) if d["relation"].name == "IS_CONTAINED"]
    )

    # Remove reflexive self-loops
    graph.remove_edges_from(redundant_reflexive(graph))
    # Remove common ancestor redundancy
    graph.remove_edges_from(redundant_common_ancestor(subgraph_contains, subgraph_overlap))
    # Remove one direction of containment
    graph.remove_edges_from(list(subgraph_contains.edges()))
    # Remove transitive redundancies for single relations 
    graph.remove_edges_from(redundant_transitive(graph))
    # Remove most specific redundancy
    graph.remove_edges_from(redundant_most_specific(subgraph_contained, subgraph_overlap))
    # Remove redundant symmetric relations 
    graph.remove_edges_from(redundant_symmetric(graph))

    # Convert back to edge list
    edges = [(s, t, d["relation"]) for s, t, d in graph.edges(data=True)]
    return edges