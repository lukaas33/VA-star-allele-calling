import networkx as nx
from .data import cache_get, cache_set
from .parse import parse_hgvs_supremal
from .utils import printProgressBar
import warnings
import algebra as va
import multiprocessing as mp
from multiprocessing.managers import SharedMemoryManager
from textwrap import wrap
from sys import getsizeof
import math

def find_context(nodes, edges):
    """Find the context (connected nodes) for a given set of nodes based on an edge list."""
    # TODO do this based on a networkx graph
    context = set()
    for node in nodes:
        for s, t, d in edges:
            if s == node or t == node:
                context.add(s)
                context.add(t)
    return context

def find_relation(args):
    """Worker for multiprocessing relations."""
    # TODO make faster
    # TODO fix progress bar resetting (integer overflow?)
    left, right, ref_chunks, sequences, count = args
    # Find relation for reconstructed variants
    reference = "".join(ref_chunks)
    lhs = va.Variant(sequences[left*3], sequences[left*3+1], sequences[left*3+2])
    rhs = va.Variant(sequences[right*3], sequences[right*3+1], sequences[right*3+2])
    relation = va.relations.supremal_based.compare(reference, lhs, rhs)
    # Print progress
    count[0] += 1
    alleles = len(sequences)//3
    printProgressBar(count[0], (alleles**2 + alleles)/2, prefix = 'Comparing:', length = 50)
    return left, right, relation

def find_relations_all(corealleles, reference_sequence, suballeles=None):
    """Find the relation between all corealleles, suballeles and variants.

    Relations are cached since they take a long time to generate.

    Returns edge list.
    """
    cache_name = "all_relations"
    if suballeles is not None: cache_name += "_sub"
    # Test if already stored
    try:
        return cache_get(cache_name)
    except:
        pass

    # Get all variants as a dictionaries of supremal strings
    all_variants = {} # Store variants, sub- and corealleles as lists of HGVS
    for coreallele in corealleles.keys():
        alleles = [corealleles[coreallele]]
        if suballeles is not None: # Include sub
            alleles += suballeles[coreallele]
        for allele in alleles:
            all_variants[allele["alleleName"]] = [] 
            for variant in allele["variants"]: # Variants for allele
                all_variants[allele["alleleName"]].append(variant["hgvs"])
                if variant["hgvs"] in all_variants.keys(): # Skip if already parsed
                    continue
                all_variants[variant["hgvs"]] = parse_hgvs_supremal([variant["hgvs"]], reference_sequence) # Store variant as supremal
            try:
                all_variants[allele["alleleName"]] = parse_hgvs_supremal(all_variants[allele["alleleName"]], reference_sequence) # Store allele as supremal
            except ValueError as e: # Fails for overlapping variants
                # TODO how to handle duplicates? And how to handle multiple variants at same position?
                error = f"{allele['alleleName']}: {e}"
                warnings.warn(error)
                del all_variants[allele["alleleName"]]

    # Parse and get relations
    relations = []
    variant_names = list(all_variants.keys())
    print(len(variant_names), "variants")
    with SharedMemoryManager() as smn:
        # TODO is there a better way to do this?
        # Store data between processes
        # Divide into chunks to allow storage in shared memory (10MB max)
        size = getsizeof(reference_sequence)   
        chunks = wrap(reference_sequence, size // math.ceil(size / 10e6))
        ref = smn.ShareableList(chunks)
        # Spread properties since variant can't be stored in shared memory
        spread = [] 
        for supremal in all_variants.values():
            spread += [supremal.start, supremal.end, supremal.sequence]
        seqs = smn.ShareableList(spread)
        count = smn.ShareableList([0])
        # Multiprocessing of relations
        with mp.Pool(mp.cpu_count()) as pool:
            # TODO change to strings
            args = ((i, j, ref, seqs, count) for i in range(len(variant_names)) for j in range(i, len(variant_names)))
            relation_pairs = pool.map(find_relation, args)
            # Store relations
            for left, right, relation in relation_pairs:
                left, right = variant_names[left], variant_names[right]
                if relation == va.Relation.CONTAINS:
                    inv_relation = va.Relation.IS_CONTAINED
                elif relation == va.Relation.IS_CONTAINED:
                    inv_relation = va.Relation.CONTAINS
                else:
                    inv_relation = relation
                relations.append((left, right, relation))
                relations.append((right, left, inv_relation))

    cache_set(relations, cache_name)
    return relations

def has_common_ancestor(graph, node1, node2):
    """Check if two nodes have a common ancestor in a directed graph."""
    if node1 not in graph.nodes() or node2 not in graph.nodes():
        return False
    # Do a parallel BFS on the directed graph and check if both have a node as a child
    # TODO possible to not do a complete BFS the first time?
        # TODO can reduce for complete graph but doesn't work in all instances
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

def redundant_common_ancestor(subgraph_contains, subgraph_overlap):
    """Returns redundant overlap relations due to common ancestor
    
    If the graph is full than the ancestor check is not needed since all relationships would exist.
    But this function can be used more generally.
    """
    # TODO use ancestors function of nx instead and do all at once
    # TODO can do this faster with lowest common ancestors function (but cannot use nx implementation)
    to_remove = []
    for s, t in subgraph_overlap.edges():
        if has_common_ancestor(subgraph_contains, s, t):
            to_remove.append((s, t))
    return to_remove

def redundant_transitive(graph):
    """Return edges redundant due to transitivity."""
    to_remove = []
    for relation in (va.Relation.IS_CONTAINED, va.Relation.EQUIVALENT, va.Relation.CONTAINS):
        subgraph = graph.edge_subgraph([edge for edge in graph.edges() if edge not in redundant_symmetric(graph)])
        subgraph = subgraph.edge_subgraph([(s, t) for s, t, d in graph.edges(data=True) if d["relation"] == relation])
        subgraph_reduced = nx.transitive_reduction(subgraph)
        to_remove += [
            (s, t, d) 
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

def redundant_equivalence(graph):
    """Remove redundant relations due to equivalence"""
    # TODO can use contracted nodes method of nx 
    # TODO can consider equivalent nodes as a component
    to_remove = []
    for s, t, d in graph.edges(data=True):
        if d["relation"].name != "EQUIVALENT":
            continue
        s_conn = set(graph[s])
        t_conn = set(graph[t])
        share_conn = s_conn & t_conn
        for shared in share_conn:
            if "*" in s: # Favour allele relations (choice arbitrary)
                to_remove.append((t, shared))
                to_remove.append((shared, t))
            else:
                to_remove.append((s, shared))
                to_remove.append((shared, s))
    return to_remove 

def prune_relations(relations):
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
    # TODO make faster
    # Cache since it can take some time
    cache_name = "pruned_relations"
    try:
        return cache_get(cache_name)
    except:
        pass

    # Create edge list in proper format for networkx
    # Filter out disjoint relations
    nodes = set()
    edges = []
    for left, right, relation in relations:
        nodes.add(left)
        nodes.add(right)
        if relation.name == "DISJOINT":
            continue
        edges.append((left, right, {"relation": relation}))
    # Construct networkx graph object from edge list
    graph = nx.DiGraph()
    graph.add_nodes_from(nodes)
    graph.add_edges_from(edges)
    # Create subgraphs for some reductions (are linked to graph)
    subgraph_contains = graph.edge_subgraph([(s, t) for s, t, d in graph.edges(data=True) if d["relation"].name == "CONTAINS"])
    subgraph_overlap = graph.edge_subgraph([(s, t) for s, t, d in graph.edges(data=True) if d["relation"].name == "OVERLAP"])
    subgraph_contained = graph.edge_subgraph([(s, t) for s, t, d in graph.edges(data=True) if d["relation"].name == "IS_CONTAINED"])
    # Remove reflexive self-loops
    graph.remove_edges_from(redundant_reflexive(graph))
    print("Reflexive removed")
    # Remove common ancestor redundancy
    graph.remove_edges_from(redundant_common_ancestor(subgraph_contains, subgraph_overlap)) # Long running
    print("Common ancestor removed")
    # Only one direction of containment will be displayed so it won't be filtered here
    # Remove transitive redundancies for single relations 
    graph.remove_edges_from(redundant_transitive(graph)) # Long running
    print("Transitive removed")
    # Remove most specific redundancy
    graph.remove_edges_from(redundant_most_specific(subgraph_contained, subgraph_overlap))
    print("Most specific removed")
    # Remove redundant symmetric relations 
    graph.remove_edges_from(redundant_symmetric(graph))
    print("Symmetric removed")
    # Remove redundant relations due to equivalence
    graph.remove_edges_from(redundant_equivalence(graph))
    print("Equivalence removed")

    # Convert back to edge list
    edges = [(s, t, d["relation"]) for s, t, d in graph.edges(data=True)]
    cache_set((nodes, edges), cache_name)
    return nodes, edges    