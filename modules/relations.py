import algebra as va
import igraph as ig
import networkx as nx
from .data import cache_get, cache_set
from .parse import parse_multi_hgvs
from .va_tools import count_relations

def find_relations(corealleles, reference_sequence):
    """Find the relation between all corealleles.

    Relations are cached since they take a long time to generate.

    Returns edge list.
    """
    # TODO use name of enum
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


def prune_relations(nodes, relations):
    """Prune relations which are redundant.

    Symmetric relations should not be displayed twice (disjoint, equivalent, overlap).
    Reflexive relations should not be displayed (equivalence to self).
    Transitive relations can be reduced (equivalence, containment).
    A combination of containment and overlap can be used to filter out some relations:
    - Most specific: if A --> B and A -- C; B -- C then B -- C is redundant
    - Common ancestor: if A --> B; A --> C and B -- C then B -- C is redundant
    Disjoint relation can be left out since they are understood as 'no relation'.
    One direction of the containment relation can be left out.

    returns list of edges.
    """
    for relation, count in count_relations(relations).items():
        print(relation, count)

    # Create edge list in proper format
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

    # Remove reflexive self-loops
    graph.remove_edges_from(nx.selfloop_edges(graph))
    # Remove common ancestor redundancy
    # TODO
    # Remove most specific redundancy
    # TODO
    # Remove one direction of containment
    graph.remove_edges_from([(s, t) for s, t, d in graph.edges(data=True) if d["relation"].name == "CONTAINS"])
    # Remove transitive redundancies for single relations TODO include equivalence
    for relation in ("IS_CONTAINED",):
        subgraph = nx.DiGraph([(s, t, d) for s, t, d in graph.edges(data=True) if d["relation"].name == relation])
        subgraph = nx.transitive_reduction(subgraph) # Reduce edges which are redundant due to transitivity of the same relation
        graph.remove_edges_from([(s, t) for s, t, d in graph.edges(data=True) if d["relation"].name == relation])
        graph.add_edges_from([(s, t, {"relation": va.Relation[relation]}) for s, t in subgraph.edges()])
    # Remove symmetric relations by treating graph as undirected
    graph = graph.to_undirected()

    # Convert back to edge list to be used elsewhere
    edges = [(s, t, d["relation"]) for s, t, d in graph.edges(data=True)]

    print()
    for relation, count in count_relations(edges).items():
        print(relation, count)
    return edges