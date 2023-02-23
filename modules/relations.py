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
    # TODO use networkx as main library?
    for relation, count in count_relations(relations).items():
        print(relation, count)
    # Construct graph object from edge list
    graph = ig.Graph(
        n=len(nodes), 
        vertex_attrs={
            "name": [node for node in nodes]
        },
        edges=[(nodes.index(relation[0]), nodes.index(relation[1])) for relation in relations], 
        edge_attrs={
            "relation": [relation[2].name for relation in relations],
        },
        directed=True
    )    

    # Remove reflexive self-loops
    graph = graph.simplify(combine_edges="random")
    # Remove disjoint relations (are implicit)
    graph.delete_edges(relation_eq="DISJOINT")
    # Remove common ancestor redundancy
    # Remove one direction of containment
    graph.delete_edges(relation_eq="CONTAINS")
    # Remove transitive redundancies for single relations TODO include equivalence
    for relation in ("IS_CONTAINED",):
        graph.delete_edges(relation_eq=relation) 
        subgraph = nx.DiGraph([edge[:2] for edge in relations if edge[2].name == relation])
        subgraph = nx.transitive_reduction(subgraph) # Reduce edges which are redundant due to transitivity of the same relation
        graph.add_edges(
            [[nodes.index(node) for node in edge] for edge in subgraph.edges()], 
            {"relation": [relation for _ in subgraph.edges()]}
        )
    # Remove symmetric relations by treating graph as undirected
    graph = graph.as_undirected(combine_edges="random")

    # Convert back to edge list to be used elsewhere
    edges = [(nodes[e.source], nodes[e.target], e['relation']) for e in graph.es]
    for relation, count in count_relations(edges).items():
        print(relation, count)
    return edges