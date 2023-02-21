import algebra as va
from .data import cache_get, cache_set
from .parse import parse_multi_hgvs

def find_relations(corealleles, reference_sequence):
    """Find the relation between all corealleles.

    Relations are cached since they take a long time to generate.

    Returns list of edges.
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
            # TODO change loop to remove all none values
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
    Transitive relations can be reduced to a tree structure (equivalence, containment).
    Disjoint relation can be left out.

    returns list of edges.
    """
    transitive = (va.Relation.EQUIVALENT, va.Relation.CONTAINS, va.Relation.IS_CONTAINED)
    symmetric = (va.Relation.EQUIVALENT, va.Relation.OVERLAP, va.Relation.DISJOINT)
    edges = []
    check_symmetric = set() 
    check_transitivity = {r: [] for r in transitive}
    for node, other, relation in relations:
        if relation is None:
            raise ValueError("Relation data is incomplete")
        # Skip trivial self equivalence (reflexivity)
        if node == other: 
            continue
        # Don't represent disjointedness explicitly, 
        # it can be seen as 'no relation'
        if relation == va.Relation.DISJOINT: 
            continue
        # Don't represent symmetric relations twice
        if relation in symmetric:
            pair = (node, other)
            inv_pair = (other, node)
            if pair in check_symmetric or inv_pair in check_symmetric: 
                continue
            check_symmetric.add(pair)
            check_symmetric.add(inv_pair)
        # Store transitive relations to prune and add later
        if relation in transitive:
            check_transitivity[relation].append((node, other, relation))
            continue
        edges.append((node, other, relation))
    # Reduce transitive graph to tree and add
    pass
    return edges