import algebra as va
import difflib
from itertools import chain, combinations

def print_seq_diff(sequence1, sequence2, start=1):
    """ Output the difference between two aligned sequences as insertions and deletions. """
    difference = difflib.ndiff(sequence1, sequence2)
    for i,s in enumerate(difference):
        if s[0]==' ': continue
        elif s[0]=='-':
            print(u'Delete "{}" on relative position {}'.format(s[-1],start+i))
        elif s[0]=='+':
            print(u'Insert "{}" on relative position {}'.format(s[-1],start+i)) 

def va_generate_subsets(variants):
    """Generate all subsets (superset) of a variant set"""
    superset = chain.from_iterable(combinations(variants, r) for r in range(1, len(variants)+1))
    return superset

def count_relations(relations):
    counts = {relationType.name: 0 for relationType in va.Relation}
    for _, _, relation in relations:
        counts[relation.name] += 1
    return counts

def count_arity(nodes, relations):
    arity = {node: {relationType.name: 0 for relationType in va.Relation} for node in nodes}
    for l_allele, r_allele, relation in relations:
        arity[l_allele][relation.name] += 1
        if relation in (va.Relation.IS_CONTAINED, va.Relation.CONTAINS): # Directional
            continue
        arity[r_allele][relation.name] += 1
    return arity
