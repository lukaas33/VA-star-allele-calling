import algebra as va
from .parse import combine_variants
import difflib
from itertools import chain, combinations

def va_minimal_overlapping_allele(reference_sequence, lhs_allele, rhs_allele):
    """Minimize the size of overlapping alleles in order to make characterization more manageable.
    
    Uses top-down linear search.

    WARNING: worse case complexity is not good for bigger variants.
    """
    # QUESTION is it valid to have both sequences be the same length
    # QUESTION binary search possible instead of linear search?
    # TODO change approach to work for more downstream overlapping variants (since area would now get large for no reason)
    #       use sliding and extending window?
    # TODO change approach to work for multiple area's
    # Start extending alleles from left until an fragment with overlap is found
    if lhs_allele.start != rhs_allele.start or lhs_allele.end != rhs_allele.end:
        # TODO solve for lengths not being equal
        raise Exception("Lengths not the same")
    min_variant_l, min_variant_r = None, None
    for slice_r in range(1, lhs_allele.end - lhs_allele.start):
        subsequence_l, subsequence_r = lhs_allele.sequence[:slice_r], rhs_allele.sequence[:slice_r]
        sub_variant_l, sub_variant_r = va.Variant(lhs_allele.start, lhs_allele.start + slice_r, subsequence_l), va.Variant(lhs_allele.start, lhs_allele.start + slice_r, subsequence_r)
        relation = va.compare(reference_sequence, {sub_variant_l}, {sub_variant_r})
        if relation == va.Relation.OVERLAP:
            min_variant_l, min_variant_r =  sub_variant_l, sub_variant_r
            break
    # See if it can be reduced more by making it smaller while overlapping
    lhs_allele, rhs_allele = min_variant_l, min_variant_r
    for slice_l in range(0, lhs_allele.end - lhs_allele.start):
        subsequence_l, subsequence_r = lhs_allele.sequence[slice_l:], rhs_allele.sequence[slice_l:]
        sub_variant_l, sub_variant_r = va.Variant(lhs_allele.start + slice_l, lhs_allele.end, subsequence_l), va.Variant(lhs_allele.start + slice_l, lhs_allele.end, subsequence_r)
        relation = va.compare(reference_sequence, {sub_variant_l}, {sub_variant_r})
        if relation != va.Relation.OVERLAP: # Starts as overlapping, continue until it isn't anymore
            break
        min_variant_l, min_variant_r =  sub_variant_l, sub_variant_r
    lhs_allele, rhs_allele = min_variant_l, min_variant_r
    return lhs_allele, rhs_allele

def va_characterize_overlap(reference_sequence, lhs, rhs):
    """For a pair of variants characterize the elements of the overlap.
    
    This extends the variant algebra which detects but does not characterize the overlap.
    Overlap cannot always be explained by shared variant positions (simple set overlap).

    WARNING: not efficient, only works for very small overlapping areas
    """
    # QUESTION is the atomics of the allele equal to the (combination of) atomics of all variants
    # Combine variants into allele (supremal representation)
    lhs_allele = combine_variants(lhs, reference_sequence)
    rhs_allele = combine_variants(rhs, reference_sequence)
    # Reduce allele size to the minimal size with overlap to make the atomics more manageable
    lhs_allele, rhs_allele = va_minimal_overlapping_allele(reference_sequence, lhs_allele, rhs_allele)
    # Find the set of minimal representations of both alleles
    l_min_set = set([tuple(min_rep) for min_rep in lhs_allele.atomics()])
    r_min_set = set([tuple(min_rep) for min_rep in rhs_allele.atomics()])
    # Characterize the overlap
    shared = l_min_set & r_min_set # Set intersection
    l_only = l_min_set - shared # Set minus
    r_only = r_min_set - shared
    return l_only, shared, r_only

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