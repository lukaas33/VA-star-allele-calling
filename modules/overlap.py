import algebra as va

raise DeprecationWarning("This module is not useful for characterizing overlap.")

def va_minimal_overlapping_allele(reference_sequence, lhs_allele, rhs_allele):
    """Minimize the size of overlapping alleles in order to make characterization more manageable.
    
    Uses top-down linear search.

    WARNING: only finds one area, if there are more the leftmost area is found.
    WARNING: complexity is limiting for bigger variants.
    """
    # QUESTION is it valid to only have both sequences be the same length
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
        if relation != va.Relation.OVERLAP: # Starts as overlapping, continue until it isn't any more
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

def combine_variants(variants, reference_sequence):
    """ Combine multiple variants into one larger variant (allele) in the supremal representation
    
    WARNING: this is redundant for comparing since the compare function already patches the variants
    """
    # TODO replace this function with calls to the va library (see compare method for the proper calls)
    # Apply variants to reference_sequence to get the observed sequence
    min_start = float('inf')
    max_end = -float('inf')
    observed_sequence = va.variants.patch(reference_sequence, variants)
    for operation in variants:
        if operation.start < min_start:
            min_start = operation.start
        elif operation.end > max_end:
            max_end = operation.end
    # Find difference and output
    allele = va.Variant(min_start, max_end, observed_sequence[min_start:max_end])
    return allele