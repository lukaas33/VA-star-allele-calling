from itertools import chain, combinations

from modules.data import reference_get, pharmvar_get
from modules.parse import parse_multi_hgvs, combine_variants
import algebra as va
import difflib

def print_seq_diff(sequence1, sequence2):
    """ Output the difference between two sequences as insertions and deletions. """
    difference = difflib.ndiff(sequence1, sequence2)
    for i,s in enumerate(difference):
        if s[0]==' ': continue
        elif s[0]=='-':
            print(u'Delete "{}" on relative position {}'.format(s[-1],i+1))
        elif s[0]=='+':
            print(u'Insert "{}" on relative position {}'.format(s[-1],i+1)) 

def va_generate_subsets(variants):
    """Generate all subsets (superset) of a variant set"""
    superset = chain.from_iterable(combinations(variants, r) for r in range(1, len(variants)+1))
    return superset

def va_variants_independent(reference, variants):
    """Test if variants within an allele are independent
    
    Variants should be independent since they could be expressed as one variant otherwise. 
    """
    # TODO implement
    # QUESTION is it enough to test if each is disjoint with the rest?

def va_minimal_overlapping_allele(reference_sequence, lhs_allele, rhs_allele):
    """Minimize the size of overlapping alleles in order to make characterization more manageable.
    
    Uses top-down linear search.

    WARNING: worse case complexity is not good for bigger variants
    """
    # QUESTION is it valid to have both sequences be the same length
    # QUESTION binary search possible instead of linear search?
    # TODO change approach to work for more downstream overlapping variants (since area would now get large for no reason)
    #       use sliding and extending window?
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
    l_minset = set([tuple(minrep) for minrep in lhs_allele.atomics()])
    r_minset = set([tuple(minrep) for minrep in rhs_allele.atomics()])
    # Characterize the overlap
    shared = l_minset & r_minset # Set intersection
    l_only = l_minset - shared # Set minus
    r_only = r_minset - shared
    return l_only, shared, r_only

def test_coreallele_containment(corealleles, suballeles, reference_sequence, coreallele_name):
    """Test if the coreallele is contained in each suballele.
    
    By definition each suballele should contain the coreallele.
    If this is not the case there is an inconsistency.
    This is detected and then characterized.
    """
    # Get the variants contained in the sub- and core allele as HGVS notation
    # QUESTION why doesn't name always match position
    coreallele = corealleles[coreallele_name]
    coreallele_variants = [variant["position"] for variant in coreallele["variants"]]
    for suballele in suballeles:
        suballele_variants = [variant["position"] for variant in suballele["variants"]]
        # Convert HGVS notation to sequence notation
        s_coreallele_variants = parse_multi_hgvs(coreallele_variants, reference_sequence)
        s_suballele_variants = parse_multi_hgvs(suballele_variants, reference_sequence)
        # Find relation between core and suballeles
        try:
            relation = va.compare(reference_sequence, s_coreallele_variants, s_suballele_variants)
        except Exception as e:
            raise ValueError(f"Could not compare variants {coreallele_variants} with {suballele_variants}")
        # Expect containment or equivalence for core and suballele
        # QUESTION is an equivalence between a sub and core allele inconsistent?
        if relation not in (va.Relation.EQUIVALENT, va.Relation.IS_CONTAINED):
            # For unexpected relationships the overlap should be characterized
            print(f"{coreallele['alleleName']}: Unexpected relationship {relation} with suballele {suballele['alleleName']}:")
            only_core, shared, only_sub = va_characterize_overlap(reference_sequence, s_coreallele_variants, s_suballele_variants)
            print(f"\tMinimal representations found in the core but not in the suballele:")
            for min_repr in only_core:
                print(f"\t{va.variants.to_hgvs(min_repr, reference_sequence)}")

def main():
    # Get the reference sequence relevant for the (current) gene of interest
    reference_sequence = reference_get()
    # List genes as symbols in Pharmvar
    genes = pharmvar_get("genes/list") 
    # All information associated with the (current) gene of interest
    gene = pharmvar_get("genes/CYP2D6") 
    # Group suballeles by core alleles and index by the star-allele notation of the core allele
    # QUESTION what is a variant group in pharmvar
    # QUESTION is an empty allele always a subset of any allele
    #   TODO test this with empty alleles like *1
    corealleles = {allele["alleleName"]: allele for allele in gene["alleles"] if allele["alleleType"] == "Core"} 
    suballeles = {coreallele: [sub_allele for sub_allele in gene["alleles"] if sub_allele["coreAllele"] == coreallele] for coreallele in corealleles.keys()}

    # TEST 1: test if all corealleles are contained in their suballeles
    for coreallele_name, suballeles in suballeles.items():
        if (coreallele_name != "CYP2D6*57"): # Test
            continue
        test_coreallele_containment(corealleles, suballeles, reference_sequence, coreallele_name)

    # TODO check if HGVS name describes position field (not always the case)
    # TODO check if position is a valid HGVS string (not always the case)
    # TODO Check if names follow a logical format
    # TODO check relations between star-alleles
    # TODO check if variants within suballele are disjoint (not always the case) and solve this

if __name__ == "__main__":
    main()