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
            print(u'Delete "{}" from position {}'.format(s[-1],i+1))
        elif s[0]=='+':
            print(u'Add "{}" to position {}'.format(s[-1],i+1)) 

def va_generate_subsets(variants):
    """Generate all subsets (superset) of a variant set"""
    superset = chain.from_iterable(combinations(variants, r) for r in range(1, len(variants)+1))
    return superset

def va_variants_independent(reference, variants):
    """Test if variants within an allele are independent
    
    Variants should be independent since they could be expressed as one variant otherwise. 
    """
    # TODO is it enough to test if each is disjoint with the rest?
    # if len(variants) == 1:
    #     return 0
    # variants = list(variants)
    # print(variants)
    # for pivot in range(1, len(variants)):
    #     lhs = variants[:pivot]
    #     rhs = variants[pivot:]
    #     relation = va.compare(reference, lhs, rhs)
    #     if relation != va.Relation.DISJOINT:
    #         print(f"{relation}: {lhs}; {rhs}")

def va_characterize_overlap(reference_sequence, lhs, rhs):
    """For a pair of variants characterize the elements of the overlap.
    
    This extends the variant algebra which detects but does not characterize the overlap.
    Overlap cannot always be explained by shared variant positions (simple set overlap).
    """
    # TODO integrate in va?
    # TODO make more efficient, runtime/memory usage is too high
    # QUESTION is the atomics of the allele equal to the (combination of) atomics of all variants
    # Combine variants into allele
    lhs_allele = combine_variants(lhs, reference_sequence)
    rhs_allele = combine_variants(rhs, reference_sequence)
    # Find the set of minimal representations of both alleles
    l_minset = set([tuple(minrep) for minrep in lhs_allele.atomics()])
    r_minset = set([tuple(minrep) for minrep in rhs_allele.atomics()])
    # Characterize the overlap
    shared = l_minset & r_minset # Set intersection
    l_only = l_minset - shared # Set minus
    r_only = r_minset - shared
    return l_only, shared, r_only

    # Find largest subset of lhs that is disjoint with the rhs 
    # These are the elements that cause the relation to be overlap instead of containment
    # Alternatively find the largest subset that is contained
    # for l_subset in va_generate_subsets(lhs):
    #     relation = va.compare(reference, l_subset, rhs)
    #     if relation == va.Relation.DISJOINT:
    #         print(l_subset)
    #         # return l_subset

    # Find the smallest subset of lts that is disjoint with rhs,
    # This explains why the lhs is not contained in the rhs
    # disjoint = []
    # for l_variant in lhs: 
    #     for l_repr in l_variant.atomics():
    #         relation = va.compare(reference, l_repr, rhs)
    #         if relation == va.Relation.DISJOINT: # Variant found which is not in the right hand side
    #             disjoint.append(l_variant)
    # if len(disjoint) == 0:
    #     raise Exception("No disjoint items found that characterizes the overlap")
    # corrected_lhs = [variant for variant in lhs if variant not in disjoint]
    # assert va.compare(reference, corrected_lhs, rhs) == va.Relation.IS_CONTAINED
    # return disjoint

    # va_overlap = []
    # for side in (lhs, rhs):
    #     other = lhs if side is rhs else rhs
    #     overlap = list(side)
    #     i = 0
    #     while i < len(overlap):
    #         removed = overlap.pop(i) # Remove variant to test overlap
    #         relation = va.compare(reference, overlap, other)
    #         if relation != va.Relation.OVERLAP: # There is no overlap anymore
    #             overlap.insert(i, removed) # Variant should be in the overlap
    #             i += 1
    #     va_overlap.append(overlap)
    # # return va_overlap
    # contained = list()
    # overlap = list(lhs)
    # i = 0
    # # while i < len(overlap):
    #     moved = overlap[i] # Try to move to contained
    #     contained.append(moved)
    #     relation = va.compare(reference, contained, rhs)
    #     if relation != va.Relation.IS_CONTAINED: # Set is no longer contained
    #         contained.pop() # Don't move element to contained
    #         i += 1 # Go to next element
    #     else: 
    #         overlap.pop(i) # Move element from overlap to contained 
    # assert va.compare(reference, contained, rhs) == va.Relation.IS_CONTAINED
    # assert va.compare(reference, overlap, rhs) == va.Relation.OVERLAP
    # print(contained)
    # print(overlap)
    # return overlap

def test_coreallele_containment(corealleles, suballeles, reference_sequence, coreallele_name):
    """Test if the coreallele is contained in each suballele."""
    # All information associated with the allele of interest
    coreallele = corealleles[coreallele_name]
    # Get the variants contained in the sub- and core allele as HGVS notation
    # QUESTION why doesn't name match position
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
        # Expect containment or equivalence
        if relation not in (va.Relation.EQUIVALENT, va.Relation.IS_CONTAINED):
            # For unexpected relationships the overlap should be characterized
            # QUESTION is an equivalence between a sub and core allele inconsistent?
            print(f"{coreallele['alleleName']}: Unexpected relationship {relation} with suballele {suballele['alleleName']}:")
            # only_core, shared, only_sub = va_characterize_overlap(reference_sequence, s_coreallele_variants, s_suballele_variants)
            # print(f"\t{only_core} is found in the core but not in the suballele")

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