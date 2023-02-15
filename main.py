import warnings
from modules.data import reference_get, pharmvar_get
from modules.parse import parse_multi_hgvs
import algebra as va

def test_variant_overlapping(variants, allele_name, reference_sequence):
    """Test overlapping variants in an variant list.
    
    This is needed since overlapping variants cannot be ordered and thus not compared.
    Because they are not interpretable without more context no relation can be determined.

    Adjacent variants are not a problem although they are strictly not HGVS proof if variants are not independent.
    """
    # TODO replace function by listening to certain error
    # TODO document in thesis and report to pharmvar
    correct = True
    sorted_variants = sorted(list(variants), key=lambda v: v.start)
    for i in range(0, len(sorted_variants)-1):
        if sorted_variants[i].end == sorted_variants[i+1].start: # Adjacency
            continue # Ignore
        elif sorted_variants[i].end > sorted_variants[i+1].start: # Overlap
            warnings.warn(f"{allele_name}: Variants {va.variants.to_hgvs([sorted_variants[i]], reference='NC000022.11', sequence_prefix=True)} and {va.variants.to_hgvs([sorted_variants[i+1]], reference='NC000022.11', sequence_prefix=True)} overlap. Cannot interpret, ignoring current allele.")
            correct = False
    return correct

def test_coreallele_containment(corealleles, suballeles, reference_sequence, coreallele_name):
    """Test if the coreallele is contained in each suballele.
    
    By definition each suballele should contain the coreallele.
    If this is not the case there is an inconsistency.
    This is detected and then characterized.
    """
    # Get the variants contained in the sub- and core allele as HGVS notation
    # hgvs has been tested to be consistent with fasta files
    # QUESTION why doesn't name always match position
    coreallele = corealleles[coreallele_name]
    coreallele_variants = [variant["hgvs"] for variant in coreallele["variants"]]
    for suballele in suballeles:
        suballele_name = suballele["alleleName"]
        # Check if name logical
        if coreallele_name not in suballele_name:
            warnings.warn(f"{coreallele_name}: naming of {suballele_name} inconsistent.")
        # Convert HGVS notation to sequence notation
        suballele_variants = [variant["hgvs"] for variant in suballele["variants"]]
        s_coreallele_variants = parse_multi_hgvs(coreallele_variants, reference_sequence, coreallele_name)
        s_suballele_variants = parse_multi_hgvs(suballele_variants, reference_sequence, suballele_name)
        # Skip uninterpretable variants
        if not test_variant_overlapping(s_coreallele_variants, coreallele_name, reference_sequence) or \
                not test_variant_overlapping(s_suballele_variants, suballele_name, reference_sequence):
            continue 
        # Find relation between core and suballeles
        try:
            relation = va.compare(reference_sequence, s_coreallele_variants, s_suballele_variants)
        except Exception as e:
            raise ValueError(f"{coreallele_name}: Could not compare variants with {suballele_name} ({coreallele_variants} and {suballele_variants})")
        # Expect containment or equivalence for core and suballele
        # QUESTION is an equivalence between a sub and core allele inconsistent?
        if relation not in (va.Relation.EQUIVALENT, va.Relation.IS_CONTAINED):
            # Not characterizing overlap since this is hard
            warnings.warn(f"{coreallele_name}: Unexpected relationship {relation} with suballele {suballele_name}")

def find_relations(corealleles, reference_sequence):
    """Find the relation between all corealleles and generate a graph structure.
    
    The graph structure should contain the relations in a minimal way, 
    e.g. removing redundant symmetric or transitive relations.

    This graph structure should later be visualized.
    """
    # Find relation for each pair of corealleles directionally
    # TODO is it needed to do the inverse test for containment?
    coreallele_names = list(corealleles.keys())
    relations = {
        coreallele: {
            coreallele2: None 
            for coreallele2 in coreallele_names
        } 
        for coreallele in coreallele_names
    }
    for i, left_coreallele in enumerate(coreallele_names):
        for right_coreallele in coreallele_names[i+1:]:
            left_hgvs = [variant["hgvs"] for variant in corealleles[left_coreallele]["variants"]]
            left_variants = parse_multi_hgvs(left_hgvs, reference_sequence)
            right_hgvs = [variant["hgvs"] for variant in corealleles[right_coreallele]["variants"]]
            right_variants = parse_multi_hgvs(right_hgvs, reference_sequence)
            relation = va.compare(reference_sequence, left_variants, right_variants)
            relations[left_coreallele][right_coreallele] = relation
            print(left_coreallele, right_coreallele, relation)

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
    # QUESTION should variants within a suballele be disjoint?
    #   TODO test this with empty alleles like *1
    corealleles = {allele["alleleName"]: allele for allele in gene["alleles"] if allele["alleleType"] == "Core"} 
    suballeles = {coreallele: [sub_allele for sub_allele in gene["alleles"] if sub_allele["coreAllele"] == coreallele] for coreallele in corealleles.keys()}

    # TEST 1: test if all corealleles are contained in their suballeles
    for coreallele_name, suballeles in suballeles.items():
        test_coreallele_containment(corealleles, suballeles, reference_sequence, coreallele_name)

    # TEST 2: find the relation between all corealleles
    # find_relations(corealleles, reference_sequence)

if __name__ == "__main__":
    main()