import warnings
from modules.data import reference_get, pharmvar_get
from modules.parse import parse_multi_hgvs
import algebra as va

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
        print(suballele["alleleName"])
        suballele_variants = [variant["hgvs"] for variant in suballele["variants"]]
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
            warnings.warn(f"{coreallele['alleleName']}: Unexpected relationship {relation} with suballele {suballele['alleleName']}")

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
        test_coreallele_containment(corealleles, suballeles, reference_sequence, coreallele_name)

    # TODO check if HGVS name describes position field (not always the case)
    # TODO check if position is a valid HGVS string (not always the case)
    # TODO Check if names follow a logical format
    # TODO check relations between star-alleles
    # TODO check if variants within suballele are disjoint (not always the case) and solve this

if __name__ == "__main__":
    main()