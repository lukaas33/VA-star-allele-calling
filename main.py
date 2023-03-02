import warnings
from modules.data import reference_get, pharmvar_get
from modules.parse import parse_multi_hgvs
from modules.graph import display_graph
from modules.relations import find_relations_all
import algebra as va

def test_naming(corealleles, suballeles):
    """Naming test, useful for checking if the data is complete."""
    numbers = []
    for coreallele_name in corealleles.keys():
        for suballele in suballeles[coreallele_name]:
            if coreallele_name not in suballele["alleleName"]:
                warnings.warn(f"{coreallele_name}: naming of {suballele['alleleName']} inconsistent.")
        number = int(coreallele_name.split("*")[1])
        numbers.append(number)
    if len(numbers) != len(set(numbers)):
        warnings.warn("There are double core alleles in the data")
    all_numbers = set(range(1, max(numbers)+1))
    if all_numbers != set(numbers):
        # QUESTION why are some numbers skipped?
        # TODO detect proper skips
        # *1 is wild type
        # *5 is full gene deletion 
        # *13, etc. There are hybrid genes which 
        # *16, etc.
        warnings.warn(f"Not all numbers present as coreallele: {sorted(list(all_numbers - set(numbers)))}")
        
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
        # Convert HGVS notation to sequence notation
        suballele_variants = [variant["hgvs"] for variant in suballele["variants"]]
        s_coreallele_variants = parse_multi_hgvs(coreallele_variants, reference_sequence, coreallele_name)
        s_suballele_variants = parse_multi_hgvs(suballele_variants, reference_sequence, suballele_name)
        # Find relation between core and suballeles
        try:
            relation = va.compare(reference_sequence, s_coreallele_variants, s_suballele_variants)
        except Exception as e:
            if "unorderable variant" in str(e): 
                # This happens when variants overlap, 
                # in this case the observed sequence cannot be derived and the variants are not interpretable
                warnings.warn(f"{coreallele_name}: Some variants overlap. Cannot interpret, ignoring current allele.")
                # Skip uninterpretable variants 
                continue
            raise ValueError(f"{coreallele_name}: Could not compare variants with {suballele_name} ({coreallele_variants} and {suballele_variants})")
        # Expect containment or equivalence for core and suballele
        # QUESTION is an equivalence between a sub and core allele inconsistent?
        if relation not in (va.Relation.EQUIVALENT, va.Relation.IS_CONTAINED):
            # Not characterizing overlap since this is hard
            warnings.warn(f"{coreallele_name}: Unexpected relationship {relation} with suballele {suballele_name}")

def main():
    # Get the reference sequence relevant for the (current) gene of interest
    reference_sequence = reference_get()
    # List genes as symbols in Pharmvar
    genes = pharmvar_get("genes/list") 
    # All information associated with the (current) gene of interest
    gene = pharmvar_get("genes/CYP2D6") 
    # Group suballeles by core alleles and index by the star-allele notation of the core allele
    # QUESTION is an empty allele always a subset of any allele
    # QUESTION should variants within a suballele be disjoint?
    # TODO use datastructure with less redundant information for consistency
    corealleles = {allele["alleleName"]: allele for allele in gene["alleles"] if allele["alleleType"] == "Core"} 
    suballeles = {coreallele: [sub_allele for sub_allele in gene["alleles"] if sub_allele["coreAllele"] == coreallele] for coreallele in corealleles.keys()}

    # TEST 0: test if naming is consistent
    # test_naming(corealleles, suballeles)

    # TEST 1: test if all corealleles are contained in their suballeles
    # for coreallele_name, suballeles in suballeles.items():
        # test_coreallele_containment(corealleles, suballeles, reference_sequence, coreallele_name)

    # TEST 2: find the relation between all corealleles and analyze
    relations = find_relations_all(corealleles, reference_sequence)
    display_graph(relations, corealleles | suballeles)

if __name__ == "__main__":
    main()