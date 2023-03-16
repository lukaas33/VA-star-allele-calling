import warnings
from modules.data import reference_get, pharmvar_get, parse_samples
from modules.graph import display_graph
from modules.compare import find_relations_all
from modules.relations import prune_relations, find_context
from modules.parse import extract_variants, to_supremal
from modules.data import cache_get, cache_set
from modules.calling import star_allele_calling, print_classification
from modules.utils import validate_relations
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
    # TODO do this on main dataset
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
    # TODO use datastructure with less redundant information and for consistency
    corealleles = {allele["alleleName"]: allele for allele in gene["alleles"] if allele["alleleType"] == "Core"} 
    corealleles |= {"CYP2D6*1": {"variants": [], "alleleName": "CYP2D6*1"}} # Add wild type TODO do this nicer
    suballeles = {coreallele: {sub_allele["alleleName"]: sub_allele for sub_allele in gene["alleles"] if sub_allele["coreAllele"] == coreallele} for coreallele in corealleles.keys()}
    # TODO check if more suballeles have an empty coreallele
    variants = {variant["hgvs"]: variant for allele in gene["alleles"] for variant in allele["variants"]}
    data = corealleles | variants
    for suballele in suballeles.values(): data = data | suballele

    # TEST 0: test if naming is consistent
    # test_naming(corealleles, suballeles)

    # TEST 1: test if all corealleles are contained in their suballeles
    # for coreallele_name, suballeles in suballeles.items():
        # test_coreallele_containment(corealleles, suballeles, reference_sequence, coreallele_name)

    # TEST 2: find the relation between all corealleles, suballeles and the contained variants
    supremal_extended = extract_variants(reference_sequence, corealleles, suballeles, cache_name="supremal_extended")
    relations_extended = find_relations_all(reference_sequence, supremal_extended, cache_name="relations_extended")	

    # TEST 2.1: validate the relations
    # validate_relations(relations_extended, variants, r"..\pharmvar-tools\data\pharmvar_5.2.19.1_CYP2D6_relations-nc.txt")

    # TEST 3: parse samples
    try:
        supremal_samples = cache_get("supremal_samples")
    except:
        samples = parse_samples() # TODO also check unphased
        supremal_samples = {}
        for name, variants in samples.items():
            try:
                supremal_samples[name] = to_supremal(variants, reference_sequence)
            except ValueError as e:
                warnings.warn(f"Could not parse sample {name}: {e}")
        cache_set(supremal_samples, "supremal_samples")
    relations_samples = find_relations_all(reference_sequence, supremal_extended, supremal_samples, cache_name="relations_samples_extended") 
    # TODO also find sample variant relations to sample? Or remove them 
    # TODO verify sample relations

    # TEST 4: determine star allele calling
    supremal_samples = {sample: value for sample, value in supremal_samples.items() if sample[:2] in ('HG', 'NA')} # Don't need sample variants for this TODO find nicer way to filter this
    classifications = {sample[:-1]: {'A': None, 'B': None} for sample in sorted(supremal_samples.keys())} 
    for sample in supremal_samples.keys():
        classification = star_allele_calling(sample, relations_samples)
        sample, phasing = sample[:-1], sample[-1]
        classifications[sample][phasing] = classification
    # print_classification(classifications)

    # TEST 5: display all samples
    # TODO display subset (pruning takes a long time)
    sample_context = find_context(["HG00337A", "HG00337B"], relations_samples, as_edges=True)

    # VISUALIZE
    _, pruned_extended = prune_relations(relations_extended, cache_name="relations_pruned_extended")
    pruned = prune_relations(pruned_extended + sample_context)
    display_graph(*pruned, data)

if __name__ == "__main__":
    main()