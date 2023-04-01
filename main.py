import warnings
from modules.data import reference_get, pharmvar_get, parse_samples
from modules.graph import display_graph
from modules.compare import find_relations_all
from modules.relations import prune_relations, find_context
from modules.parse import extract_variants, to_supremal
from modules.data import cache_get, cache_set
from modules.calling import star_allele_calling, print_classification, sort_types
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
        
def test_coreallele_containment(corealleles, suballeles, relations_extended):
    """Test if the coreallele is contained in each suballele.
    
    By definition each suballele should contain the coreallele.
    There are cases where the suballele is equal to the coreallele.
    If neither is the case there is an inconsistency.
    This is detected and printed.
    """
    for core in corealleles.keys():
        for sub in suballeles[core].keys():
            for left, right, relation in relations_extended:
                if left == core and right == sub: # Can do this since the set is full and thus contains both directions
                    if relation == va.Relation.IS_CONTAINED:
                        break
                    if relation == va.Relation.EQUIVALENT:
                        break
                    warnings.warn(f"{core} is not contained in {sub}, but it is {relation.name}")
                    break

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
    corealleles |= {"CYP2D6*1": {"variants": [], "alleleName": "CYP2D6*1", "function": "normal function"}} # Add wild type TODO do this nicer
    suballeles = {coreallele: {sub_allele["alleleName"]: sub_allele for sub_allele in gene["alleles"] if sub_allele["coreAllele"] == coreallele} for coreallele in corealleles.keys()}
    # TODO check if more suballeles have an empty coreallele
    variants = {variant["hgvs"]: variant for allele in gene["alleles"] for variant in allele["variants"]}
    data = corealleles | variants
    for suballele in suballeles.values(): data = data | suballele
    functions = {a: d["function"] for a, d in data.items() if "function" in d}

    # TEST 0: test if naming is consistent
    # test_naming(corealleles, suballeles)

    # TEST 1: find the relation between all corealleles, suballeles and the contained variants
    supremal_extended = extract_variants(reference_sequence, corealleles, suballeles, cache_name="supremal_extended")
    relations_extended = find_relations_all(reference_sequence, supremal_extended, cache_name="relations_extended")	
    _, pruned_extended = prune_relations(relations_extended, cache_name="relations_pruned_extended")

    # TEST 1.1: validate the relations
    # validate_relations(relations_extended, variants, r"..\pharmvar-tools\data\pharmvar_5.2.19_CYP2D6_relations-nc.txt")
    # validate_relations(pruned_extended, variants, r"..\pharmvar-tools\data\pharmvar_5.2.19_CYP2D6_relations-nc-reduced.txt")

    # TEST 2: parse samples
    samples_source = parse_samples(reference_sequence) # TODO also check unphased # TODO cache
    try:
        # TODO solve: UserWarning: Could not parse sample NA18526A: unorderable variants
        supremal_samples = cache_get("supremal_samples")
    except:
        supremal_samples = {}
        for name, variants in samples_source.items():
            if variants is None:
                continue
            try:
                supremal_samples[name] = to_supremal(variants, reference_sequence)
            except ValueError as e:
                warnings.warn(f"Could not parse sample {name}: {e}")
        cache_set(supremal_samples, "supremal_samples")
    personal_variants = {variant: value for variant, value in supremal_samples.items() if sort_types(variant) == 3} 
    for p_variant in list(personal_variants.keys()): # TODO do this nicer
        for variant in supremal_extended.keys():
            if sort_types(variant) != 3:
                continue
            rel = va.relations.supremal_based.compare(reference_sequence, supremal_extended[variant], supremal_samples[p_variant])
            if rel == va.Relation.EQUIVALENT: # Remove personal variants which already exist
                del personal_variants[p_variant]
    samples = {sample: value for sample, value in supremal_samples.items() if sort_types(sample) == 4} 
    # Find all relations with samples
    # TODO simplify 
    # TODO run again to remove non-personal variants
    # TODO why are personal variants not included in the relations in this way?
    #   relations_samples = find_relations_all(reference_sequence, supremal_extended | personal_variants, supremal_samples, cache_name="relations_samples_extended") 
    relations_samples = find_relations_all(reference_sequence, supremal_extended, samples, cache_name="relations_samples_extended")
    relations_samples += find_relations_all(reference_sequence, samples, personal_variants, cache_name="relations_samples_personal")
    relations_samples += find_relations_all(reference_sequence, personal_variants, cache_name="relations_personal")

    # TODO verify sample relations

    # TEST 3: check if cores are contained in subs, if variants are contained in alleles and if personal variants are contained in samples
    test_coreallele_containment(corealleles, suballeles)
    exit()

    # TEST 4: determine star allele calling
    pruned_samples = prune_relations(pruned_extended + relations_samples, cache_name="relations_pruned_samples_extended")
    classifications = {sample[:-1]: {'A': None, 'B': None} for sample in sorted(samples_source.keys()) if sort_types(sample) == 4} 
    for sample in samples_source.keys():
        if sort_types(sample) != 4:
            continue
        sample_source, phasing = sample[:-1], sample[-1]
        classification = star_allele_calling(sample, *pruned_samples, functions)
        classifications[sample_source][phasing] = classification
    print_classification(classifications)

    # TEST 5: display some samples
    # TODO only show context of samples?
    sample_context = find_context(["NA21105A"], pruned_samples[1], as_edges=True)

    # VISUALIZE some context with information of interest
    context = pruned_extended
    pruned_nodes, pruned_edges = prune_relations(context + sample_context)
    display_graph(pruned_nodes, pruned_edges, data)

if __name__ == "__main__":
    main()