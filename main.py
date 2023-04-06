import warnings
from modules.data import reference_get, pharmvar_get, parse_samples
from modules.graph import display_graph
from modules.compare import find_relations_all
from modules.relations import prune_relations, find_context
from modules.parse import extract_variants, to_supremal
from modules.data import cache_get, cache_set, api_get
from modules.calling import star_allele_calling, print_classification, sort_types, classify_region
from modules.utils import validate_relations, validate_calling
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

def test_coreallele_containment(suballeles, relations_extended):
    """Test if the coreallele is contained in each suballele.
    
    By definition each suballele should contain the coreallele.
    There are cases where the suballele is equal to the coreallele.
    If neither is the case there is an inconsistency.
    This is detected and printed.
    """
    for core in suballeles.keys():
        for sub in suballeles[core].keys():
            for left, right, relation in relations_extended: # TODO use different ds here
                if left == core and right == sub: # Can do this since the set is full and thus contains both directions
                    if relation in (va.Relation.IS_CONTAINED, va.Relation.EQUIVALENT):
                        break
                    warnings.warn(f"{core} is not contained in {sub}, but it is {relation.name}")
                    break
            else:
                warnings.warn(f"Could not find relation for {core} and {sub}")

def test_variant_containment(corealleles, suballeles, relations_extended):
    """Test if the variants are contained in the coreallele or suballele.
    
    If not there are some interactions happening and this should be investigated.
    Variants that are not simply contained should not be viewed as atomic.
    """
    for coreallele, core_values in corealleles.items():
        for allele, values in list(suballeles[coreallele].items()) + [(coreallele, core_values)]:
            variants = [variant["hgvs"] for variant in values["variants"]]
            for variant in variants:
                for left, right, relation in relations_extended: # TODO use different ds here
                    if left == variant and right == allele: # Can do this since the set is full and thus contains both directions
                        if relation in (va.Relation.IS_CONTAINED, va.Relation.EQUIVALENT):
                            break
                        warnings.warn(f"{variant} is not contained in {allele}, but instead {relation.name}")
                        break
                else:
                    warnings.warn(f"Could not find relation for {variant} and {allele}")

def test_personal_variant_containment(samples_source, relations_samples):
    """Test if the variants found in a sample are contained in the entire allele
    
    Similar to the test_variant_containment but for the personal variants.
    """
    for sample, variants in samples_source.items():
        for hgvs in variants.keys():
            for left, right, relation in relations_samples: # TODO use different ds here
                if left == sample and right == hgvs:
                    if relation in (va.Relation.CONTAINS, va.Relation.EQUIVALENT):
                        break
                    warnings.warn(f"{hgvs} is not contained in {sample}, but instead {relation.name}")
                    break
            else:
                warnings.warn(f"Could not find relation for {sample} and {hgvs}")

def test_functional_annotation(suballeles, functions):
    """ Test if the functional annotation is consistent between core and suballele."""
    for core in suballeles.keys():
        for sub in suballeles[core].keys():
            if functions[core] != functions[sub]:
                warnings.warn(f"Function of {core} and {sub} is not consistent: {functions[core]} and {functions[sub]}")

def test_central_personal_variants(personal_variants, relations):
    """Test if personal variants are central, have a lot of relations.
    
    When personal variants have many relations to samples it would be expected that they are described in the literature.
    """
    for personal_variant in personal_variants:
        count = 0   
        for left, right, relation in relations: # TODO use different ds here
            if relation == va.Relation.DISJOINT:
                continue
            if left != personal_variant: # Data contains two directions
                continue
            if sort_types(right) != 4: 
                # warnings.warn(f"{left} has relation {relation.name} with {right} which is not a sample")
                continue
            if relation != va.Relation.IS_CONTAINED:
                continue
            count += 1
        if count > 1:
            warnings.warn(f"{personal_variant} has {count} relations while 1 is expected")

def test_core_annotation(corealleles, functions):
    """Test if the coreallele variants have an impact annotation."""
    for core in corealleles.keys():
        for variant in corealleles[core]["variants"]:
            variant = variant["hgvs"]
            if variant not in functions.keys():
                warnings.warn(f"{variant} is in {core} but has no function annotation")
            elif functions[variant] == "":
                warnings.warn(f"{variant} is in {core} but has 'no change' function annotation")

def test_variant_annotation(functions):
    """Find the correct annotation for the variants with no function annotation.
    
    Tests if variants are in an intron.
    If one representation of the variant is in an exon this is taken as the annotation (worst case).
    """
    # TODO store these annotations
    for variant, function in functions.items():
        if function is not None:
            continue
        # Find equivalent representations of the variant
        data = api_get(f"https://mutalyzer.nl/api/normalize/{variant}")
        if "equivalent_descriptions" not in data.keys():
            warnings.warn(f"Could not find equivalent descriptions for {variant}")
            continue
        data = data["equivalent_descriptions"]
        if all([t for t in data.keys() if t in {'c', 'n'}]):
                    warnings.warn(f"Unhandled types: {set(data.keys()) }")
                    continue
        if 'c' in data.keys():
            for repr, _ in data['c']:
                # print('\t', repr, classify_region(repr))
                if classify_region(repr) == 'exon':
                    print(variant, 'exon')
                    break
            else: # Not exon, do second check
                if 'n' in data.keys():
                    for repr in data['n']:
                        # print('\t', repr, classify_region(repr))
                        if classify_region(repr) == 'exon':
                            print(variant, 'exon')
                            break
                    else: # Not exon
                        print(variant, 'intron')
            

def main():
    # Get the reference sequence relevant for the (current) gene of interest
    reference_sequence = reference_get()
    # List genes as symbols in Pharmvar
    genes = pharmvar_get("genes/list") 
    # All information associated with the (current) gene of interest
    gene = pharmvar_get("genes/CYP2D6") 
    # Group suballeles by core alleles and index by the star-allele notation of the core allele
    # QUESTION should variants within a suballele be disjoint?
    # TODO use datastructure with less redundant information and for consistency
    # TODO add wild type relations (will be disjoint with all other variants)
    corealleles = {allele["alleleName"]: allele for allele in gene["alleles"] if allele["alleleType"] == "Core"} 
    corealleles |= {"CYP2D6*1": {"variants": [], "alleleName": "CYP2D6*1", "function": "normal function"}} # Add wild type TODO do this nicer
    suballeles = {coreallele: {sub_allele["alleleName"]: sub_allele for sub_allele in gene["alleles"] if sub_allele["coreAllele"] == coreallele} for coreallele in corealleles.keys()}
    variants = {variant["hgvs"]: variant for allele in gene["alleles"] for variant in allele["variants"]}
    # Save information of alleles and variants
    data = corealleles | variants
    for suballele in suballeles.values(): data = data | suballele
    # Get functional annotation of alleles and impact of variants
    functions = {a: d["function"] for a, d in data.items() if "function" in d}
    functions |= {a: d["impact"] for a, d in data.items() if "impact" in d}

    # TEST 1: test if naming is consistent
    # test_naming(corealleles, suballeles)

    # Find the relation between all corealleles, suballeles and the contained variants
    supremal_extended = extract_variants(reference_sequence, corealleles, suballeles, cache_name="supremal_extended")
    relations_extended = find_relations_all(reference_sequence, supremal_extended, cache_name="relations_extended")	
    _, pruned_extended = prune_relations(relations_extended, cache_name="relations_pruned_extended")

    # TEST 2: validate the relations
    # validate_relations(relations_extended, variants, r"..\pharmvar-tools\data\pharmvar_5.2.19_CYP2D6_relations-nc.txt")
    # validate_relations(pruned_extended, variants, r"..\pharmvar-tools\data\pharmvar_5.2.19_CYP2D6_relations-nc-reduced.txt")

    # TEST 3: check if the functional annotations are consistent
    # test_functional_annotation(suballeles, functions)
    # test_core_annotation(corealleles, functions)
    test_variant_annotation(functions)
    exit()

    # parse samples
    samples_source = parse_samples(reference_sequence) # TODO also check unphased 
    try:
        supremal_samples = cache_get("supremal_samples")
    except:
        supremal_samples = {}
        for sample, variants in samples_source.items():
            if variants == {}: # TODO handle differently
                continue
            for hgvs, variant in variants.items(): # Find supremal for individual variants
                supremal_samples[hgvs] = to_supremal([variant], reference_sequence)
            try:
                supremal_samples[sample] = to_supremal(list(variants.values()), reference_sequence) # Try to find supremal for sample
            except ValueError as e:
                warnings.warn(f"Could not parse sample {sample}: {e}")
        cache_set(supremal_samples, "supremal_samples")

    # Filter out non-personal variants (are already in dataset)
    for sample, p_variants in samples_source.items():
        for p_variant in list(p_variants.keys()):
            if sort_types(p_variant) != 5:
                continue
            for variant in supremal_extended.keys():
                if sort_types(variant) != 3:
                    continue
                rel = va.relations.supremal_based.compare(reference_sequence, supremal_extended[variant], supremal_samples[p_variant])
                if rel == va.Relation.EQUIVALENT: # Remove personal variants which already exist
                    del supremal_samples[p_variant] # Already exists in dataset
                    for sample2 in samples_source.keys(): # Remove everywhere
                        if p_variant in samples_source[sample2].keys():
                            del samples_source[sample2][p_variant]
                    break

    # Split into personal variants and samples
    personal_variants = {variant: value for variant, value in supremal_samples.items() if sort_types(variant) == 5} 
    samples = {sample: value for sample, value in supremal_samples.items() if sort_types(sample) == 4} 
    # Find all relations with samples
    # TODO simplify 
    relations_samples = find_relations_all(reference_sequence, supremal_extended, samples, cache_name="relations_samples_extended")
    relations_samples += find_relations_all(reference_sequence, supremal_extended, personal_variants, cache_name="relations_personal_extended")
    relations_samples += find_relations_all(reference_sequence, samples, personal_variants, cache_name="relations_samples_personal")
    relations_samples += find_relations_all(reference_sequence, personal_variants, cache_name="relations_personal")

    # TEST 4: check if relations are consistent with atomic variants
    # test_variant_containment(corealleles, suballeles, relations_extended)
    # test_personal_variant_containment(samples_source, relations_samples)
    # test_central_personal_variants(personal_variants.keys(), find_relations_all(reference_sequence, samples, personal_variants, cache_name="relations_samples_personal"))
    # test_central_personal_variants(personal_variants.keys(), relations_samples)

    # Determine star allele calling
    # TODO simplify
    pruned_samples = prune_relations(pruned_extended + relations_samples, cache_name="relations_pruned_samples_extended")
    classifications = {sample[:-1]: {'A': None, 'B': None} for sample in sorted(samples_source.keys()) if sort_types(sample) == 4} 
    for sample in samples_source.keys():
        if sort_types(sample) != 4:
            continue
        sample_source, phasing = sample[:-1], sample[-1]
        classification = star_allele_calling(sample, *pruned_samples, functions)
        classifications[sample_source][phasing] = classification
    # print_classification(classifications, detail_level=0)
    exit()

    # TEST 5 validate star allele calling
    validate_calling(classifications, r"data\bastard.txt")

    # display some samples
    # TODO only show context of samples?
    sample_context = find_context([], pruned_samples[1], as_edges=True)

    # VISUALIZE some context with information of interest
    context = pruned_extended
    pruned_nodes, pruned_edges = prune_relations(context + sample_context)
    # pruned_nodes, pruned_edges = pruned_samples
    display_graph(pruned_nodes, pruned_edges, data)

if __name__ == "__main__":
    main()