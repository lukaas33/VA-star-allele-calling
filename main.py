import warnings
from modules.data import reference_get, pharmvar_get, parse_samples
from modules.graph import display_graph
from modules.compare import find_relations_all
from modules.relations import prune_relations, find_context, redundant_reflexive
from modules.parse import extract_variants, to_supremal
from modules.data import cache_get, cache_set, api_get
from modules.calling import star_allele_calling_all, print_classification, sort_types, impact_position
from modules.other_sources import is_silent_mutalyzer, is_silent_entrez, find_id_hgvs
from modules.utils import validate_relations, validate_calling
from modules.assets.generate_images import *
import algebra as va

# TODO move checks to own file
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

def _test_pharmvar_annotation(variant, function, classification):
    """Test if a prediction of the function agrees with the pharmvar annotation."""
    if function is None: # No expectation
        # QUESTION: can make a prediction of function here?
        return True
    elif function == '': # No change
        if not classification['exon']: # No change because not in exon
            return True
        if not classification['non-synonymous']: # No change because synonymous
            return True
        warnings.warn(f"{variant} is annotated as 'no change' but may not be silent: {classification}")
    elif function == 'splice defect': # No protein level change expected
        if not classification['non-synonymous']: # No protein change because synonymous
            return True
        if classification['splicing']:
            return True
        warnings.warn(f"{variant} is annotated as 'splice defect' but may not be silent: {classification}")
    else: # Explicit change on protein level TODO split further
        if classification['exon']: # is in exon
            return True
        if classification['non-synonymous']: # is non-synonymous
            return True
        warnings.warn(f"{variant} is annotated as '{function}' but may be silent: {classification}")
    return False

def test_variant_annotation_mutalyzer(variants, functions):
    """Check if PharmVar annotations are consistent with the equivalent variant descriptions."""
    for variant in variants:
        function = functions[variant]
        classification = is_silent_mutalyzer(variant)
        _test_pharmvar_annotation(variant, function, classification)

def test_variant_annotation_entrez(variants, ids, functions):
    """Check if PharmVar annotations are consistent with the Entrez annotations."""
    for variant in variants:
        # Find classification 
        classification = is_silent_entrez(variant, ids)
        # Compare against PharmVar annotation
        _test_pharmvar_annotation(variant, functions[variant], classification)

def test_get_id(variants, ids, reference):
    """Test if the getting of the id works"""
    # TODO further investigate cases None/None
    # TODO check for personal variants
    for variant in variants:
        found_id = find_id_hgvs(variant, reference)
        id = ids[variant] 
        if id == '': id = None
        if found_id != id:
            warnings.warn(f"{variant} has id {id} but {found_id} was found")

def test_variant_annotation_position(variants, supremals, functions):
    """Test if the position of the variant is consistent with the annotation."""
    for variant in variants:
        classification = impact_position(variant, supremals, None)
        if functions[variant] is None or \
                functions[variant] == '' or \
                functions[variant] == 'splice defect': # No expectation on exon or intron
            continue
        else: # Protein effect
            if classification:
                warnings.warn(f"{variant} is annotated as '{functions[variant]}' but is in intron")

def main():
    # Get the reference sequence relevant for the (current) gene of interest
    reference_name = "NC_000022.11"
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
    ids = {variant["hgvs"]: variant["rsId"] for allele in gene["alleles"] for variant in allele["variants"]}
    # test_functional_annotation(suballeles, functions)
    # test_core_annotation(corealleles, functions)
    # test_variant_annotation_mutalyzer(variants, functions)
    # test_get_id(variants, ids, reference_sequence) # TODO run this
    # test_variant_annotation_entrez(variants, ids, functions) # TODO run this
    # test_variant_annotation_position(variants, supremal_extended, functions)

    # parse samples
    samples_source = parse_samples(reference_sequence) 
    try:
        supremal_samples = cache_get("supremal_samples")
    except:
        supremal_samples = {}
        for sample, variants in samples_source.items():
            if variants == []: # TODO handle differently (as empty variant?)
                continue
            # Parse individual variants
            for hgvs, variant in variants: # Find supremal for individual variants
                if hgvs in supremal_samples:
                    continue
                supremal_samples[hgvs] = to_supremal([variant], reference_sequence)
            # Parse alleles (variants together)
            try:
                supremal_samples[sample] = to_supremal([v[1] for v in variants], reference_sequence) # Try to find supremal for sample
            except ValueError as e:
                if "unorderable variants" in str(e):
                    warnings.warn(f"Could not parse sample {sample} due to double/overlapping variants")
                    print(*variants)
                    # TODO how to handle for unphased
                else:
                    raise e
        cache_set(supremal_samples, "supremal_samples")

    # Filter out non-personal variants (are already in dataset and will appear in relations)
    # TODO simplify or do earlier
    for sample, p_variants in samples_source.items():
        for p_variant, _ in p_variants:
            for variant in supremal_extended.keys():
                if sort_types(variant) != 3:
                    continue
                if p_variant in supremal_samples:
                    continue
                rel = va.relations.supremal_based.compare(reference_sequence, supremal_extended[variant], supremal_samples[p_variant])
                if rel != va.Relation.EQUIVALENT: # Remove personal variants which already exist in Pharmvar
                    continue
                del supremal_samples[p_variant] # Already exists in dataset
                for sample2 in samples_source.keys(): # Remove everywhere TODO needed?
                    for i, p_variant in enumerate(samples_source[sample2]):
                        if p_variant[0] != p_variant:
                            continue
                        del samples_source[sample2][i]
                break

    # Split into personal variants and samples
    personal_variants = {variant: value for variant, value in supremal_samples.items() if sort_types(variant) == 5} 
    samples = {sample: value for sample, value in supremal_samples.items() if sort_types(sample) == 4} 
    # Find all relations with samples
    relations_samples = find_relations_all(reference_sequence, supremal_extended, samples, cache_name="relations_samples_extended")
    relations_samples += find_relations_all(reference_sequence, samples, personal_variants, cache_name="relations_samples_personal")
    relations_samples += find_relations_all(reference_sequence, supremal_extended, personal_variants, cache_name="relations_personal_extended")
    relations_samples += find_relations_all(reference_sequence, personal_variants, cache_name="relations_personal")

    # TEST 4: check if relations are consistent with atomic variants
    # test_variant_containment(corealleles, suballeles, relations_extended)
    # test_personal_variant_containment(samples_source, relations_samples)
    # test_central_personal_variants(personal_variants.keys(), find_relations_all(reference_sequence, samples, personal_variants, cache_name="relations_samples_personal"))
    # test_central_personal_variants(personal_variants.keys(), relations_samples)

    # Determine star allele calling
    pruned_samples = prune_relations(pruned_extended + relations_samples, cache_name="relations_pruned_samples_extended")
    phased_samples = [sample for sample in samples_source.keys() if sample[-1] in ("A", "B")] 
    unphased_samples = [sample for sample in samples_source.keys() if sample[-1] not in ("A", "B")] # TODO use types for this (this is error prone)
    # calling_phased = star_allele_calling_all(phased_samples, *pruned_samples, functions, supremal_extended | supremal_samples)
    # print_classification(calling_phased, detail_level=0)
    calling_unphased = star_allele_calling_all(unphased_samples, *pruned_samples, functions, supremal_extended | supremal_samples, phased=False)
    print_classification(calling_unphased, detail_level=0)

    # TEST 5 validate star allele calling
    # validate_calling(calling_phased, r"data\bastard.txt")
    validate_calling(calling_unphased, r"data\bastard.txt")
    exit()

    # VISUALIZE some context with information of interest
    # TODO only show context of samples?
    sample_context = find_context(["HG00337"], pruned_samples[1], as_edges=True)
    context = pruned_extended
    pruned_nodes, pruned_edges = prune_relations(context + sample_context)
    display_graph(pruned_nodes, pruned_edges, data)

    # Generate images
    # nodes, edges, positions = image_reduction_equivalence(relations_extended)
    # display_graph(nodes, edges, data, positions=positions, default_layout="preset")

if __name__ == "__main__":
    main()