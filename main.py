import warnings
import argparse
from modules.data import reference_get, pharmvar_get, parse_samples
from modules.graph import display_graph
from modules.compare import find_relations_all
from modules.relations import prune_relations, find_context, redundant_reflexive
from modules.parse import extract_variants, to_supremal
from modules.data import cache_get, cache_set, api_get
from modules.calling import star_allele_calling_all, find_type, Type, impact_position, relevance
from modules.other_sources import is_silent_mutalyzer, get_annotation_entrez, find_id_hgvs, get_personal_ids, get_personal_impacts
from modules.utils import validate_relations, validate_calling, make_samples_unphased, validate_alternative_calling
from modules.assets.generate_images import image_configs
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
            if find_type(right) != Type.SAMPLE: 
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
        classification = get_annotation_entrez(variant, ids)
        # Compare against PharmVar annotation
        _test_pharmvar_annotation(variant, functions[variant], classification)

def test_get_id(variants, ids, reference):
    """Test if the getting of the id works"""
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

def test_extended_simplified(samples, pruned_samples_simple, supremal_simple, pruned_samples_extended, supremal_extended, supremal_samples, functions, reference):
    for phase in ('het', 'hom', 'all', 'A', 'B'):
        # TODO add inferred phasing experiment
        sel_samples = [sample for sample in samples.keys() if sample.split('_')[1] == phase] 
        sel_calling_simple = star_allele_calling_all(sel_samples, *pruned_samples_simple, functions, supremal_simple | supremal_samples, reference, detail_level=1)
        sel_calling_extended = star_allele_calling_all(sel_samples, *pruned_samples_extended, functions, supremal_extended | supremal_samples, reference, detail_level=1)
        if sel_calling_simple != sel_calling_extended:
            warnings.warn(f"Calling with {phase} is not the same for simple and extended")
            for sel_sample in sel_samples:
                sample, phase = sel_sample.split('_')
                if sel_calling_simple[sample] != sel_calling_extended[sample]:
                    print(sel_sample, *sel_calling_simple[sample][phase], "vs", *sel_calling_extended[sample][phase])
        else:
            print(f"Calling with {phase} variants is the same for simple and extended relations")

def main(text, visual, example, select, interactive, phased, unphased, detail, download):
    # Get the reference sequence relevant for the (current) gene of interest
    reference_name = "NC_000022.11"
    reference_sequence = reference_get(reference_name)
    reference = {"name": reference_name, "sequence": reference_sequence}
    # List genes as symbols in Pharmvar
    genes = pharmvar_get("genes/list") 
    # All information associated with the (current) gene of interest
    gene = pharmvar_get("genes/CYP2D6") 
    # Group suballeles by core alleles and index by the star-allele notation of the core allele
    # QUESTION should variants within a suballele be disjoint?
    # TODO use datastructure with less redundant information
    corealleles = {allele["alleleName"]: allele for allele in gene["alleles"] if allele["alleleType"] == "Core"} 
    corealleles |= {"CYP2D6*1": {"variants": [], "alleleName": "CYP2D6*1", "function": "normal function"}} # Add wild type so that the suballeles are added
    suballeles = {coreallele: {sub_allele["alleleName"]: sub_allele for sub_allele in gene["alleles"] if sub_allele["coreAllele"] == coreallele} for coreallele in corealleles.keys()}
    variants = {variant["hgvs"]: variant for allele in gene["alleles"] for variant in allele["variants"]}
    # Save information of alleles and variants
    data = corealleles | variants
    for suballele in suballeles.values(): data = data | suballele
    # Get functional annotation of alleles and impact of variants
    ids = {variant["hgvs"]: variant["rsId"] for allele in gene["alleles"] for variant in allele["variants"]}
    functions = {a: d["function"] for a, d in data.items() if "function" in d}
    functions |= {a: d["impact"] for a, d in data.items() if "impact" in d}

    # TEST 1: test if naming is consistent
    # test_naming(corealleles, suballeles)

    # Find the relation between all corealleles, suballeles and the contained variants
    supremal_extended = extract_variants(reference_sequence, corealleles, suballeles, cache_name="supremal_extended")
    relations_extended = find_relations_all(reference_sequence, supremal_extended, cache_name="relations_extended")	
    pruned_extended = prune_relations(relations_extended, cache_name="relations_pruned_extended")
    pruned_extended[0].add("CYP2D6*1") # Add since it won't be found in the relations
    supremal_simple = extract_variants(reference_sequence, corealleles, cache_name="supremal")
    relations_simple = find_relations_all(reference_sequence, supremal_simple, cache_name="relations_simple")
    pruned_simple = prune_relations(relations_simple, cache_name="relations_pruned_simple")
    pruned_simple[0].add("CYP2D6*1") # Add since it won't be found in the relations

    # TEST 2: validate the relations
    # validate_relations(relations_extended, variants, r"..\pharmvar-tools\data\pharmvar_5.2.19_CYP2D6_relations-nc.txt")
    # validate_relations(pruned_extended[1], variants, r"..\pharmvar-tools\data\pharmvar_5.2.19_CYP2D6_relations-nc-reduced.txt")

    # TEST 3: check if the functional annotations are consistent
    # test_functional_annotation(suballeles, functions)
    # test_core_annotation(corealleles, functions)
    # test_variant_annotation_mutalyzer(variants, functions)
    # test_get_id(variants, ids, reference_sequence) 
    # test_variant_annotation_entrez(variants, ids, functions) 
    # test_variant_annotation_position(variants, supremal_extended, functions)

    # parse samples
    samples_phased = parse_samples("data/samples", reference_sequence, phased=True) 
    samples_unphased = parse_samples("data/samples_unphased", reference_sequence) 
    # TODO move to function
    try:
        supremal_samples = cache_get("supremal_samples")
    except:
        supremal_samples = {}
        for sample, variants in (samples_phased | samples_unphased).items():
            if variants == {}: # No supremal for empty, will be disjoint with everything, can be ignored
                supremal_samples[sample] = None
                continue
            # Parse alleles (variants together)
            try:
                supremal_samples[sample] = to_supremal(list(variants.values()), reference_sequence) # Try to find supremal for sample
            except ValueError as e:
                # TODO can assume that overlapping variants are in different phases?
                # TODO use different placement for the overlapping variants?
                if "unorderable variants" in str(e):
                    warnings.warn(f"Could not parse sample {sample} due to double/overlapping variants")
                    print(*variants.keys())
                else:
                    raise e
            # Parse individual variants
            for hgvs, variant in list(variants.items()): # Find supremal for individual variants
                if hgvs in supremal_samples:
                    continue
                supremal_v = to_supremal([variant], reference_sequence)
                # Filter out variants that are in the database (not personal)
                # Relations with these variants will show up later
                for v in supremal_extended:
                    if find_type(v) != Type.VAR:
                        continue
                    rel = va.relations.supremal_based.compare(reference_sequence, supremal_v, supremal_extended[v])
                    if rel == va.Relation.EQUIVALENT: # Not personal
                        del variants[hgvs] # Can delete since already parsed as allele
                        break
                else: # Personal variant is saved individually
                    supremal_samples[hgvs] = supremal_v
        cache_set(supremal_samples, "supremal_samples")

    # Split into personal variants and samples
    personal_variants = {variant: value for variant, value in supremal_samples.items() if find_type(variant) == Type.P_VAR} 
    samples = {sample: value for sample, value in supremal_samples.items() if find_type(sample) == Type.SAMPLE} 

    # TEST 4: check if more information can be found about personal variants.
    ids |= get_personal_ids(personal_variants, reference, cache_name="ids_personal")
    functions |= get_personal_impacts(personal_variants, ids, reference, cache_name="impacts_personal")
 
    # Find all relations with samples
    relations_samples_extended = find_relations_all(reference_sequence, supremal_extended, samples, cache_name="relations_samples_extended")
    relations_samples_extended += find_relations_all(reference_sequence, samples, personal_variants, cache_name="relations_samples_personal")
    relations_samples_extended += find_relations_all(reference_sequence, supremal_extended, personal_variants, cache_name="relations_personal_extended")
    relations_samples_extended += find_relations_all(reference_sequence, personal_variants, cache_name="relations_personal")
    # Simplified
    relations_samples_simple = find_relations_all(reference_sequence, supremal_simple, samples, cache_name="relations_samples_simple")
    relations_samples_simple += find_relations_all(reference_sequence, samples, personal_variants, cache_name="relations_samples_personal")
    relations_samples_simple += find_relations_all(reference_sequence, supremal_simple, personal_variants, cache_name="relations_personal_simple")
    relations_samples_simple += find_relations_all(reference_sequence, personal_variants, cache_name="relations_personal")

    # TEST 5: check if relations are consistent with atomic variants
    # test_variant_containment(corealleles, suballeles, relations_extended)
    # test_personal_variant_containment(samples_source, relations_samples)
    # test_central_personal_variants(personal_variants.keys(), find_relations_all(reference_sequence, samples, personal_variants, cache_name="relations_samples_personal"))
    # test_central_personal_variants(personal_variants.keys(), relations_samples)

    # Simplify sample relations
    pruned_samples_extended = prune_relations(relations_samples_extended + relations_extended, cache_name="relations_pruned_samples_extended")
    pruned_samples_simple = prune_relations(relations_samples_simple + relations_simple, cache_name="relations_pruned_samples_simple")

    # TEST 6: test if star allele based on corealleles is the same as calling with suballeles
    # test_extended_simplified(samples, pruned_samples_simple, supremal_simple, pruned_samples_extended, supremal_extended, supremal_samples, functions, reference)

    # EXPERIMENT 1: Determine star allele calling for phased samples
    if unphased and phased: 
        raise ValueError("Either phased or unphased should be True")
    if phased:
        # TODO change used by detail?
        calling = star_allele_calling_all(samples_phased, *pruned_samples_extended, functions, supremal_extended | supremal_samples, reference, detail_level=detail, reorder=text is not None)
    
    # EXPERIMENT 2: Determine star allele calling for unphased samples
    # EXPERIMENT 2.1: use all variants in single allele
    # sel_samples = [sample for sample in samples_unphased.keys() if sample.split('_')[1] == 'all'] 
    # calling_unphased = star_allele_calling_all(sel_samples, *pruned_samples_simple, functions, supremal_extended | supremal_samples, reference, detail_level=1)
    # for sample, line in calling_unphased.items(): print(f"{sample}: {'+'.join(line['all'])}/")
    # EXPERIMENT 2.2: use homozygous variants alleles
    # sel_samples = [sample for sample in samples_unphased.keys() if sample.split('_')[1] == 'hom'] 
    # calling_unphased = star_allele_calling_all(sel_samples, *pruned_samples_simple, functions, supremal_extended | supremal_samples, reference, detail_level=1)
    # for sample, line in calling_unphased.items(): print(f"{sample}: {'+'.join(line['hom'])}/")
    # EXPERIMENT 2.3: use heterozygous variants alleles
    # sel_samples = [sample for sample in samples_unphased.keys() if sample.split('_')[1] == 'het'] 
    # calling_unphased = star_allele_calling_all(sel_samples, *pruned_samples_simple, functions, supremal_extended | supremal_samples, reference, detail_level=1)
    # for sample, line in calling_unphased.items(): print(f"{sample}: {'+'.join(line['het'])}/")

    # EXPERIMENT 3: unphased star allele calling and trying to infer phasing
    if unphased:
        # TODO change by detail?
        # TODO use extended
        calling = star_allele_calling_all(samples_unphased, *pruned_samples_simple, functions, supremal_extended | supremal_samples, reference, phased=False, detail_level=detail, reorder=text is not None)

        # TEST 7: validate alternative callings
        # validate_alternative_calling(r"results\calling\calling_alt_all.txt", r"data\bastard.txt")

    # Output as text
    if text and visual or text and interactive or visual and interactive:
        raise ValueError("Only one of text, visual or interactive can be True") 
    if text:
        # TODO filter selection
        for sample, line in calling.items(): 
            if type(select) == list and not any([sample in s for s in select]): 
                continue
            print(f"{sample}: {'+'.join(line['A'])}/{'+'.join(line['B'])}")

        # TEST 8: validate calling
        validate_calling(calling, r"data\bastard.txt") # compare to M&J method

    # VISUALISATION 1: Visualise a specific calling and its context
    if visual:
        # TODO change by detail
        if type(select) != list or len(select) != 1: # TODO handle multiple?
            raise Exception("Only one sample can be selected for visualisation")
        select = select[0]
        sample, phase = select.split('_')
        # Check the relevance of the extra variant
        variants_relevance = relevance(select, *pruned_samples_extended, functions, supremal_extended | supremal_samples, reference) # TODO remove relevance?
        group = variants_relevance.keys() if len(variants_relevance.keys()) > 3 else None
        # Mark the star-allele calling
        marked_calling = None
        if phase in "AB": 
            marked_calling = calling[sample][phase]
        else:
            marked_calling = calling[sample]['A'] + calling[sample]['B']
        # Find extend context (including core alleles)
        nodes, edges = find_context({select,}, pruned_samples_extended[1], extend=True, extended=set(), directional=True)
        # Find homozygous
        sel_samples = [sample for sample in samples_unphased.keys() if sample.split('_')[1] == 'hom'] 
        sel_calling = star_allele_calling_all(sel_samples, *pruned_samples_extended, functions, supremal_extended | supremal_samples, reference, detail_level=4)
        homozygous = set([allele for allele in sel_calling[sample]['hom'] if allele != "CYP2D6*1"])
        homozygous, _ = find_context(homozygous, pruned_samples_extended[1], extend=True, extended=set(), directional=True)
        # TODO taxi edges?
        display_graph(nodes, edges, data, functions, default_layout="breadthfirst", auto_download=select if download else None, relevance=variants_relevance, marked_calling=marked_calling, group_variants=group, sample=select, homozygous=homozygous)
        
    # VISUALISATION 2: Show all relations of PharmVar
    if interactive:
        # TODO include relevance?
        edges = pruned_extended[1]
        if type(select) == list:
            edges.extend(find_context(set(select), pruned_samples_extended[1])[1])
        nodes = set([edge[0] for edge in edges] + [edge[1] for edge in edges])
        display_graph(nodes, edges, data, functions)

    # VISUALISATION 3: Generate images for report
    if example:
        config = image_configs[example]
        # TODO allow for multiple
        nodes = config["selection"]
        if "positions" in config:
            positions = config["positions"]
        else:
            positions = None
        if "edges" in config:
            edges = config["edges"]
        else:
            edges = [] # TODO allow for edge selection
            for s, t, r in pruned_extended[1]:
                if s not in nodes or t not in nodes:
                    continue
                edges.append((s, t, r))
        display_graph(nodes, edges, data, functions, default_layout="preset", positions=positions, auto_download=example)

if __name__ == "__main__":
    arguments_parser = argparse.ArgumentParser(description='Star allele calling')
    group_output = arguments_parser.add_mutually_exclusive_group()
    group_output.add_argument('-v', '--visual', action=argparse.BooleanOptionalAction, help="Output calling as image")
    group_output.add_argument('-i', '--interactive', action=argparse.BooleanOptionalAction, help="Run an interactive visualisation")
    group_output.add_argument('-t', '--text', action=argparse.BooleanOptionalAction, help="Output calling as text") 
    group_output.add_argument('-e', '--example', type=str, help="Output specific example image") 
    group_input = arguments_parser.add_mutually_exclusive_group()
    group_input.add_argument('-p', '--phased', action=argparse.BooleanOptionalAction, help="Phased star allele calling")
    group_input.add_argument('-u', '--unphased', action=argparse.BooleanOptionalAction, help="Unphased star allele calling")
    arguments_parser.add_argument('-s', '--select', type=str, nargs='+', default=None, help='Selection of samples to call')
    arguments_parser.add_argument('--detail', type=int, default=0, help="Output detail level") 
    arguments_parser.add_argument('--download', action=argparse.BooleanOptionalAction, help="Download image")
    arguments = vars(arguments_parser.parse_args())
    main(**arguments)