import warnings
import argparse
from modules.data import reference_get, pharmvar_get
from modules.graph import display_graph
from modules.compare import find_relations_all
from modules.relations import prune_relations, find_context, redundant_reflexive
from modules.parse import extract_variants, to_supremal, parse_samples, samples_to_supremal
from modules.calling import star_allele_calling_all, find_type, Type
from modules.other_sources import get_personal_ids, get_personal_impacts
from modules.assets.generate_images import image_configs
import algebra as va
import modules.tests as tests

def main(text, visual, example, select, interactive, phased, unphased, detail, download):
    # Get the reference sequence relevant for the (current) gene of interest
    print("Get reference sequence...")
    reference_name = "NC_000022.11"
    reference_sequence = reference_get(reference_name)
    reference = {"name": reference_name, "sequence": reference_sequence}
    # List genes as symbols in Pharmvar
    genes = pharmvar_get("genes/list") 
    # All information associated with the (current) gene of interest
    print("Get gene information...")
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
    print("Find relations...")
    supremal_extended = extract_variants(reference_sequence, corealleles, suballeles, cache_name="supremal_extended")
    relations_extended = set(find_relations_all(reference_sequence, supremal_extended, cache_name="relations_extended"))
    pruned_extended = prune_relations(relations_extended, cache_name="relations_pruned_extended")
    pruned_extended[0].add("CYP2D6*1") # Add since it won't be found in the relations
    supremal_simple = extract_variants(reference_sequence, corealleles, cache_name="supremal_simple")
    relations_simple = set(find_relations_all(reference_sequence, supremal_simple, cache_name="relations_simple"))
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
    # make_samples_unphased(reference)
    print("Parse samples...")
    samples_phased = parse_samples("data/samples", reference, phased=True, cache_name="samples_phased") 
    samples_unphased = parse_samples("data/samples_unphased", reference, phased=False, cache_name="samples_unphased") 
    supremal_samples, homozygous = samples_to_supremal(samples_phased, samples_unphased, reference, supremal_extended, "supremal_samples")

    # Split into personal variants and samples
    # TODO fix for personal
    personal_variants = {variant: value for variant, value in supremal_samples.items() if find_type(variant) == Type.P_VAR} 
    samples = {sample: value for sample, value in supremal_samples.items() if find_type(sample) == Type.SAMPLE} 

    # TEST 4: check if more information can be found about personal variants.
    ids |= get_personal_ids(personal_variants, reference, cache_name="ids_personal")
    functions |= get_personal_impacts(personal_variants, ids, reference, cache_name="impacts_personal")
 
    # Find all relations with samples
    print("Find relations with samples...")
    relations_samples_extended = set(find_relations_all(reference_sequence, supremal_extended, samples, cache_name="relations_samples_extended"))
    relations_samples_extended |= set(find_relations_all(reference_sequence, samples, personal_variants, cache_name="relations_samples_personal"))
    relations_samples_extended |= set(find_relations_all(reference_sequence, supremal_extended, personal_variants, cache_name="relations_personal_extended"))
    relations_samples_extended |= set(find_relations_all(reference_sequence, personal_variants, cache_name="relations_personal_ex"))
    # Simplified
    relations_samples_simple = set(find_relations_all(reference_sequence, supremal_simple, samples, cache_name="relations_samples_simple"))
    relations_samples_simple |= set(find_relations_all(reference_sequence, samples, personal_variants, cache_name="relations_samples_personal"))
    relations_samples_simple |= set(find_relations_all(reference_sequence, supremal_simple, personal_variants, cache_name="relations_personal_simple"))
    relations_samples_simple |= set(find_relations_all(reference_sequence, personal_variants, cache_name="relations_personal_si"))

    # Manually remove unparsable TODO fix the parsing of these and remove this
    unparsable_samples = ("HG00373_all","HG01680_all","NA18526_all" ,"NA18632_all","NA19095_all","NA19908_all","NA20289_all","NA20296_all")
    for s in unparsable_samples:
        if s in samples:
            del samples[s]
        if s in supremal_samples:
            del supremal_samples[s]
    to_remove_ex = set()
    for s, t, r in relations_samples_extended:
        if s in unparsable_samples or t in unparsable_samples:
            to_remove_ex.add((s, t, r))
    relations_samples_extended -= to_remove_ex
    to_remove_si = set()
    for s, t, r in relations_samples_simple:
        if s in unparsable_samples or t in unparsable_samples:
            to_remove_si.add((s, t, r))
    relations_samples_simple -= to_remove_si

    # TEST 5: check if relations are consistent with atomic variants
    # test_variant_containment(corealleles, suballeles, relations_extended)
    # test_personal_variant_containment(samples_source, relations_samples)
    # test_central_personal_variants(personal_variants.keys(), find_relations_all(reference_sequence, samples, personal_variants, cache_name="relations_samples_personal"))
    # test_central_personal_variants(personal_variants.keys(), relations_samples)

    # Simplify sample relations
    print("Simplify sample relations...")
    pruned_samples_extended = prune_relations(relations_samples_extended | relations_extended, cache_name="relations_pruned_samples_extended")
    pruned_samples_simple = prune_relations(relations_samples_simple | relations_simple, cache_name="relations_pruned_samples_simple")

    # TEST 6: test if star allele based on corealleles is the same as calling with suballeles
    # test_extended_simplified(samples, pruned_samples_simple, supremal_simple, pruned_samples_extended, supremal_extended, supremal_samples, functions, reference)

    # EXPERIMENT 1: Determine star allele calling for phased samples
    if unphased and phased: 
        raise ValueError("Either phased or unphased should be True")
    if phased:
        print("Calling phased...")
        # TODO change used by detail?
        calling = star_allele_calling_all(samples_phased, *pruned_samples_extended, functions, supremal_extended | supremal_samples, reference, detail_level=detail, reorder=text is not None)
    
    # EXPERIMENT 2: Determine star allele calling for unphased samples
    # EXPERIMENT 2.1: use all variants in single allele
    # sel_samples = [sample for sample in samples_unphased.keys() if sample.split('_')[1] == 'all'] 
    # calling_unphased = star_allele_calling_all(sel_samples, *pruned_samples_extended, functions, supremal_extended | supremal_samples, reference, detail_level=detail)
    # for sample, line in calling_unphased.items(): print(f"{sample}: {','.join(line['all'])}/")
    # EXPERIMENT 2.2: use homozygous variants alleles
    # sel_samples = [sample for sample in samples_unphased.keys() if sample.split('_')[1] == 'hom'] 
    # calling_unphased = star_allele_calling_all(sel_samples, *pruned_samples_extended, functions, supremal_extended | supremal_samples, reference, detail_level=detail)
    # for sample, line in calling_unphased.items(): print(f"{sample}: {','.join(line['hom'])}/")

    # EXPERIMENT 3: unphased star allele calling and trying to infer phasing
    if unphased:
        # TODO change by detail?
        print("Calling unphased...")
        calling = star_allele_calling_all(samples_unphased, *pruned_samples_extended, functions, supremal_extended | supremal_samples, reference, homozygous=homozygous, phased=False, detail_level=detail, reorder=text is not None)

    # Statistics
    # tests.statistics(corealleles, suballeles, relations_extended, pruned_extended[1], calling)

    # TEST 7: validate alternative callings
    # tests.validate_alternative_calling(r"results\calling\calling_alt.txt", r"results/calling/calling_phased.txt")
    # also check suballele callings
    tests.validate_alternative_calling(r"results\calling\calling_alt_sub+.txt", r"results/calling/calling_phased_sub.txt")

    # TEST 8: validate alternative calling method on simulated data
    # tests.test_alternative_callings(supremal_extended, reference, relations_extended, functions)

    # Output as text
    if text and visual or text and interactive or visual and interactive:
        raise ValueError("Only one of text, visual or interactive can be True") 
    if text:
        print(f"Outputting text (detail={detail})...")
        for sample, line in calling.items(): 
            if type(select) == list and not any([sample in s for s in select]): 
                continue
            print(f"{sample}: {','.join(line['A'])}/{','.join(line['B'])}")

        # TEST 8: validate calling
        # validate_calling(calling, r"data\bastard.txt") # compare to M&J method

    # VISUALISATION 1: Visualise a specific calling and its context
    if visual:
        print("Visualising calling...")
        # TODO change by detail
        if type(select) != list or len(select) != 1: # TODO handle multiple?
            raise Exception("Only one sample can be selected for visualisation")
        select = select[0]
        sample, phase = select.split('_')
        # Group extra variants
        group = set()
        for var1, var2, _ in pruned_samples_extended[1]:
            if var1 != select and var2 != select:
                continue
            var = var2 if var1 == select else var1
            if find_type(var) not in (Type.VAR, Type.P_VAR):
                continue            
            group.add(var)
        if len(group) <= 5:
            group = None
        # Mark the star-allele calling
        marked_calling = None
        if phase in "AB": 
            marked_calling = calling[sample][phase]
        else:
            marked_calling = calling[sample]['A'] + calling[sample]['B']
        # Find extended context
        # TODO leave out overlap? (can mess up dagre)
        nodes, edges = find_context({select,}, pruned_samples_extended[1], extend=True, extended=set(), directional=True, overlap=False)
        # Find homozygous
        sel_samples = [sample for sample in samples_unphased.keys() if sample.split('_')[1] == 'hom'] 
        sel_calling = star_allele_calling_all(sel_samples, *pruned_samples_extended, functions, supremal_extended | supremal_samples, reference, detail_level=4)
        homozygous_alleles = set([allele for allele in sel_calling[sample]['hom'] if allele != "CYP2D6*1"])
        homozygous_alleles = find_context(homozygous_alleles, pruned_samples_extended[1], directional=True, overlap=False, extend=True, extended=set())[0]
        homozygous_alleles = set((a for a in homozygous_alleles if find_type(a) in (Type.SUB, Type.CORE)))
        # TODO taxi edges?
        display_graph(nodes, edges, data, functions, default_layout="dagre", auto_download=select if download else None, relevance=None, marked_calling=marked_calling, group_variants=group, sample=select, homozygous=homozygous[sample] | homozygous_alleles)
        
    # VISUALISATION 2: Show all relations of PharmVar
    if interactive:
        print("Interactive map...")
        edges = set(pruned_extended[1])
        if type(select) == list:
            edges |= find_context(set(select), pruned_samples_extended[1])[1]
        nodes = set([edge[0] for edge in edges] + [edge[1] for edge in edges])
        display_graph(nodes, edges, data, functions)

    # VISUALISATION 3: Generate images for report
    if example:
        config = image_configs[example]
        print(f"Visualising example {example}...")
        homozygous = None if "homozygous" not in config else config["homozygous"]
        # TODO allow for multiple
        nodes = config["selection"]
        layout = "cose-bilkent"
        if "layout" in config:
            layout = config["layout"]
        if "positions" in config:
            positions = config["positions"]
            layout = "preset"
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
        display_graph(nodes, edges, data, functions if config["color"] else None, default_layout=layout, positions=positions, auto_download=example, homozygous=homozygous)

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
    arguments_parser.add_argument('--detail', type=int, default=1, help="Output detail level") 
    arguments_parser.add_argument('--download', action=argparse.BooleanOptionalAction, help="Download image")
    arguments = vars(arguments_parser.parse_args())
    main(**arguments)