from modules.utils import count_relations, count_arity
from modules.data import cache_get, cache_set, api_get
import warnings
from modules.calling import star_allele_calling_all, find_type, Type
import algebra as va
from modules.other_sources import is_silent_mutalyzer, get_annotation_entrez, find_id_hgvs, get_personal_ids, get_personal_impacts, impact_position, relevance
from itertools import combinations
from modules.parse import to_supremal
from modules.relations import prune_relations

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

def test_alternative_callings(supremals, reference, relations_ref, functions):
    """Test alternative calling method on simulated data"""
    variants = {"NC_000022.11:g.42126611C>G", "NC_000022.11:g.42127941G>A", "NC_000022.11:g.42130692G>A"}
    alleles = {}
    i = 0
    for r in range(0, len(variants)+1):
        for comb in combinations(variants, r):
            i += 1
            for t, vs in (("all", variants), ("hom", comb)):
                alleles[f"HG{i}_{t}"] = []
                print(f"HG{i}_{t}:", end=' ')
                for o in vs:
                    print(f"{o}", end=', ')
                    alleles[f"HG{i}_{t}"].extend(va.variants.parse_hgvs(o))
                print()
                if len(alleles[f"HG{i}_{t}"]) == 0:
                    alleles[f"HG{i}_{t}"] = None
                    continue
                alleles[f"HG{i}_{t}"] = to_supremal(alleles[f"HG{i}_{t}"], reference["sequence"])
    try:
        relations = cache_get("TEST_relations")
    except:
        relations = []
        for o_allele, o_supremal in alleles.items():
            c = 0
            for r_allele, r_supremal in supremals.items():
                if o_supremal is None or r_supremal is None:
                    continue
                rel = va.relations.supremal_based.compare(reference["sequence"], o_supremal, r_supremal)
                relations.append((o_allele, r_allele, rel))
                if rel == va.Relation.CONTAINS:
                    rel = va.Relation.IS_CONTAINED
                elif rel == va.Relation.IS_CONTAINED:
                    rel = va.Relation.CONTAINS
                relations.append((r_allele, o_allele, rel))
                c += 1
                print(c, len(supremals))
        cache_set(relations, "TEST_relations")
    supremals = supremals | alleles
    relations += relations_ref 
    pruned = prune_relations(relations)
    star_allele_calling_all(alleles.keys(), *pruned, functions, supremals, reference, phased=False, detail_level=1)

def statistics(corealleles, suballeles, relations, pruned_relations, callings):
    print("Core alleles:", len(corealleles.keys()))
    var_core = list(set([var["hgvs"] for core in corealleles.keys() for var in corealleles[core]["variants"]]))
    print("Core variants:", len(var_core), var_core[:4])
    all_sub = [sub for core in corealleles.keys() for sub in suballeles[core]]
    print("Suballeles:", len(all_sub))
    var_all = list(set([var["hgvs"] for core in corealleles.keys() for sub in suballeles[core] for var in suballeles[core][sub]["variants"]]))
    print("All variants:", len(var_all), var_all[:5])
    print("Alleles:", len(corealleles.keys()) + len(all_sub))
    print("Theoretical relation count:", (len(corealleles.keys()) + len(all_sub) + len(var_all))**2)
    print("Variants:", len(set(var_core + var_all)))
    print("Relations before pruning")
    print(count_relations(relations))
    print("Relations after pruning")
    print(count_relations(pruned_relations))
    print("Most common nodes:")
    arity = count_arity(corealleles.keys(), pruned_relations)
    print(sorted([(a["total"], n) for n, a in arity.items()], reverse=True)[:5])
    counts = {}
    calling_count = {"hom": 0, "het": 0, "1def": 0, "2def": 0}
    for sample, calling in callings.items():
        if sample == "NA18526": # skip unparsable TODO fix
            continue
        for phase in "AB":
            name = f"{sample}_{phase}"
            if calling[phase][0] not in counts:
                counts[calling[phase][0]] = 0
            counts[calling[phase][0]] += 1
        if calling['A'][0] == calling['B'][0] == "CYP2D6*1":
            calling_count["2def"] += 1
        elif calling['A'][0] == "CYP2D6*1" or calling['B'][0] == "CYP2D6*1":
            calling_count["1def"] += 1
        elif calling['A'] == calling['B']:
            calling_count["hom"] += 1
        else:
            calling_count["het"] += 1
    assert sum(counts.values()) == 240
    counts = list(counts.items()) 
    counts.sort(key=lambda x: x[1], reverse=True)
    print(counts[:10])
    assert sum(calling_count.values()) == 120
    print(calling_count)

  
def validate_alternative_calling(calling_filename, validate_filename):
    """Validate if alternative calling matches with the M&J method"""
    # TODO integrate differently (not using file)
    validate = {}
    with open(validate_filename, 'r', encoding="utf-16") as file:
        for line in file:
            # convert to format of classifications
            sample, calling = line.rstrip().split(': ')
            calling = calling.split('/')
            # calling.sort(key=lambda x: int(x.split('*')[1]))
            # calling = ["CYP2D6" + c for c in calling]
            validate[sample] = calling
    to_find = set(validate.keys())

    def end_sample(sample, prev, found, count, totals, to_find):
        if found != -1:
            print(f"{sample} has the correct calling as the {found}th allele out of {count} alternatives ({'preferred' if found == 1 else 'not preferred'})")
            if found == 1: totals['preferred'] += 1
            else: totals['not_preferred'] += 1
        else:
            if "CYP2D6*?" in prev:
                print(f"{sample} has an unparsable calling.")
                totals['unparsable'] += 1
            else:
                print(f"{sample} has an incorrect calling. Last was {prev} (out of {count}). Calling should be {validate[sample]}")
                totals['incorrect'] += 1
        to_find.remove(sample)

    # TODO make into single variable
    totals = {
        'preferred': 0,
        'not_preferred': 0,
        'incorrect': 0,
        "unparsable": 0
    }
    with open(calling_filename, 'r', encoding="utf-16") as file:
        line = next(file).strip()
        sample = line
        prev = None
        found = -1
        count = 0
        for line in file:
            line = line.strip()
            if not line: 
                continue
            if line.startswith('CYP2D6'):
                line = line.split(' ')[0] # Handle annotations
                prev = line
                count += 1
                if line.split('/') == validate[sample]:
                    found = count
            else:
                end_sample(sample, prev, found, count, totals, to_find)
                sample = line
                found = -1
                count = 0
                prev = None
        end_sample(sample, prev, found, count, totals, to_find)

    print(f"{totals['preferred']} samples correct as preferred allele")
    print(f"{totals['not_preferred']} samples not correct preferred allele")
    print(f"{totals['incorrect']} samples incorrect")
    print(f"{totals['unparsable']} samples unparsable")
    if len(to_find) > 0:
        print(f"{len(to_find)} not in alternative callings:")
        for sample in to_find:
            print(f"\t{sample} should be {validate[sample]}")
    print(f"{sum(totals.values()) + len(to_find)} total")


def validate_relations(data, variants, filename):
    """Validate if relations match with the M&J method"""
    # TODO move
    data = set([
        (left, right, rel) 
            for left, right, rel in data 
            if rel != va.Relation.DISJOINT \
                and left != right]
    )
    variants = {v["variantId"]: v["hgvs"] for v in variants.values()}
    ref = set()
    wrong_ref = set()
    not_in_ref = set(data)
    not_in_data = set()
    with open(filename) as file:
        for line in file:
            # Convert to same edge notation as data
            edge = line.rstrip().split(' ')
            edge[2] = va.Relation[edge[2].upper()] 
            for i in range(2):
                # Convert variant_id notation to HGVS
                if 'variant_' in edge[i]: 
                    id = edge[i].split('variant_')[1]
                    if id not in variants.keys():
                        warnings.warn(f"Variant {id} not found in variants")
                        wrong_ref.add(tuple(edge))
                        break
                    edge[i] = variants[id]
            else: # No break, can continue
                # Find reversed edge
                reversed = [edge[1], edge[0], edge[2]] # Reverse
                if edge[2] == va.Relation.CONTAINS: reversed[2] = va.Relation.IS_CONTAINED
                elif edge[2] == va.Relation.IS_CONTAINED: reversed[2] = va.Relation.CONTAINS
                # Check for orientation that is in data
                edge = tuple(edge)
                reversed = tuple(reversed)
                if edge in data and reversed in data: # Both in data
                    ref.add(edge)
                    ref.add(reversed)
                elif edge in data: # Specific orientation in data
                    ref.add(edge)
                elif reversed in data: # Specific orientation in data
                    ref.add(reversed)
                else:
                    not_in_data.add(edge) # In ref but not in data
                # Remove from data to see what is left at the end
                if edge in not_in_ref:
                    not_in_ref.remove(edge) 
                if reversed in not_in_ref: 
                    not_in_ref.remove(reversed)
    # Relations test
    if len(not_in_data) > 0:
        print("These relations are in the reference but not in the data")
        for n in not_in_data:
            print('\t', n)
    if len(not_in_ref) > 0:
        print("These relations are in the data but not in the reference")
        for n in not_in_ref:
            print('\t', n)
    if len(wrong_ref) > 0:
        print("Reference contains these alleles that are wrong")
        for w in wrong_ref:
            print('\t', w)
    # Count test
    count_ref = count_relations(ref)
    count_data = count_relations(data)
    for pair in zip(count_ref.items(), count_data.items()):
        if pair[0][1] != pair[1][1]:
            print(f"Difference in relation count for {pair[0][0]}: {pair[0][1]} in ref vs {pair[1][1]} in data")
    
    print("No (further) differences found between the reference and the data relations")

def validate_calling(callings, validate_filename, soft=False):
    """Validate if classifications match with the M&J method"""
    # TODO move
    n_errors = 0
    with open(validate_filename, 'r') as validate:
        for line in validate:
            # convert to format of classifications
            sample, reference = line.strip().split(' ')
            reference = ["CYP2D6" + c for c in reference.split('/')]
            calling = [] # Get priority answer (assumes filtering has occurred already with representation method)  
            for c in callings[sample].values():
                calling.extend(c)
            if soft: # Soft check, check if prediction contains correct answers
                if all([r in calling for r in reference]):
                    continue
            else: # Hard check, exactly equal
                if reference == calling or reference[::-1] == calling:
                    continue
            print(f"Sample {sample} was predicted as {calling} but should be {reference}")
            n_errors += 1
    if n_errors > 0:
        print(f"{n_errors} errors found in the classifications")