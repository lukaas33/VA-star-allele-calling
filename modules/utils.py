import algebra as va
import difflib
from itertools import chain, combinations
from .data import cache_get
import warnings
import os

def print_seq_diff(sequence1, sequence2, start=1):
    """ Output the difference between two aligned sequences as insertions and deletions. """
    difference = difflib.ndiff(sequence1, sequence2)
    for i,s in enumerate(difference):
        if s[0]==' ': continue
        elif s[0]=='-':
            print(u'Delete "{}" on relative position {}'.format(s[-1],start+i))
        elif s[0]=='+':
            print(u'Insert "{}" on relative position {}'.format(s[-1],start+i)) 

def va_generate_subsets(variants):
    """Generate all subsets (superset) of a variant set"""
    superset = chain.from_iterable(combinations(variants, r) for r in range(1, len(variants)+1))
    return superset

def count_relations(relations):
    counts = {relationType.name: 0 for relationType in va.Relation}
    for _, _, relation in relations:
        counts[relation.name] += 1
    counts["CONTAINS"] = counts["IS_CONTAINED"]
    counts["total"] = sum(counts.values())
    return counts

def count_arity(nodes, relations):
    # TODO not accurate since it doesn't account for unreduced relations
    arity = {node: {relationType.name: 0 for relationType in va.Relation} for node in nodes}
    for l_allele, r_allele, relation in relations:
        if l_allele not in arity or r_allele not in arity:
            continue
        arity[l_allele][relation.name] += 1
        # Find inverse
        if relation.name == "IS_CONTAINED": # Directional
            arity[r_allele]["CONTAINS"] += 1
            continue
        if relation.name == "CONTAINS": # Directional
            arity[r_allele]["IS_CONTAINED"] += 1
            continue
        arity[r_allele][relation.name] += 1
    for node in arity:
        arity[node]["total"] = sum(arity[node].values())
    return arity

def printProgressBar (iteration, total, prefix = '', suffix = '', decimals = 1, length = 100, fill = '█', printEnd = "\r"):
    """
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        length      - Optional  : character length of bar (Int)
        fill        - Optional  : bar fill character (Str)
        printEnd    - Optional  : end character (e.g. "\r", "\r\n") (Str)
    """
    # TODO use tqdm library instead
    # TODO display ETA
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    print(f'\r{prefix} |{bar}| {percent}% {suffix}', end = printEnd, flush=True)
    # Print New Line on Complete
    if iteration == total: 
        print()


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
                

def make_samples_unphased():
    """Make all samples unphased"""
    # TODO retire
    directory = "data/samples"
    new_directory = "data/samples_unphased"
    if not os.path.exists(new_directory):
        os.mkdir(new_directory)
    for filename in os.listdir(directory):
        with open(os.path.join(directory, filename), 'r') as file:
            data = file.read()
        data = data.replace("1|1", "1/1")
        data = data.replace("0|1", "0/1")
        data = data.replace("1|0", "0/1")
        with open(os.path.join(new_directory, filename), 'w') as file:
            file.write(data)


def validate_alternative_calling(calling_filename, validate_filename):
    """Validate if alternative calling matches with the M&J method"""
    # TODO integrate differently (not using file)
    validate = {}
    with open(validate_filename, 'r') as file:
        for line in file:
            # convert to format of classifications
            sample, calling = line.rstrip().split(' ')
            calling = calling.split('/')
            calling.sort(key=lambda x: int(x.split('*')[1]))
            calling = ["CYP2D6" + c for c in calling]
            validate[sample] = calling
    to_find = set(validate.keys())

    def end_sample(sample, prev, how, found, count, totals):
        if found != -1:
            if how == '(hom)':
                print(f"{sample} has the correct calling due to homozygous alleles")
                totals['hom'] += 1
            else:
                print(f"{sample} has the correct calling as the {found}th allele out of {count} alternatives ({'preferred' if found == 1 else 'not preferred'})")
                if found == 1: totals['preferred'] += 1
                else: totals['not_preferred'] += 1
        else:
            print(f"{sample} has an incorrect calling. Last was {prev} (out of {count}). Calling should be {validate[sample]}")
            totals['incorrect'] += 1

    # TODO make into single variable
    totals = {
        'hom': 0,
        'preferred': 0,
        'not_preferred': 0,
        'incorrect': 0
    }
    with open(calling_filename, 'r', encoding="utf-16") as file:
        sample, prev, how = None, None, None
        found = -1
        count = 0
        for _line in file:
            line = _line.strip()
            if not line: 
                continue
            if line.startswith('CYP2D6'):
                line = line.split(' ')[0] # Handle annotations
                prev = line
                count += 1
                if line.split('/') == validate[sample]:
                    found = count
            else:
                sample = line.split(' ')[0]
                how = line.split(' ')[1]
                to_find.remove(sample)
                if prev is None:
                    continue
                end_sample(sample, prev, how, found, count, totals)
                found = -1
                count = 0
        end_sample(sample, prev, how, found, count, totals)

    print(f"{totals['hom']} samples correct due to homozygous alleles")
    print(f"{totals['preferred']} samples correct as preferred allele")
    print(f"{totals['not_preferred']} samples not correct preferred allele")
    print(f"{totals['incorrect']} samples incorrect")
    if len(to_find) > 0:
        print(f"{len(to_find)} not in alternative callings:")
        for sample in to_find:
            print(f"\t{sample} should be {validate[sample]}")
    print(f"{sum(totals.values()) + len(to_find)} total")