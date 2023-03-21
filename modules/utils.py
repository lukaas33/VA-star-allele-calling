import algebra as va
import difflib
from itertools import chain, combinations
from .data import cache_get
import warnings

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
    return counts

def count_arity(nodes, relations):
    # TODO not accurate since it doesn't account for unreduced relations
    arity = {node: {relationType.name: 0 for relationType in va.Relation} for node in nodes}
    for l_allele, r_allele, relation in relations:
        arity[l_allele][relation.name] += 1
        # Find inverse
        if relation.name == "IS_CONTAINED": # Directional
            arity[r_allele]["CONTAINS"] += 1
            continue
        if relation.name == "CONTAINS": # Directional
            arity[r_allele]["IS_CONTAINED"] += 1
            continue
        arity[r_allele][relation.name] += 1
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
    data = set([row for row in data if row[2] != va.Relation.DISJOINT and row[0] != row[1]])
    variants = {v["variantId"]: v["hgvs"] for v in variants.values()}
    ref = set()
    wrong_ref = set()
    with open(filename) as file:
        for line in file:
            edge = line.rstrip().split(' ')
            for i in range(2):
                # Convert variant notation to HGVS
                if 'variant_' in edge[i]: 
                    id = edge[i].split('variant_')[1]
                    if id not in variants.keys():
                        warnings.warn(f"Variant {id} not found in variants")
                        break
                    edge[i] = variants[id]
            else: # No break, can continue
                # Convert to same edge notation as data
                edge[2] = va.Relation[edge[2].upper()]
                reversed = [edge[1], edge[0], edge[2]]
                if edge[2] == va.Relation.CONTAINS: reversed[2] = va.Relation.IS_CONTAINED
                elif edge[2] == va.Relation.IS_CONTAINED: reversed[2] = va.Relation.CONTAINS
                # Store in reference
                ref.add(tuple(edge))
                ref.add(tuple(reversed))
    # Find differences
    not_in_ref = data - ref
    not_in_data = ref - data

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
    print("No (further) differences found between the reference and the data relations")