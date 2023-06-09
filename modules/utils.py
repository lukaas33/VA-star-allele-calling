import algebra as va
import difflib
from itertools import chain, combinations
from .data import api_get
import warnings
import os
import vcf

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
                

def make_samples_unphased(reference):
    """Make all samples unphased"""
    # TODO retire
    directory = "data/samples"
    new_directory = "data/samples_unphased"
    if not os.path.exists(new_directory):
        os.mkdir(new_directory)
    for filename in os.listdir(directory):
        with open(os.path.join(directory, filename), 'r') as file:
            data = file.read()
        # Find variants with overlapping variants
        variants = []
        info = []
        with open(os.path.join(directory, filename), 'r') as file:
            reader = vcf.Reader(file)
            for record in reader:
                var = va.Variant(record.start, record.end, record.ALT[0].sequence)
                if var in variants: # Ignore duplicates here, will not be changed
                    continue
                info.append((var, record))
        overlapping = []
        for i, (var, record) in enumerate(info):
            for var2, record2 in info[i+1:]:
                phase, phase2 = record.samples[0].data[0], record2.samples[0].data[0]
                id, id2 = record.ID, record2.ID
                if max(var.start, var2.start) <= min(var.end, var2.end)-1: # Overlap in positions
                    if phase == phase2: # Same phase, cannot coexist
                        continue
                    # Overlapping in two phases can exist
                    overlapping.append((record, record2))
        # Combine overlapping
        if len(overlapping) > 0:
            for record, record2 in overlapping:
                # print(f"{filename} has overlapping variants {record.ID} and {record2.ID}")
                for line in data.split('\n'):
                    if record.ID in line: # Change notation of this line
                        # Find reference sequence
                        pos1, pos2 = record.start, record2.start
                        end1, end2 = record.end, record2.end
                        ref1, ref2 = record.REF, record2.REF
                        change1, change2 = record.ALT[0].sequence, record2.ALT[0].sequence
                        # print(pos1, end1, ref1, change1)
                        # print(pos2, end2, ref2, change2)
                        min_pos = min(pos1, pos2)
                        max_pos = max(end1, end2)
                        ref = reference["sequence"][min_pos:max_pos]
                        if len(ref) > len(ref1):
                            change1 = change1 + ref[len(ref1):]
                        if len(ref) > len(ref2):
                            change2 = change2 + ref[len(ref2):]
                        # Find change and scale to reference
                        columns = line.split('\t')
                        # print('\t'.join(columns[:5]))
                        columns[1] = str(min(pos1, pos2)+1) # Lowest position
                        columns[2] = record.ID + '+' + record2.ID
                        columns[3] = ref # Reference sequence
                        columns[4] = f"{change1},{change2}" # Combine sequences TODO fix for differing positions
                        columns[-1] = "1/2" # Notation heterozygous unphased
                        # TODO should also combine other fields but not necessary for now
                        data = data.replace(line, '\t'.join(columns))
                        # print('\t'.join(columns[:5]))
                        # print()
                    elif record2.ID in line: # Remove line
                        data = data.replace(line, "")
        # Replace all phased with unphased
        data = data.replace("1|1", "1/1")
        data = data.replace("0|1", "0/1")
        data = data.replace("1|0", "0/1")
        with open(os.path.join(new_directory, filename), 'w') as file:
            file.write(data)

temp_cache = {}
def normalise(hgvs, ref="NC_000022.11"):
    global temp_cache
    hgvs = f"{ref}:g.{hgvs}" if ':' not in hgvs else hgvs
    if hgvs in temp_cache:
        return temp_cache[hgvs]
    else:
        lookup = api_get(f"https://mutalyzer.nl/api/normalize/{hgvs}")
        temp_cache[hgvs] = lookup['normalized_description']
        return lookup['normalized_description']

def change_ref(hgvs, from_ref="NC_000022.11", to_ref="NM_000106.6"):
    "Change the reference of a hgvs"
    raise NotImplementedError("Not implemented yet")
    inv = {"A": "T", "T": "A", "C": "G", "G": "C"}
    hgvs = f"{from_ref}:g.{hgvs}" if ':' not in hgvs else hgvs
    lookup = api_get(f"https://mutalyzer.nl/api/normalize/{hgvs}")
    for eq, _ in lookup['equivalent_descriptions']['c']:
        if to_ref in eq:
            eq = eq.split('(')[1].split(')')[0] + eq.split('(')[1].split(')')[1]
            return eq
