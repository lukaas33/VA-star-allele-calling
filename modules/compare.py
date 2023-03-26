import algebra as va
from .utils import printProgressBar
from .data import cache_get, cache_set
import algebra as va
import multiprocessing as mp
from multiprocessing.managers import SharedMemoryManager
from textwrap import wrap
from sys import getsizeof
import math

def find_relation(args):
    """Worker for multiprocessing relations.
    
    Finds all relations for a given left variant.
    This is faster than finding the relation for a pair since sharing memory costs time.
    """
    start, left, ref_chunks, sequences, count = args
    # Find relation for reconstructed variants
    relations = []
    reference = "".join(ref_chunks)
    lhs = va.Variant(sequences[left*3], sequences[left*3+1], sequences[left*3+2])
    for right in range(start, len(sequences)//3):
        rhs = va.Variant(sequences[right*3], sequences[right*3+1], sequences[right*3+2])
        relation = va.relations.supremal_based.compare(reference, lhs, rhs)
        relations.append((left, right, relation)) 
    # Print progress
    count[0] -= 1.0
    printProgressBar(count[1] - count[0], count[1], prefix = 'Comparing:', length = 50)
    return relations

def find_relations_all(reference_sequence, right_variants, left_variants={}, cache_name=None):
    """Find the relation between all corealleles, suballeles and variants.

    Variants should be stored as a dictionary of supremal representations.
    If left_variants is None, it is considered to be all left variants (all against each other).   
    Can be cached to avoid repeating.

    Returns edge list.
    """
    # TODO apply reduction when comparing to calculate fewer relations
    try:
        if cache_name: return cache_get(cache_name)
    except:
        pass
    # Parse and get relations
    relations = []
    all_variants = left_variants | right_variants
    variant_names = list(all_variants.keys())
    with SharedMemoryManager() as smn:
        # Store data between processes
        # Divide into chunks to allow storage in shared memory (10MB max)
        # TODO is there a better way to do this?
        size = getsizeof(reference_sequence)   
        chunks = wrap(reference_sequence, size // math.ceil(size / 10e6))
        ref = smn.ShareableList(chunks)
        # Spread properties since variant can't be stored in shared memory
        spread = [] 
        for supremal in all_variants.values():
            print(supremal)
            spread += [supremal.start, supremal.end, supremal.sequence]
        seqs = smn.ShareableList(spread)
        total = len(left_variants)
        if total == 0: total = len(right_variants)
        count = smn.ShareableList([total, total])
        # Multiprocessing of relations
        with mp.Pool(mp.cpu_count()) as pool:
            n = len(right_variants)
            if len(left_variants) == 0: # Only calculate one direction since the other can be found 
                args = ((i, i, ref, seqs, count) for i in range(n)) 
            else: 
                start = len(left_variants)
                args = ((start, i, ref, seqs, count) for i in range(start)) # Exclude some from the right side
            relations_2D = pool.map(find_relation, args)
            # Store relations
            for relations_1D in relations_2D:
                for i, j, relation in relations_1D:
                    left, right = variant_names[i], variant_names[j]
                    if relation == va.Relation.CONTAINS:
                        inv_relation = va.Relation.IS_CONTAINED
                    elif relation == va.Relation.IS_CONTAINED:
                        inv_relation = va.Relation.CONTAINS
                    else:
                        inv_relation = relation
                    relations.append((left, right, relation))
                    relations.append((right, left, inv_relation))

    if cache_name: cache_set(relations, cache_name)
    return relations