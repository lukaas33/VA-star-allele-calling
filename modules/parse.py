import algebra as va
import warnings

def parse_multi_hgvs(hgvs_lst, reference):
    """Wrapper for parsing a list of hgvs variants and combining them into a set of variants (allele).

    Also includes preprocessing of data which involves fixing HGVS position notation and checking for wrong sets of variants.
    """        
    variant_lst = []
    for hgvs in hgvs_lst:
        hgvs = fix_hgvs_position(hgvs) # HGVS preprocessing
        try:
            variant_lst += va.variants.parse_hgvs(hgvs, reference=reference) # Reference needed to handle insertions
        except: 
            raise ValueError(f"HGVS string '{hgvs}' could not be parsed.")
    variant_set = set(variant_lst)
    # Test if any variants were equivalent
    if len(variant_set) != len(variant_lst): 
        raise ValueError(f"Double variant positions in {hgvs_lst}")
    variant_set = fix_variant_overlapping(variant_lst) # Variant preprocessing
    return variant_set  

def fix_hgvs_position(hgvs):
    """Fix HGVS notation problems with the position notation 
    
    Problems occur for deletions where a range is not specified but the (normally redundant) deleted area is.

    These notations is not HGVS compliant: 
    https://varnomen.hgvs.org/recommendations/DNA/variant/deletion/
    But they are present in the Pharmvar database so they are preprocessed here.
    """
    # TODO integrate change in va parser?
    if "del" in hgvs:
        area = hgvs.split("del")[1]
        position = hgvs.split("del")[0].split('.')[-1]
        if len(area) > 1 and '_' not in position:
            range = f"{position}_{int(position)+len(area)-1}"
            corrected_hgvs = hgvs.replace(position, range)
            warnings.warn(f"Illegal but interpretable HGVS notation used '{corrected_hgvs}'. Correcting and continuing.")
            return corrected_hgvs
    return hgvs

def fix_variant_overlapping(variants):
    """Fix overlapping variants in an variant list.
    
    This is needed since overlapping variants cannot be ordered and thus not compared.
    """
    # QUESTION why can't unorderable variants be compared?
    # TODO Consider adjacent (>=) as overlapping? Adjacent HGVS-proof but doesn't cause va problems.
    sorted_variants = sorted(list(variants), key=lambda v: v.start)
    for i in range(0, len(sorted_variants)-1):
        if sorted_variants[i].end > sorted_variants[i+1].start:
            combined_sequence = sorted_variants[i].sequence[:(1 + sorted_variants[i+1].start - sorted_variants[i].start)]
            combined_sequence += sorted_variants[i+1].sequence[(sorted_variants[i].end - sorted_variants[i+1].start):]
            combined_variant = va.Variant(sorted_variants[i].start, sorted_variants[i+1].end, combined_sequence)
            warnings.warn(f"Variants {sorted_variants[i]} and {sorted_variants[i+1]} overlap. Combining into {combined_variant} and continuing.")
            sorted_variants[i] = None
            sorted_variants[i+1] = combined_variant

    non_overlapping_variants = set([var for var in sorted_variants if var is not None])
    return non_overlapping_variants

def combine_variants(variants, reference_sequence):
    """ Combine multiple variants into one larger variant (allele)
    
    WARNING: this is redundant for comparing since the compare function already patches the variants
    """
    # QUESTION: is this a correct method?
    # TODO make representation more minimal
    #       using supremal representation gives a longer sequence
    # TODO replace this function with calls to the va library (see compare method for the proper calls)
    # Apply variants to reference_sequence to get the observed sequence
    min_start = float('inf')
    max_end = -float('inf')
    observed_sequence = va.variants.patch(reference_sequence, variants)
    for operation in variants:
        if operation.start < min_start:
            min_start = operation.start
        elif operation.end > max_end:
            max_end = operation.end
    # Find difference and output
    allele = va.Variant(min_start, max_end, observed_sequence[min_start:max_end])
    return allele