import algebra as va
import warnings

def parse_multi_hgvs(hgvs_lst, reference_sequence, allele_name):
    """Wrapper for parsing a list of hgvs variants and combining them into a set of variants (allele).

    Sometimes variants have the same hgvs representation but a different position value,
    this is because there a multiple ways to optimally place the variant which have the same hgvs representation
    For variant algebra this is not relevant so these are removed by using a set.
    """        
    variants = set()
    for hgvs in hgvs_lst:
        try:
            for variant in va.variants.parse_hgvs(hgvs, reference=reference_sequence): # Reference needed to handle insertions
                variants.add(variant)
        except: 
            raise ValueError(f"HGVS string '{hgvs}' could not be parsed.")
    return variants

def combine_variants(variants, reference_sequence):
    """ Combine multiple variants into one larger variant (allele) in the supremal representation
    
    WARNING: this is redundant for comparing since the compare function already patches the variants
    """
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