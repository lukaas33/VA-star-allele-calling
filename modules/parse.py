import algebra as va
import warnings

def parse_multi_hgvs(hgvs_lst, reference_sequence, allele_name):
    """Wrapper for parsing a list of hgvs variants and combining them into a set of variants (allele).

    Also includes preprocessing of data which involves fixing HGVS position notation and checking for wrong sets of variants.
    """        
    variant_lst = []
    for hgvs in hgvs_lst:
        hgvs = fix_hgvs_position(hgvs, allele_name) # HGVS preprocessing
        try:
            variant_lst += va.variants.parse_hgvs(hgvs, reference=reference_sequence) # Reference needed to handle insertions
        except: 
            raise ValueError(f"HGVS string '{hgvs}' could not be parsed.")
    variant_set = set(variant_lst)
    # Test if any variants were equivalent
    if len(variant_set) != len(variant_lst): 
        warnings.warn(f"{allele_name}: Double variant hgvs: {[{va.variants.to_hgvs([var], reference='NC000022.11', sequence_prefix=True)} for var in variant_lst if variant_lst.count(var) > 1]}")
    return variant_set  

def fix_hgvs_position(hgvs, allele_name):
    """Fix HGVS notation problems with the position notation 
    
    Problems occur for deletions where a range is not specified but the (normally redundant) deleted area is.

    These notations is not HGVS compliant: 
    https://varnomen.hgvs.org/recommendations/DNA/variant/deletion/
    But they are present in the Pharmvar database so they are preprocessed here.
    """
    # TODO remove since position field is not hgvs
    if "del" in hgvs:
        area = hgvs.split("del")[1]
        position = hgvs.split("del")[0].split('.')[-1]
        if len(area) > 1 and '_' not in position:
            range = f"{position}_{int(position)+len(area)-1}"
            corrected_hgvs = hgvs.replace(position, range)
            warnings.warn(f"{allele_name}: Illegal but interpretable HGVS notation used '{corrected_hgvs}'. Correcting and continuing.")
            return corrected_hgvs
    return hgvs

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