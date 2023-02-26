import algebra as va

def parse_multi_hgvs(hgvs_lst, reference_sequence):
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
    return list(variants)

