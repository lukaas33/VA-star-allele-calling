import algebra as va

def parse_hgvs_supremal(hgvs_lst, reference_sequence):
    """Parse multiple hgvs as supremal representation. 
    
    Useful for comparing variants faster.
    """
    # Convert to variants
    variants = [] 
    for hgvs in hgvs_lst:    
        variants += va.variants.parse_hgvs(hgvs, reference=reference_sequence)
    # Convert to supremal
    if len(variants) > 1: # Need to patch
        observed = va.variants.patch(reference_sequence, variants)
        spanning = va.relations.supremal_based.spanning_variant(reference_sequence, observed, variants)
        supremal = va.relations.supremal_based.find_supremal(reference_sequence, spanning)
        return supremal
    else:
        supremal = va.relations.supremal_based.find_supremal(reference_sequence, variants[0])
        return supremal