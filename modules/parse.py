import algebra as va
import warnings
from .data import cache_get, cache_set

def  to_supremal(variants, reference_sequence):
    """Convert a list of variants to a supremal representation.
    
    Faster because it accounts for single variations
    """
    if len(variants) > 1: # Need to patch
        observed = va.variants.patch(reference_sequence, variants)
        spanning = va.relations.supremal_based.spanning_variant(reference_sequence, observed, variants)
        supremal = va.relations.supremal_based.find_supremal(reference_sequence, spanning)
        return supremal
    else:
        # QUESTION: is the supremal representation the same as the variant?
        supremal = va.relations.supremal_based.find_supremal(reference_sequence, variants[0])
        return supremal

def parse_hgvs_supremal(hgvs_lst, reference_sequence):
    """Parse multiple hgvs as supremal representation. 
    
    Useful for comparing variants faster.
    """
    # Convert to variants
    variants = [] 
    for hgvs in hgvs_lst:    
        variants += va.variants.parse_hgvs(hgvs, reference=reference_sequence)
    # Convert to supremal
    return to_supremal(variants, reference_sequence)

    
def extract_variants(reference_sequence, corealleles, suballeles=None, cache_name=None):
    """Find supremal representations for all variants in the core and suballeles."""
    # TODO use pharmvar data directly?
    try:
        if cache_name: return cache_get(cache_name)
    except:
        pass
    all_variants = {} # Store variants, sub- and corealleles as supremal
    for coreallele in corealleles.keys():
        alleles = [corealleles[coreallele]]
        if suballeles is not None: # Include sub
            alleles += suballeles[coreallele].values()
        for allele in alleles:
            if len(allele["variants"]) == 0:
                warnings.warn(f"Empty variant for {allele['alleleName']}")
                continue
            all_variants[allele["alleleName"]] = [] 
            for variant in allele["variants"]: # Variants for allele
                all_variants[allele["alleleName"]].append(variant["hgvs"])
                if variant["hgvs"] in all_variants.keys(): # Skip if already parsed
                    continue
                all_variants[variant["hgvs"]] = parse_hgvs_supremal([variant["hgvs"]], reference_sequence) # Store variant as supremal
            try:
                all_variants[allele["alleleName"]] = parse_hgvs_supremal(all_variants[allele["alleleName"]], reference_sequence) # Store allele as supremal
            except ValueError as e: # Fails for overlapping variants
                # TODO how to handle duplicates? And how to handle multiple variants at same position?
                error = f"{allele['alleleName']}: {e}"
                warnings.warn(error)
                del all_variants[allele["alleleName"]]
    if cache_name: cache_set(all_variants, cache_name)
    return all_variants
