import algebra as va
import warnings
from .data import cache_get, cache_set

def to_supremal(variants, reference_sequence):
    """Convert a list of variants to a supremal representation.
    
    Faster because it accounts for single variations
    """
    if not isinstance(variants, list):
        raise Exception("Variants should be a list")
    if len(variants) == 0: # No variants
        raise Exception("Variant list empty") # TODO empty variant?
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
            vs = [] 
            for variant in allele["variants"]: # Variants for allele
                vs.append(variant["hgvs"])
                if variant["hgvs"] in all_variants.keys(): # Skip if already parsed
                    continue
                all_variants[variant["hgvs"]] = parse_hgvs_supremal([variant["hgvs"]], reference_sequence) # Store variant as supremal
            try:
                all_variants[allele["alleleName"]] = parse_hgvs_supremal(vs, reference_sequence) # Store allele as supremal
            except ValueError as e: # Fails for overlapping variants
                if "unorderable variants" in str(e):
                    # TODO how to handle duplicates? And how to handle multiple variants at same position?
                    warnings.warn(f"Could not parse sample {allele['alleleName']} due to double/overlapping variants")
                elif "empty" in str(e):
                    all_variants[allele["alleleName"]] = None # Will be disjoint with everything can leave out
                else:
                    raise e
    if cache_name: cache_set(all_variants, cache_name)
    return all_variants
