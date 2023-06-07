import algebra as va
import warnings
from .data import cache_get, cache_set
import os
import vcf
from .utils import normalise
from modules.calling import find_type, Type

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


def parse_samples(directory, reference, phased=False, cache_name=None):
    """ Parse sample VCF files as variant objects."""
    # TODO detect if phased
    # TODO handle mixed phasing
    if cache_name:
        try:
            return cache_get(cache_name)
        except:
            pass
    samples = {}
    for filename in os.listdir(directory):
        sample_name = filename.split('.')[0]
        if phased:
            allele = {"A": {}, "B": {}} # First and second allele
        else:
            allele = {"hom": {}, "het": {}, "all": {}} # Homozygous, heterozygous and all variants
        with open(os.path.join(directory, filename), 'r') as file:
            reader = vcf.Reader(file)
            # QUESTION: what are the filter, quality, format, info fields?
            for record in reader:
                # Validate if reference of vcf files is on the reference sequence
                if record.REF != reference["sequence"][record.start:record.end]:
                    raise ValueError("Reference sequence does not match")
                phasing = record.samples[0].data[0]
                # Alternative variants known to be heterozygous are stored in the ALT field
                if not phased and phasing == '1/2': 
                    # TODO handle
                    continue
                # Check ALT field. Can have multiple values for different phases, etc.
                if len(record.ALT) > 1: 
                    raise ValueError("Multiple ALT alleles not supported") # TODO handle different alt values?
                # Check multiple samples
                if len(record.samples) > 1:
                    raise ValueError("Multiple samples not supported") # TODO handle
                # Create variant with Zero based half-open positions 
                # TODO parse insT and >T variants differently?
                variant = va.Variant(record.start, record.end, record.ALT[0].sequence)
                hgvs = normalise(va.variants.to_hgvs([variant])).split(':')[1].split('.')[1]
                # Store variant in correct alleles
                if phased:
                    if '|' not in phasing: 
                        raise ValueError("Sample is not phased")
                    # Add to correct phased allele
                    # 0|1 or 1|0 shows that the variant is one of the two alleles, and can be used to group variants into two alleles
                    # 1|1 
                    for p, phase in zip("AB", phasing.split('|')):
                        if phase == '1':
                            allele[p][hgvs] = variant
                else:
                    if '|' in phasing:
                        raise ValueError("Sample is phased")
                    # Add to unphased allele
                    # 1/1 is homozygous, here perfect phasing information is available
                    # 1/0 or 0/1 is heterozygous and no phasing information is available
                    allele["all"][hgvs] = variant
                    if phasing == "1/1":
                        allele["hom"][hgvs] = variant
                    elif phasing == "1/0" or phasing == "0/1": # TODO handle 1/0 differently?
                        allele["het"][hgvs] = variant
                    else:
                        raise ValueError("Unknown phasing", phasing)   
            if not phased: del allele["het"] # Remove as not used
            # Add alleles to sample
            for p in allele:
                samples[f"{sample_name}_{p}"] = allele[p]
    if cache_name: cache_set(samples, cache_name)
    return samples

def samples_to_supremal(samples_phased, samples_unphased, reference, supremal_extended, cache_name=None):
    if cache_name:
        try:
            return cache_get(cache_name)
        except:
            pass
    supremal_samples = {}
    homozygous = {sample.split('_')[0]: set(samples_unphased[sample].keys()) for sample in samples_unphased if sample.split('_')[1] == 'hom'}
    for sample, variants in (samples_phased | samples_unphased).items():
        if variants == {}: # No supremal for empty, will be disjoint with everything, can be ignored
            supremal_samples[sample] = None
            continue
        # Parse alleles (variants together)
        try:
            supremal_samples[sample] = to_supremal(list(variants.values()), reference["sequence"]) # Try to find supremal for sample
        except ValueError as e:
            # TODO can assume that overlapping variants are in different phases?
            # TODO use different placement for the overlapping variants?
            if "unorderable variants" in str(e):
                warnings.warn(f"Could not parse sample {sample} due to double/overlapping variants")
                print(*variants.keys())
            else:
                raise e
        # Parse individual variants
        for hgvs, variant in list(variants.items()): # Find supremal for individual variants
            if hgvs in supremal_samples:
                continue
            supremal_v = to_supremal([variant], reference["sequence"] )
            # Filter out variants that are in the database (not personal)
            # Relations with these variants will show up later
            for v in supremal_extended:
                if find_type(v) != Type.VAR:
                    continue
                rel = va.relations.supremal_based.compare(reference["sequence"], supremal_v, supremal_extended[v])
                if rel == va.Relation.EQUIVALENT: # Not personal
                    del variants[hgvs] # Can delete since already parsed as allele
                    # Change homozygous hgvs
                    for _, hs in homozygous.items():
                        if hgvs in hs:
                            hs -= {hgvs,}
                            hs |= {v,}
                    break
            else: # Personal variant is saved individually
                supremal_samples[hgvs] = supremal_v
    if cache_name: cache_set((supremal_samples, homozygous), "supremal_samples")
    return supremal_samples, homozygous