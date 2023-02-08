from modules.data import reference_get, pharmvar_get
import algebra as va

def parse_multi_hgvs(hgvs_lst, reference):
    """Wrapper for parsing a list of hgvs variants and combining them into a set of variants.
    """        
    variant_lst = []
    for hgvs in hgvs_lst:
        hgvs = fix_hgvs_position(hgvs) # HGVS preprocessing
        try:
            variant_lst += va.variants.parse_hgvs(hgvs, reference=reference) # Reference needed to handle insertions
        except: 
            raise ValueError(f"HGVS string '{hgvs}' could not be parsed.")
        
    variant_set = set(variant_lst)
    if len(variant_set) != len(variant_lst):
        raise ValueError(f"Double variants contained in {hgvs_lst}")
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
            hgvs = hgvs.replace(position, range)
    return hgvs

def va_characterize_overlap(reference, lhs, rhs):
    """For a pair of variants characterize the elements of the overlap.
    
    This extends the variant algebra which detects but does not characterize the overlap.
    Overlap cannot always be explained by shared variant positions.
    """
    # TODO integrate change in va?
    # TODO make more efficient?
    # DOESNT WORK
    # Get all atomic operations
    # operations_lhs = [set(operations) for variant in lhs for operations in variant.atomics()]
    # operations_rhs = [set(operations) for variant in rhs for operations in variant.atomics()]
    # # Find membership
    # for var1 in operations_lhs:
    #     for var2 in operations_rhs:
    #         pass
    # shared = operations_lhs | operations_rhs # Operations in both sets
    # only_left = operations_lhs - shared # Operations only in left set
    # only_right = operations_rhs - shared # Operations only in right set
    # return only_left, shared, only_right


# Get the reference sequence relevant for the (current) gene of interest
reference_sequence = reference_get()

# List genes as symbols in Pharmvar
genes = pharmvar_get("genes/list") 
# All information associated with the (current) gene of interest
gene = pharmvar_get("genes/CYP2D6") 
# Group suballeles by core alleles and index by the star-allele notation of the core allele
corealleles = {allele["alleleName"]: allele for allele in gene["alleles"] if allele["alleleType"] == "Core"} 
suballeles = {coreallele: [sub_allele for sub_allele in gene["alleles"] if sub_allele["coreAllele"] == coreallele] for coreallele in corealleles.keys()}
for coreallele_name, suballeles in suballeles.items():
    if (coreallele_name != "CYP2D6*57"): # Test
        continue
    # All information associated with the allele
    coreallele = corealleles[coreallele_name]
    # Get the variants contained in the sub- and core allele as HGVS notation
    coreallele_variants = [variant["position"] for variant in coreallele["variants"]]
    for suballele in suballeles:
        suballele_variants = [variant["position"] for variant in suballele["variants"]]
        # Convert HGVS notation to sequence notation
        s_coreallele_variants = parse_multi_hgvs(coreallele_variants, reference_sequence)
        s_suballele_variants = parse_multi_hgvs(suballele_variants, reference_sequence)
        # Find relation between core and suballeles
        try:
            relation = va.compare(reference_sequence, s_coreallele_variants, s_suballele_variants)
        except Exception as e:
            if str(e) == "unorderable variants": # TODO fix these variants
                continue
            raise ValueError(f"Could not compare variants {coreallele_variants} + {suballele_variants}")
        # Expect containment or equivalence
        if relation not in (va.Relation.EQUIVALENT, va.Relation.IS_CONTAINED):
            # For unexpected relationships the overlap should be characterized
            # only_core, shared, only_sub = va_characterize_overlap(reference_sequence, s_coreallele_variants, s_suballele_variants)
            print(f"{coreallele['alleleName']}: Unexpected relationship {relation} with suballele {suballele['alleleName']}:")
            print(f"Core: {coreallele_variants}")
            print(f"sub: {suballele_variants}")

# QUESTION what is a variant group
# QUESTION why doesn't name match position
# QUESTION is an equivalence between a sub and core allele wrong

# TODO! Check if core allele equivalent or contained in each suballele
# TODO check if HGVS name describes position field (not always the case)
# TODO check if position is a valid HGVS string (not always the case)
# TODO Check if names follow a logical format
# TODO check relations between star-alleles
# TODO check if variants within suballele are disjoint (not always the case) and solve this