from modules.data import reference_get, pharmvar_get
import algebra as va

def parse_multi_hgvs(hgvs_lst, reference):
    """Wrapper for parsing a list of hgvs variants and combining them into a set of variants.
    """        
    variant_lst = []
    for hgvs in hgvs_lst:
        hgvs = fix_hgvs_position(hgvs)
        try:
            variant_lst += va.variants.parse_hgvs(hgvs, reference=reference) # Reference needed to handle insertions
        except: 
            raise Exception(f"HGVS string '{hgvs}' could not be parsed.")
        
    print(variant_lst)
    variant_set = set(variant_lst)
    if len(variant_set) != len(variant_lst):
        print(f"Double variants contained in {hgvs_lst}")
    return variant_set


def fix_hgvs_position(hgvs):
    """Fix HGVS notation problems with the position notation 
    
    Problems occur for deletions where a range is not specified but the (normally redundant) deleted area is.

    This notation is not HGVS compliant: 
    https://varnomen.hgvs.org/recommendations/DNA/variant/deletion/
    But it is present in the Pharmvar database
    """
    # TODO integrate change in va parser?
    if "del" in hgvs:
        area = hgvs.split("del")[1]
        position = hgvs.split("del")[0].split('.')[-1]
        if len(area) > 1 and '_' not in position:
            range = f"{position}_{int(position)+len(area)-1}"
            hgvs = hgvs.replace(position, range)
    return hgvs

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
    # All information associated with the allele
    coreallele = corealleles[coreallele_name]
    if coreallele_name != "CYP2D6*19":
        continue
    print(coreallele_name)
    # Get the variants contained in the sub- and core allele as HGVS notation
    coreallele_variants = [variant["position"] for variant in coreallele["variants"]]
    for suballele in suballeles:
        print(suballele["alleleName"])
        suballele_variants = [variant["position"] for variant in suballele["variants"]]
        # Convert HGVS notation to sequence notation
        s_coreallele_variants = parse_multi_hgvs(coreallele_variants, reference_sequence)
        s_suballele_variants = parse_multi_hgvs(suballele_variants, reference_sequence)
        # Find relation between core and suballeles
        try:
            relation = va.compare(reference_sequence, s_coreallele_variants, s_suballele_variants)
        except:
            raise ValueError(f"Could not compare variants {coreallele_variants} + {suballele_variants}")
        # Expect containment or equivalence
        if relation not in (va.Relation.EQUIVALENT, va.Relation.IS_CONTAINED):
            # print(coreallele["alleleName"], coreallele_variants)
            # print(suballele["alleleName"], suballele_variants)
            print(f"Unexpected relationship between core allele {coreallele['alleleName']} and suballele {suballele['alleleName']}: {relation}")

# QUESTION what is a variant group
# QUESTION why doesn't name match position
# QUESTION is an equivalence between a sub and core allele wrong

# TODO! fix error on 19.002

# TODO Check if core allele equivalent or contained in each suballele
# TODO check if HGVS name describes position field (not always the case)
# TODO Check if names follow a logical format
# TODO check if position is a valid HGVS string
# TODO check relations between star-alleles