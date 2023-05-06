import re
import algebra as va
from easy_entrez import EntrezAPI
from easy_entrez.parsing import parse_dbsnp_variants
from .data import api_get, cache_get, cache_set
import warnings

# TODO don't use easy_entrez, use entrez directly
# TODO use variation service instead of entrez? 
#       https://api.ncbi.nlm.nih.gov/variation/v0/
entrez_api = EntrezAPI(
    'va-star-allele-calling',
    'lucas@vanosenbruggen.com', # TODO hide
    api_key="52181bc8fa51eacb5b70448a9e6fd6ae8209", # TODO hide
    return_type='json'
)

def is_silent_mutalyzer(variant):
    """Characterize a variant as silent or not.
     
    More specifically the variant will be classified as being in an exon or not.
    and as being synonymous or not. 
    If not in an exon it will always be considered synonymous.
    Will assume the worst case scenario which means that if the variant can be described as an exon it will be considered as such.

    Based on the Mutalyzer API which finds online annotations for the variant.
    """
    raise DeprecationWarning("This method is not reliable since it does not find all annotations (bug in mutalyzer!")
    def classify_region(variant):
        """Classify region that a variant is in as UTR, intron or exon based on HGVS."""
        position = variant.split(':')[1].split('.')[1]
        if position[0] == '-': # Left from coding region
            return "5'UTR" # TODO use enums
        elif position[0] == '*':
            return "3'UTR"
        elif re.match(r"[0-9]{1,}[-+][0-9]{1,}", position):
            return "intron"
        return "exon" 
    classification = {'exon': False, 'non-synonymous': False, "splicing": False} # If not proven differently
    # Find equivalent representations of the variant
    data = api_get(f"https://mutalyzer.nl/api/normalize/{variant}") # TODO make faster with single call
    if "equivalent_descriptions" not in data.keys():
        warnings.warn(f"No equivalent descriptions found for {variant}.")
        # No annotations available, must assume worst case
        classification['exon'] = True
        classification['non-synonymous'] = True
        return classification
    data = data["equivalent_descriptions"]
    if any(((t not in {'c', 'n'}) for t in data.keys())): # Must be c or n or both
        raise Exception(f"Unhandled types: {set(data.keys()) }")
    if 'c' in data.keys():
        for nucleotide, protein in data['c']: # Check equivalent coding representations
            if classify_region(nucleotide) == 'exon': # Possibly in coding region
                classification['exon'] = True # Assume worst
                if '=' not in protein: # not synonymous and in exon
                    classification['non-synonymous'] = True
            # Don't need to check if synonymous since this only affects the protein sequence in exons
    else: 
        # n must be present (earlier check)
        # so only variants in non-coding area present
        pass
    classification['splicing'] = True # Must assume since no information otherwise
    # QUESTION how to detect splice/transcription factor variants?
    #           maybe by considering UTR and outside ORF differently?
    # QUESTION can assume that when some equivalent representations are known, all are known?
    return classification # All were intronic, outside ORF or UTR

def find_id_hgvs(variant, reference=None):
    """Find the reference snp id of a variant based n hgvs."""
    # TODO make more general
    chromosome = re.findall(r"NC_0*([0-9]*)\.g", variant)[0]
    va_variant = va.variants.parse_hgvs(variant) if reference is None else va.variants.parse_hgvs(variant, reference=reference)
    # TODO smarter range 
    position = f"{va_variant[0].start - 25}:{va_variant[0].end + 25}" # Larger since position of target must be entirely in range 
    # Lookup ids around the position
    result = entrez_api.search(
        {"chromosome": chromosome, "organism": 'human', "position": position},
        database='snp',
        max_results=1000 # Should not be limiting
    )
    # Check if any of these ids is the same as the variant
    ids = ["rs" + id for id in result.data['esearchresult']['idlist']]
    if len(ids) >= 1:
        # Check matches
        # TODO make faster with single call
        result = entrez_api.fetch(ids, database='snp', max_results=len(ids), 
        )  
        ids = []
        try:
            variants_data = parse_dbsnp_variants(result)
        except KeyError as e:
            warnings.warn(f"{variant} error: {e}. Trying again.") 
            # Seems to happen randomly, trying again should fix it.
            return find_id_hgvs(variant, reference) 
        # Find id that matches HGVS 
        for id, hgvs_lst in variants_data.summary.HGVS.items(): 
            if variants_data.preferred_ids[id] != id: # Skip merged snps
                continue
            # Check different representations 
            # TODO only check one?  
            for other in hgvs_lst:
                if other.split(':')[0] != variant.split(':')[0]: # Same reference
                    continue
                # Compare with va since format may be different
                try:
                    va_other = va.variants.parse_hgvs(other, reference=reference)
                except ValueError as e:
                    warnings.warn(f"{other} could not be parsed: {e}. Continuing with next hgvs.")
                    continue
                relation = va.compare(reference, va_variant, va_other)
                if relation == va.Relation.EQUIVALENT:
                    ids.append(id)
                    break
    if len(ids) > 1:
        raise Exception(f"{variant} multiple ids found: {ids}")
    elif len(ids) == 0:
        return None 
    return ids[0]

def get_annotation_entrez(variant, id):
    """Get impact of a variant based on entrez annotations.
    """
    # QUESTION are other transcripts and isoforms relevant?
    # Find id of variant
    if id is None or id == '': # No information available
        raise Exception(f"{variant} no id found.")
    # Do lookup on entrez
    # TODO make faster with single call
    result = entrez_api.fetch([id], database='snp', max_results=1)
    try:
        variants_data = parse_dbsnp_variants(result)
    except KeyError as e: 
        warnings.warn(f"{variant} error: {e}. Trying this again")
        # Try again, this error seems to occur randomly
        return get_annotation_entrez(variant, id) 
    for hgvs in variants_data.summary.HGVS[0]:
        # TODO don't hardcode references
        if hgvs.split(':')[0] == "NM_000106.6": # Transcript
            print(hgvs)
        elif hgvs.split(':')[0] == "NP_000097.3": # Protein
            print(hgvs)
    raise NotImplementedError("Not completely implemented yet")
    exit()
    # Convert possible annotations to boolean
    # TODO check only correct variant (NM_000106.6)
    consequences = list(variants_data.coordinates.consequence)
    if len(consequences) != 1:
        raise Exception(f"{variant} wrong number of consequences found: {consequences}")
    # TODO filter?
    # QUESTION what is None in this case?
    return consequences[0]


def entrez_consequence_binary(consequences):
    """Convert entrez consequence to binary.
    
    https://www.ncbi.nlm.nih.gov/variation/docs/glossary/
    """
    raise DeprecationWarning("Not used any more since it may be too simplistic")
    classification = {"exon": False, "non-synonymous": False, "splicing": False}
    for consequence in consequences[0].split(','): 
        if consequence == "coding_sequence_variant": # Change in exon
            classification["exon"] = True
        elif consequence == "missense_variant": # Change in protein
            classification["exon"] = True 
            classification["non-synonymous"] = True
        elif consequence == "stop_gained": # Early stop mutation
            classification["non-synonymous"] = True
            classification["exon"] = True
        elif consequence == "frameshift_variant": # Frameshift 
            classification["non-synonymous"] = True
            classification["exon"] = True
        elif consequence == "inframe_deletion": # Deletes amino acids 
            classification["non-synonymous"] = True
            classification["exon"] = True
        elif consequence == "inframe_insertion": # Inserts amino acids 
            classification["non-synonymous"] = True
            classification["exon"] = True
        elif consequence == "splice_acceptor_variant": # Splice defect
            classification["splicing"] = True
        elif consequence == "splice_donor_variant": # Splice defect
            classification["splicing"] = True
        elif consequence == "synonymous_variant": # No impact on protein
            pass
        elif consequence == "intron_variant": # Not in exon
            pass
        elif "upstream" in consequence: # Outside ORF (TODO check if this is correct)
            pass
        elif "downstream" in consequence: # Outside ORF (TODO check if this is correct)
            pass
        elif "UTR_variant" in consequence: # In UTR (TODO check if this is correct)
            pass
        else:
            # TODO handle other possible consequences    
            raise Exception(f"Unknown consequence {consequence}")
    return classification

def get_personal_ids(personal_variants, reference, cache_name=None):
    """Get ids of personal variants."""
    try: # Check if already calculated
        if cache_name is not None: return cache_get(cache_name)
    except:
        pass
    personal_ids = {}
    for personal_variant in personal_variants:
        hgvs = f"{reference['name']}:g.{personal_variant}"
        id = find_id_hgvs(hgvs, reference["sequence"])
        # TODO handle None
        personal_ids |= {personal_variant: id}
    if cache_name is not None: cache_set(personal_ids, cache_name)
    return personal_ids

def get_personal_impacts(personal_variants, ids, reference, cache_name=None):
    """Get Entrez impact of personal variants."""
    try: # Check if already calculated
        if cache_name is not None: return cache_get(cache_name)
    except:
        pass
    personal_impacts = {}
    for personal_variant in personal_variants:
        hgvs = f"{reference['name']}:g.{personal_variant}"
        id = ids[personal_variant]
        # TODO handle None
        impact = get_annotation_entrez(hgvs, id)
        personal_impacts |= {personal_variant: impact}
    if cache_name is not None: cache_set(personal_impacts, cache_name)
    return personal_impacts