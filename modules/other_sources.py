import re
import algebra as va
from easy_entrez import EntrezAPI
from easy_entrez.parsing import parse_dbsnp_variants
from .data import api_get
import warnings

entrez_api = EntrezAPI(
    'va-star-allele-calling',
    'lucas@vanosenbruggen.com',
    return_type='json'
    # TODO add api key for more requests (10/s instead of 3/s)
)

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

def is_silent_mutalyzer(variant):
    """Characterize a variant as silent or not.
     
    More specifically the variant will be classified as being in an exon or not.
    and as being synonymous or not. 
    If not in an exon it will always be considered synonymous.
    Will assume the worst case scenario which means that if the variant can be described as an exon it will be considered as such.

    Based on the Mutalyzer API which finds online annotations for the variant.
    """
    # WARNING this method is not reliable since it does not find all annotations (bug in mutalyzer!)
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

def find_id_hgvs(variant, reference):
    """Find the reference snp id of a variant based n hgvs."""
    # WARNING: not useful, is consistent with PharmVar and does not find a id when it is None
    chromosome = re.findall(r"NC_0*([0-9]*)\.", variant)[0]
    va_variant = va.variants.parse_hgvs(variant, reference=reference)
    position = f"{va_variant[0].start - 10}:{va_variant[0].end + 10}" # Larger since position of target must be entirely in range TODO smarter range 
    # Lookup ids at the same position
    result = entrez_api.search(
        {"chromosome": chromosome, "organism": 'human', "position": position},
        database='snp',
        max_results=1000 # Should not be limiting
    )
    # Check if any of these ids is the same as the variant
    ids = ["rs" + id for id in result.data['esearchresult']['idlist']]
    if len(ids) >= 1:
        result = entrez_api.fetch( # TODO make faster with single call
            ids, 
            database='snp',
            max_results=len(ids), 
        )  
        # Find id that matches HGVS 
        ids = []
        try:
            variants_data = parse_dbsnp_variants(result)
        except KeyError as e:
            warnings.warn(f"{variant} error: {e}") # TODO solve
            return None
        for id, hgvs_lst in variants_data.summary.HGVS.items(): 
            if variants_data.preferred_ids[id] != id: # Skip merged
                continue
            for other in hgvs_lst:
                if "NC" not in other or ".11" not in other: # Same reference TODO do nicer
                    continue
                # Compare with va since format may be different
                try:
                    va_other = va.variants.parse_hgvs(other, reference=reference)
                except ValueError as e:
                    warnings.warn(f"{other} could not be parsed: {e}")
                    continue
                relation = va.compare(reference, va_variant, va_other) # TODO fix this taking too long
                if relation == va.Relation.EQUIVALENT:
                    ids.append(id)
                    break
    if len(ids) > 1:
        warnings.warn(f"{variant} multiple ids found: {ids}")
        return None # TODO handle
    elif len(ids) == 0:
        warnings.warn(f"{variant} no ids found")
        return None # TODO handle
    return ids[0]
    


def is_silent_entrez(variant, ids):
    """Characterize a variant as silent or not.
    
    Similar to mutalyzer method but using the entrez API.
    """
    # TODO fix uid error that occurs sometimes
    classification = {"exon": False, "non-synonymous": False, "splicing": False}
    # QUESTION are these default values correct? (not if data is incomplete)
    # Find id of variant
    # TODO possible for unknown/personal variants?
    id = ids[variant]
    if id is None or id == '': # No information available, assume worst
        warnings.warn(f"{variant} has no id: '{id}'")
        classification['exon'] = True
        classification['non-synonymous'] = True
        classification['splicing'] = True
        # TODO is this correct?
        return classification
    # Do lookup on entrez
    result = entrez_api.fetch( # TODO make faster with single call
        [id], 
        database='snp',
        max_results=1, 
    )
    variants_data = parse_dbsnp_variants(result)
    #  Convert possible annotations to boolean
    consequences = list(variants_data.coordinates.consequence)
    if len(consequences) != 1: # Wrong
        raise Exception(f"Unexpected number of consequences: {consequences}")
    for consequence in consequences[0].split(','): 
        if consequence == "coding_sequence_variant": # Change in exon
            classification["exon"] = True
        elif consequence == "missense_variant": # Change in protein
            classification["exon"] = True 
            classification["non-synonymous"] = True
        elif consequence == "stop_gained": # Early stop mutation
            classification["non-synonymous"] = True
            classification["exon"] = True
        elif consequence == "frameshift_variant": # Frameshift (http://purl.obolibrary.org/obo/SO_0001589)
            classification["non-synonymous"] = True
            classification["exon"] = True
        elif consequence == "inframe_deletion": # Deletes amino acids (http://purl.obolibrary.org/obo/SO_0001822)
            classification["non-synonymous"] = True
            classification["exon"] = True
        elif consequence == "inframe_insertion": # Inserts amino acids (http://purl.obolibrary.org/obo/SO_0001821)
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
            raise Exception(f"Unknown consequence for {variant} ({id}): {consequence}")
    return classification