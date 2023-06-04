import re
import algebra as va
from easy_entrez import EntrezAPI
from easy_entrez.parsing import parse_dbsnp_variants
from .data import api_get, cache_get, cache_set
import warnings
from .calling import find_type, Type

# https://api.ncbi.nlm.nih.gov/variation/v0/
entrez_api = EntrezAPI(
    # TODO hide
    'va-star-allele-calling',
    'lucas@vanosenbruggen.com', 
    api_key="52181bc8fa51eacb5b70448a9e6fd6ae8209",
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
    # TODO use variation service instead of entrez? 
    # TODO make more general
    chromosome = re.findall(r"NC_0*([0-9]*)\.", variant)[0]
    va_variant = va.variants.parse_hgvs(variant) if reference is None else va.variants.parse_hgvs(variant, reference=reference)
    # TODO smarter range 
    position = f"{va_variant[0].start - 50}:{va_variant[0].end + 50}" # Larger since position of target must be entirely in range 
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

    Returns a set of ontology terms that describe the impact of the variant.
    An empty set is returned if no impact annotations are found.
    """
    # Find id of variant
    if id is None or id == '': # No information available
        raise Exception(f"{variant} no id found.")
    # Do lookup on entrez
    consequences = set()
    result = api_get("https://api.ncbi.nlm.nih.gov/variation/v0/refsnp/" + id[2:])
    for annotations in result['primary_snapshot_data']["allele_annotations"]:
        for annotation in annotations["assembly_annotation"]:
            # TODO why are there multiple allele/assembly annotations?
            for gene in annotation["genes"]:
                for rna in gene["rnas"]:
                    if rna["id"] != "NM_000106.6": # Transcript
                        # TODO don't hardcode references
                        # QUESTION are other transcripts relevant?
                        continue
                    if rna["product_id"] != "NP_000097.3":
                        # TODO don't hardcode references
                        # TODO needed after transcript check?
                        # QUESTION are other products relevant?
                        continue
                    for ontology in rna["sequence_ontology"]:
                        # TODO why would there be more
                        consequences.add(ontology["name"])
                    if "protein" not in rna.keys():
                        continue
                    for ontology in rna["protein"]["sequence_ontology"]:
                        # TODO why would there be more
                        consequences.add(ontology["name"])
    # TODO convert to pharmvar format?
    return consequences


def severity_GO(consequences):
    """Classify the GO effects of a variant by severity.
    
    GO terms:
    https://www.ncbi.nlm.nih.gov/variation/docs/glossary/

    Severity levels:
    0: unknown
    1: likely benign
    2: possibly harmful
    3: likely harmful
    TODO use enums
    """
    max_severity = 0
    for consequence in consequences:
        severity = 0
        if consequence == "coding_sequence_variant": # Change in exon 
            severity = 1
        elif consequence == "missense_variant": # Change in protein
            severity = 2
        elif consequence == "stop_gained": # Early stop mutation
            severity = 3
        elif consequence == "frameshift_variant": # Frameshift 
            severity = 3
        elif consequence == "inframe_deletion": # Deletes amino acids 
            severity = 2
        elif consequence == "inframe_insertion": # Inserts amino acids 
            severity = 2
        elif consequence == "splice_acceptor_variant": # Splice defect
            severity = 3
        elif consequence == "splice_donor_variant": # Splice defect
            severity = 3
        elif consequence == "synonymous_variant": # No impact on protein
            severity = 1
        elif consequence == "intron_variant": # Not in exon
            severity = 1
        elif "upstream" in consequence: # Outside ORF 
            severity = 1
        elif "downstream" in consequence: # Outside ORF 
            severity = 1
        elif "UTR_variant" in consequence: # In UTR 
            severity = 1
        else:
            # TODO handle other possible consequences    
            raise Exception(f"Unknown consequence {consequence}")
        max_severity = max(max_severity, severity)
    return max_severity

def severity_pharmvar(impact):
    """Classify the PharmVar effects of a variant by severity.

    HGVS notation:
    https://varnomen.hgvs.org/
    
    Severity levels:
    0: unknown
    1: likely benign
    2: possibly harmful
    3: likely harmful
    TODO use enums
    """
    if impact is None:
        # TODO should this be unknown?
        return 1
    if impact == "":
        return 0 
        # TODO should this be likely benign?
    if impact == "splice defect": # Splice defect
        return 3
    if re.match(r"[A-Z][0-9]*fs", impact): # Frameshift
       return 3
    if re.match(r"[A-Z][0-9]*X", impact): # Early stop
       return 3 
    if re.match(r"[A-Z][0-9]*([A-Z]|del|ins)", impact) or re.match(r"[0-9]*_[0-9]*[del|ins][A-Z]*", impact): # Missense
        return 2
    raise Exception(f"Unknown impact {impact}")

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
        if id is None:
            # TODO handle None
            raise Exception(f"{personal_variant} no id found.")
        # print(hgvs, id)
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


def impact_position(supremal):
    """Check if a variant is silent based on positions of variants.
    
    This is a different approach to the other methods since it doesn't rely on online sources.
    Therefore it can be used on novel variants.
    Instead it uses the position of a variant to see if it can influence the protein.
    This is based on the intron-exon borders.
    """
    raise DeprecationWarning("Will not be finished since annotations are used")
    # TODO do earlier for all personal variants (and more?)
    # TODO find exons dynamically for any gene
    # Exon borders (one-based; closed end) 
    # Source? https://www.ncbi.nlm.nih.gov/genome/gdv/browser/gene/?id=1565
    # QUESTION are there more exons sometimes?
    # QUESTION Does alternative splicing occur?
    # QUESTION are these fixed?
    exons = [ 
        (42126499, 42126752),
        (42126851, 42126992),
        (42127447, 42127634),
        (42127842, 42127983),
        (42128174, 42128350),
        (42128784, 42128944),
        (42129033, 42129185),
        (42129738, 42129909),
        (42130612, 42130810)
    ]    
    # Length of splice site
    # Part of intron that is used for recognition by spliceosome
    # Source: https://www.ncbi.nlm.nih.gov/gene/1565
    # TODO get dynamically
    # QUESTION is this correct?
    splice_sites = ("GT", "AG") 
    # Functions for checking exon influence
    # Check if a range is entirely in an exon
    in_exon = lambda s, e: any(
        start <= s and e <= end # Interval is contained in exon
        for start, end in exons)
    # Check if a range overlaps with an exon 
    overlap_exon = lambda s, e: any(
        max(s, start) <= min(e, end) # Interval overlaps with exon (highest begin is lowe than lowest end)
        for start, end in exons) 
    # Check if a range overlaps with a splice site
    overlap_splice_site = lambda s, e: any(
        s <= (start - len(splice_sites[0]) <= e) or # interval contains left splice site
        s <= (end + len(splice_sites[1])) <= e
        for start, end in exons)
    # Check if a mutation disturbs triplets by itself
    # difference (deleted or inserted) between area of influence and sequence is not a multiple of 3 
    frameshift = lambda v: (abs(v.end - v.start - len(v.sequence)) % 3) != 0 # TODO is this correct?
    # Check if supremal representation of variant (covers different placements) can influence the exon
    start = supremal.start + 1 # One-based position
    end = supremal.end # Closed end position
    if overlap_splice_site(start, end):
        return "possible splice defect" # QUESTION is this correct or can splice sites occur in a different way as well?
    if in_exon(start, end) or overlap_exon(start, end): 
        if frameshift(supremal): 
            return "possible frameshift" # QUESTION is this correct
        # TODO split further into synonymous, early stop and missense?
        return "possible missense" 
    return None # Intronic and thus certainly silent TODO correct value?

def relevance(sample, nodes, edges, functions, supremals, reference): 
    """Determine if extra variants may be relevant for calling."""
    raise DeprecationWarning("Not used")
    def overlap(a1, a2): max(a1.start, a2.start) <= min(a1.end, a2.end)
    # Get all called corealleles
    # TODO Only cores needed?
    alleles = star_allele_calling_all([sample], nodes, edges, functions, supremals, reference, detail_level=1)
    alleles = alleles[sample.split("_")[0]][sample.split("_")[1]]
    # Find extra variants: variants in the sample but not in any of the called alleles (pruned edges)
    # TODO do this on a graph?
    variants = []
    for edge in edges:
        if edge[0] == sample and find_type(edge[1]) in (Type.VAR, Type.P_VAR):
            variants.append(edge[1])
        elif edge[1] == sample and find_type(edge[0]) in (Type.VAR, Type.P_VAR):
            variants.append(edge[0])
    if len(variants) == 0: # No variants found
        return {}
    # Check the relevance of each variant
    variants_relevance = {}
    for variant in variants:
        # Check if variant possibly interferes with any allele (overlaps with supremal)
        interferes = any([overlap(supremals[variant], supremals[allele]) for allele in alleles if allele != "CYP2D6*1"])
        # Find the impact of the variant
        impact = functions[variant]
        if find_type(variant) == Type.P_VAR: 
            severity = severity_GO(impact)
        elif find_type(variant) == Type.VAR:
            severity = severity_pharmvar(impact)
        if severity == 2 and not interferes: # Change on protein level but doesn't interfere with any allele
            severity = 1
        variants_relevance[variant] = severity != 1 # Only not relevant if benign and doesn't interfere
        # TODO be less conservative with unknown impact?
    return variants_relevance
