import algebra as va
import networkx as nx
import warnings
from .data import api_get
import re
from easy_entrez import EntrezAPI
from easy_entrez.parsing import parse_dbsnp_variants

entrez_api = EntrezAPI(
    'va-star-allele-calling',
    'lucas@vanosenbruggen.com',
    return_type='json'
    # TODO add api key for more requests (10/s instead of 3/s)
)

all_functions = ['unknown function', 'uncertain function', 'normal function', 'decreased function', 'no function']

def sort_function(f):
    """Sort function annotation based on severity."""
    # QUESTION: what is the difference between uncertain and unknown?
    if f not in all_functions:
        raise Exception("Unknown function: " + f)
    functions = all_functions[3:]
    if f not in functions: 
        return 0 
    return functions.index(f) + 1

def sort_types(v):
    """Sort variant types in order of specificity (somewhat arbitrary)."""
    if '*' in v:
        if '.' in v:
            return 2 # Suballele
        else:
            return 1 # Core allele
    elif v[:2] in ('HG', 'NA'):
        return 4 # Sample
    else:
        if ':' in v:
            return 3 # Variant
        return 5 # Personal variant
    
def find_contained_alleles(start, cont_graph, eq_graph, matches, visited):
    """Recursively the equivalent and contained sub- and corealleles from a given start node."""
    if start in visited: return # Already visited
    visited.add(start) # Needed to avoid equivalence loop and to avoid doubles
    if sort_types(start) in (1, 2): # Core or suballele
        matches.append(start)
    if sort_types(start) == 1: 
        return # Stop here
    # Find equivalent alleles
    if start in eq_graph.nodes():
        for match in eq_graph[start]: 
            find_contained_alleles(match, cont_graph, eq_graph, matches, visited) # Add equivalent alleles and maybe traverse
    # Find contained alleles
    if start in cont_graph.nodes():
        for match, _ in cont_graph.in_edges(start):
            find_contained_alleles(match, cont_graph, eq_graph, matches, visited) # Add contained alleles and maybe traverse

def find_best_match(sample, matches, functions):
    """Find the best matches from a list of matches.
    
    Also return a measure of certainty.
    """
    sorted_matches = []
    # Equivalent alleles are the most certain
    if len(matches["equivalent"]) > 0:
        if len(matches["equivalent"]) > 1: # The same as some allele
            raise Exception(f"{sample}: This should not happen, multiple equivalent matches found: {matches['equivalent']}.")
        sorted_matches.append(matches["equivalent"])
    # Next step is the contained alleles
    if len(matches["contained"]) > 0:
        # In the case of multiple contained corealleles, a choice is made based on functional annotation
        annotated = {f: [m for m in matches["contained"] if functions[m] == f] for f in all_functions}
        # If the annotation is no function, this allele is most disruptive and is chosen
        if len(annotated["no function"]) > 0:	
            sorted_matches.append(annotated["no function"])
        # If uncertain or unknown functions are present, the ordering after no function cannot be determined
        # An unknown function could be no function and thus no allele should be preferred.
        if len(annotated["unknown function"]) > 0 or len(annotated["uncertain function"]) > 0:
            # Put all in the same position of ordering
            sorted_matches.append(
                annotated["decreased function"] +
                annotated["normal function"] +
                annotated["unknown function"] +
                annotated["uncertain function"]
            )
        else: # Ordering can be determined between normal and decreased function
            # More disruptive alleles are chosen first
            if len(annotated["decreased function"]) > 0:
                sorted_matches.append(annotated["decreased function"])
            if len(annotated["normal function"]) > 0:
                sorted_matches.append(annotated["normal function"])

    core_matches = [[m for m in match if sort_types(m) == 1] for match in sorted_matches]
    if sum([len(m) for m in core_matches]) == 0: # No coreallele found
        # Catch-all is the default allele
        # *1 will never be in contained so it is not prioritized over others
        sorted_matches.append(["CYP2D6*1"]) # QUESTION is this valid?
    else: 
        # Check if answer is ambiguous (has multiple corealleles)
        for core_match in core_matches:
            if len(core_match) == 0: # Skip
                continue
            if len(core_match) > 1: # Ambiguous
                # TODO return multiple equally likely answers
                raise Exception(f"{sample}: This should not happen, multiple corealleles found: {core_match}.")
            else: # No ambiguity for best answer
                break
    # Look at variants to determine the certainty of the calling
    # TODO how to express this to the user
    n_protein = len(matches['variants']['protein'])
    n_uncertain = len(matches['variants']['uncertain'])
    if n_protein + n_uncertain > 0:
        # warnings.warn(f"{sample}: classification of {sorted_matches[0]} is not certain due to {n_protein} variants that may affect protein function ({matches['variants']['protein']}) and {n_uncertain} variants with uncertain effect on protein function ({matches['variants']['uncertain']}).")
        pass
    return sorted_matches


def star_allele_calling(sample, nodes, edges, functions):
    """Determine star allele calling for a sample based on va relations.
    
    Find based on pruned graph containing relations between samples and alleles and between themselves.
    """
    # QUESTION: is it needed to look at suballeles for calling?
    # QUESTION: is it needed to look at individual variants for calling?
    eq_graph = nx.Graph()
    eq_graph.add_nodes_from(nodes)
    eq_graph.add_edges_from([(left, right) for left, right, relation in edges if relation == va.Relation.EQUIVALENT])
    cont_graph = nx.DiGraph()
    cont_graph.add_nodes_from(nodes)
    cont_graph.add_edges_from([(left, right) for left, right, relation in edges if relation == va.Relation.IS_CONTAINED])

    # All information of a calling
    matches = {
        "equivalent": [], 
        "contained": [], 
        "variants": {
            "protein": [], 
            "uncertain": []
        }
    }
    # Find equivalent alleles
    contained = []
    if sample in eq_graph.nodes(): matches["equivalent"] = [match for match in eq_graph[sample] if sort_types(match) != 5]
    # Find contained alleles
    find_contained_alleles(sample, cont_graph, eq_graph, contained, set())
    # Check if any contain another, keep most specific
    # TODO add to initial getting?
    for node1 in set(contained):
        for node2 in contained:
            if node1 == node2: # Skip self 
                continue
            if sort_types(node1) == 1 and sort_types(node2) == 2: # Skip containment of core in suballeles as this is expected
                continue
            if node1 in nx.ancestors(cont_graph, node2): # Node1 is contained in node2
                break # Skip this allele
        else: # No issues 
            if node1 in matches["equivalent"]: # Skip equivalent alleles
                continue
            matches["contained"].append(node1) # Keep match
    # Split variants into types
    variants = [match for match, _ in cont_graph.in_edges(sample) if sort_types(match) in (3, 5)]
    for variant in variants:
        complete_variant = variant
        if sort_types(variant) == 5: # Personal variant (not equal to any in the database)
            complete_variant = "NC_000022.11:g." + variant
        if is_noise(complete_variant, functions): # No evidence of relevance
            matches["variants"]["uncertain"].append(variant)
        else: # Evidence of protein impact
            matches["variants"]["protein"].append(variant)
    # Filter and return
    return find_best_match(sample, matches, functions) 

def star_allele_calling_all(samples, nodes, edges, functions):
    classifications = {sample[:-1]: {'A': None, 'B': None} for sample in sorted(samples)} 
    for sample in samples:
        sample_source, phasing = sample[:-1], sample[-1]
        classification = star_allele_calling(sample, nodes, edges, functions)
        classifications[sample_source][phasing] = classification
    return classifications

def matches_core(match):
    """Print the core allele from a match."""
    if sort_types(match) == 1:
        return match
    elif sort_types(match) == 2:
        return match[:-4] # QUESTION is this the best way to get the core?
    else:
        raise Exception(f"Unexpected match type: {match}")
        # QUESTION needed to also handle variants?

def classify_region(variant):
    """Classify region that a variant is in as UTR, intron or exon."""
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

def is_silent_position(variant, allele):
    """Check if a variant is silent based on positions of variants.
    
    This is a different approach to the other methods since it doesn't rely on online sources.
    Instead it uses the position of a variant to see if it can undo the effect of a different variant.
    """
    pass

# Check if a string is a protein mutation
protein_mutation = lambda s: re.match(r"([A-Z][0-9]{1,}([A-Z]|fs|del)|[0-9]{1,}_[0-9]{1,}(ins|dup)[A-Z]{1,}(x2){0,})", s) # TODO change to match all possible values instead of observed

def is_noise(variant, functions): 
    """Determine if a variant is noise.
    
    Noise is defined as not being relevant for calling.
    The variant is noise if it has no impact on the protein.

    This approach uses the online annotations but falls back on a sequence based approach.
    """
    if variant in functions: # PharmVar variant
        function = functions[variant] 
    else: # Personal variant
        function = None # No annotation known
    # Check the PharmVar impact annotation for the variant
    if function == '': # Explicit no change
        # TODO handle exceptions here
        return True 
    elif function == 'splice defect': # Change in expression
        return False
    elif function == None: # Not known in PharmVar
        # TODO handle this
        pass
    elif protein_mutation(function): # Explicit change on protein level
        return False
    else:
        raise Exception(f"Unhandled function: {function}")

    return False # Default is not noise

def print_classification(classifications, detail_level=0):
    """Print the classification of samples.
    
    Different detail levels are available.
    0: Only print best core match
    1: Simplify to core matches
    2: Print all matches
    """
    for sample, classification in classifications.items():
        selected_alleles = {'A': [], 'B': []}
        for phase, ordered_alleles in classification.items():
            for alleles in ordered_alleles:
                for allele in alleles: # TODO do this nicer
                    if detail_level in (0, 1) and sort_types(allele) != 1: # Only core
                        continue
                    selected_alleles[phase].append(allele)
                    if detail_level == 0: # Only display best
                        break
                else: # No break 
                    continue
                break
        # Print
        print(f"{sample}: {','.join(selected_alleles['A'])}/{','.join(selected_alleles['B'])}")
