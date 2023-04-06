import algebra as va
import networkx as nx
import warnings
from .data import api_get
import re

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
    # if len(matches["variants"]["personal"]) > 0:
        # warnings.warn(f"{sample}: classification may not be exact due to personal variants: {matches['variants']['personal']}.")
        # pass
    # if len(matches["variants"]["uncertain"]) > 0:
        # warnings.warn(f"{sample}: classification may not be exact due to extra variants which cannot be determined as noise: {matches['variants']['uncertain']}.")
        # pass
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
            "personal": [], 
            "noise": [], 
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
        if sort_types(variant) == 5: # Personal variant (not equal to any in the database)
            # TODO determine noise or uncertain based on sequence?
            matches["variants"]["personal"].append(variant)
        elif sort_types(variant) == 3: # Variant (equal to some in the database)
            if is_noise(variant, functions): # Noise, not relevant for calling
                matches["variants"]["noise"].append(variant)
            else: # May be relevant for calling
                matches["variants"]["uncertain"].append(variant)
        else:
            raise Exception(f"Unknown variant type: {variant}.")
    # Filter and return
    return find_best_match(sample, matches, functions) 

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

def is_silent(variant):
    """Determine if a variant could affect the protein sequence.
    
    Does this by checking if the variant equivalent descriptions.
    Equivalents are found by using the Mutalyzer API which account for different placements of the variant.
    
    If one of the equivalent descriptions affects the amino acid sequence, the variant is not silent.
    All coding representations must be in an UTR or intron to be considered non-coding.
    This is erring on the side of caution.
    """
    # TODO merge with test (REDO)
    # TODO what about splice variants?
    # Find equivalent representations of the variant
    data = api_get(f"https://mutalyzer.nl/api/normalize/{variant}") 
    # TODO make faster (do one call?/calculate this locally?)
    if "equivalent_descriptions" not in data: # TODO why is this?
        warnings.warn(f"Variant {variant} had no equivalent descriptions.")
        return False # Assume worst
    if 'c' not in data['equivalent_descriptions']: # Can ignore non-coding equivalents as these are always silent
        # n means non-coding, not in an ORF
        return True
    for nucleotide, protein in data["equivalent_descriptions"]['c']: # Check coding equivalents
        if '=' not in protein: # Not silent, affects aminoacids
            if classify_region(nucleotide) == "exon": # In exon
                # TODO should return False even if in intron?
                return False
    # TODO check for other representations?
    return True

def is_noise(variant, functions): 
    """Determine if a variant is noise.
    
    Noise is defined as not being relevant for calling.
    The variant is noise if it has no impact on the protein.
    Else it can be classified as 'uncertain'.

    This approach uses the pharmvar annotation but falls back on a sequence based approach.
    """    
    # Check the PharmVar impact annotation for the variant
    if variant not in functions:
        raise Exception(f"Variant {variant} had no impact annotation.")
    if functions[variant] == '': # Explicit no change
        # TODO confirm that this is the correct interpretation
        return True 
    elif functions[variant] == 'splice defect': # Change in expression
        return False
    elif functions[variant] == None: # Not known
        # TODO check meaning of None
        print(variant)
        if is_silent(variant): # TODO handle
            pass
        else:
            pass
    else: # Explicit change on protein level
        # TODO check non intronic
        # TODO check pattern
        return False

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
