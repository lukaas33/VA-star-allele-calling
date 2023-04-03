import algebra as va
import networkx as nx
import warnings

def sort_function(f):
    """Sort function annotation based on severity."""
    # QUESTION: what is the difference between uncertain and unknown?
    all_functions = ['normal function', 'unknown function', 'uncertain function', 'decreased function', 'no function']
    if f not in all_functions:
        raise Exception("Unknown function: " + f)
    functions = [all_functions[0]] + all_functions[2:]
    if f not in functions: 
        # TODO Where to put uncertain functions, now below normal
        #   but maybe they should be above normal?
        #   or maybe it should be considered as normal?
        #   or it can be handled differently
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
    if len(matches["equivalent"]) > 1: # The same as some allele
        raise Exception(f"{sample}: This should not happen, multiple equivalent matches found: {matches['equivalent']}.")
    sorted_matches.extend(matches["equivalent"])
    # Next step is the contained alleles
    # In the case of multiple contained corealleles, a choice is made based on functional annotation
    # *1 will never be in contained so it is not prioritized
    matches["contained"].sort(key=lambda c: sort_function(functions[c]), reverse=True)
    sorted_matches.extend(matches["contained"])
    core_matches = [m for m in sorted_matches if sort_types(m) == 1]
    # Check if answer is ambiguous
    # TODO implement choice here
    if len(core_matches) >= 2 and \
            sort_function(functions[core_matches[0]]) == sort_function(functions[core_matches[1]]):
        raise Exception(f"{sample}: At least two equally likely core matches found: {core_matches[0]} and {core_matches[1]}.")
    # Catch-all is the default allele
    # QUESTION is this valid?
    if len(core_matches) == 0:
        sorted_matches.append("CYP2D6*1")
    # Look at variants to determine the certainty of the calling
    # TODO how to express this to the user
    if len(matches["variants"]["personal"]) > 0:
        # warnings.warn(f"{sample}: classification may not be exact due to personal variants: {matches['variants']['personal']}.")
        pass
    if len(matches["variants"]["uncertain"]) > 0:
        warnings.warn(f"{sample}: classification may not be exact due to extra variants which cannot be determined as noise: {matches['variants']['uncertain']}.")
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
            if is_noise(variant, cont_graph): # Noise, not relevant for calling
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

def is_noise(variant, cont_graph): 
    """Determine if a variant is noise.
    
    Noise is determined as not being relevant for calling.
    """
    # TODO use pharmvar information on variant
    # Check based on occurrence in core alleles
    for connected in cont_graph[variant]:
        if sort_types(connected) == 1: # Is contained in a core allele
            return False
    return True

def print_classification(classifications, detail_level=0):
    """Print the classification of samples.
    
    Different detail levels are available.
    0: Only print best core match
    1: Simplify to core matches
    2: Print all matches
    """
    for sample, classification in classifications.items():
        selected_alleles = {'A': [], 'B': []}
        for phase, alleles in classification.items():
            for allele in alleles:
                if detail_level in (0, 1) and sort_types(allele) != 1: # Only core
                    continue
                selected_alleles[phase].append(allele)
                if detail_level == 0: # Only display best
                    break
        # Print
        print(f"{sample}: {','.join(selected_alleles['A'])}/{','.join(selected_alleles['B'])}")
