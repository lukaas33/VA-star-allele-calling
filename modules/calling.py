import algebra as va
import networkx as nx

def sort_function(f):
    """Sort function annotation based on severity."""
    # QUESTION: what is the difference between uncertain and unknown?
    all_functions = ['normal function', 'unknown function', 'uncertain function', 'decreased function', 'no function']
    if f not in all_functions:
        raise Exception("Unknown function: " + f)
    functions = ['normal function', 'decreased function', 'no function']
    if f not in functions: 
        return 0 # TODO Where to put uncertain functions?
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
    
def find_contained_alleles(start, cont_graph, eq_graph, matches, visited=set()):
    """Recursively the equivalent and contained sub- and corealleles from a given start node."""
    if start in visited: return # Already visited
    visited.add(start) # Needed to avoid equivalence loop and to avoid doubles
    if sort_types(start) in (1, 2): # Core or suballele
        matches.append(start)
    if sort_types(start) == 1: 
        return # Stop here

    # Find equivalent alleles
    if start in eq_graph.nodes():
        equivalent = list(eq_graph[start])
        if len(equivalent) > 1: # check for multiple equal connections
            # Should only exist for cores which are not expanded
            # raise Exception(f"This should not happen, multiple equivalent alleles of {start} found: {equivalent}."
            pass
        for match in equivalent: 
            find_contained_alleles(match, cont_graph, eq_graph, matches) # Add equivalent alleles and maybe traverse
    # Find contained alleles
    if start in cont_graph.nodes():
        for match, _ in cont_graph.in_edges(start):
            find_contained_alleles(match, cont_graph, eq_graph, matches)

    # for node, _ in cont_graph.in_edges(start):
    #     if sort_types(node) == 2: # Found suballele
    #         matches.append(node)
    #     elif sort_types(node) == 1: # Core allele
    #         # Check if this allele is contained in another 
    #         for other in matches:
    #             if other == node: # Skip self
    #                 continue
    #             if sort_types(other) == 2: # Skip suballeles as this is expected
    #                 continue
    #             if node in nx.ancestors(cont_graph, other): # Is contained in one of the matches
    #                 return # Skip this allele
    #         matches.append(node)
    #         return # Don't traverse further since this doesn't add any information
    #     # Go deeper to find others
    #     find_contained_alleles(node, cont_graph, matches)
    # return

def find_best_match(matches, functions):
    """Find the best matches from a list of matches.
    
    Also return a measure of certainty.
    """
    # TODO how to express certainty (depend on extra variants?)
    best_matches = []
    if len(matches["equivalent"]) > 0: # The same as some allele 
        if len(matches["equivalent"]) > 1:
            raise Exception("This should not happen, multiple equivalent matches found.")
        # TODO merge this case (since extra variants will change the certainty)
        best_matches.append((matches["equivalent"][0], 10))
    if len(matches["indirect"]) > 0: # Contains some allele
        # TODO select best match from multiple
        for m in matches["indirect"]:
            certainty = 5
            # CYP2D6*1 is the default, prioritize other alleles
            # TODO is this valid?
            if matches_core(m) == "CYP2D6*1":
                best_matches.append((m, certainty))
                continue
            # Function influences certainty 
            # deleterious mutations are less likely to be reversed
            # QUESTION is this valid
            certainty += (sort_function(functions[m]) + 1) / 4
            best_matches.append((m, certainty))
    if len(best_matches) == 0: # Return default
        # QUESTION: is this valid
        best_matches.append(("CYP2D6*1", 0))
    # # Also add core alleles
    # extended_best_matches = best_matches[:]
    # for match, certainty in best_matches:
    #     core = matches_core(match)
    #     # Less certain because of extra variants
    #     if core not in [m[0] for m in extended_best_matches]:
    #         extended_best_matches.append((core, certainty)) 
    # Sort based on certainty 
    best_matches.sort(key=lambda c: c[1], reverse=True)
    # TODO remove 'certainty' 
    return best_matches

index = 7
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

    matches = {"equivalent": [], "contained": [], "variants": []}
    # Find equivalent alleles
    contained = []
    matches["equivalent"] = [match for match in eq_graph[sample] if sort_types(match) != 5]
    # Find contained alleles
    find_contained_alleles(sample, cont_graph, eq_graph, contained)
    # Check if any contain another, keep most specific
    for i, node1 in enumerate(contained):
        for node2 in contained[i+1:]:
            if node1 == node2: # Double 
                break # Don't keep
            if sort_types(node1) == 1 and sort_types(node2) == 2: # Skip containment suballeles as this is expected
                continue
            if node1 in nx.ancestors(cont_graph, node2): # Node1 is contained in node2
                break # Skip this allele
        else: # No issues 
            if node1 in matches["equivalent"]: # Skip equivalent alleles
                continue
            matches["contained"].append(node1) # Keep match
    print(sample)
    print(matches)
    global index
    index -= 1
    if index == 0:
        exit()
    # Find extra variants
    matches["variants"] = [match for match, _ in cont_graph.in_edges(sample) if sort_types(match) in (3, 5)]

    # # STEP 1: Trivial matching
    # if sample in eq_graph.nodes():
    #     matches["equivalent"] = [match for match in eq_graph[sample] if sort_types(match) != 5]
    # # # STEP 2: Matching by containment
    # # if sample in cont_graph.nodes():
    # #     # TODO split personal variants
    # #     for match, _ in cont_graph.in_edges(sample):
    # #         if sort_types(match) in (3, 5):
    # #             matches["variants"].append(match)
    # #         elif sort_types(match) in (1, 2):
    # #             matches["contained"].append(match) # TODO can leave out?
    # #         else:
    # #             raise Exception(f"Unexpected match type: {match}")
    # # STEP 3: matching by indirect containment (more detail)
    # for variant in [sample] + matches["equivalent"]: # Also find contained in equivalent alleles
    #     if variant in cont_graph.nodes():
    #         find_contained_alleles(variant, cont_graph, matches["indirect"])
    # # Filter and return
    # return find_best_match(matches, functions) 

def matches_core(match):
    """Print the core allele from a match."""
    if sort_types(match) == 1:
        return match
    elif sort_types(match) == 2:
        return match[:-4] # QUESTION is this the best way to get the core?
    else:
        raise Exception(f"Unexpected match type: {match}")
        # QUESTION needed to also handle variants?

def print_classification(classifications, detail_level=0):
    """Print the classification of samples.
    
    Different detail levels are available.
    0: Only print best match without certainty score
    1: Simplify to core matches
    2: Print all matches
    """
    for sample, classification in classifications.items():
        if detail_level == 0:
            print(f"{sample}:", end=' ')
        for i, key in enumerate(classification.keys()):
            staralleles = classification[key]
            if detail_level != 0:
                print(f"{sample}{key}:", end='\n')
            for starallele, certainty in staralleles:
                if detail_level in (0, 1): # Simplify to core matches
                    if sort_types(starallele) == 2:
                        continue
                elif detail_level == 2: # Print all matches
                    pass
                if detail_level == 0: # Only print best match based on certainty
                    print(f"{starallele}", end=('/' if i == 0 else ''))
                    # Check if there are equally likely matches (cannot print one in that case)
                    optimal = [c for a, c in staralleles if c == certainty and sort_types(a) == 1]
                    if len(optimal) > 1:
                        raise Exception(f"This should not happen, multiple equally likely core alleles found: {optimal}")
                    break
                print(f"{starallele} ({round(certainty, 2)})", end='\n')
        print()