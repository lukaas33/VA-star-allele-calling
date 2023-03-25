import algebra as va
import networkx as nx

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
        return 3 # Variant
    
def find_contained_alleles(start, cont_graph, matches):
    """Find the first contained cores from a given start node."""
    for node, _ in cont_graph.in_edges(start):
        if sort_types(node) in (1, 2): # Found core, don't go deeper (don't want core alleles contained in core alleles)
            # Check if this core is contained in another (possible via other path)
            for other in matches:
                if node in nx.ancestors(cont_graph, other):
                    return
            matches.append(node)
        else: # Go deeper to find others
            find_contained_alleles(node, cont_graph, matches)

def find_best_match(matches):
    """Find the best matches from a list of matches.
    
    Also return a measure of certainty.
    """
    # TODO how to express certainty here (depend on extra variants?)
    best_matches = []
    if len(matches["equivalent"]) > 0: # Best match
        if len(matches["equivalent"]) > 1:
            raise Exception("This should not happen, multiple equivalent matches found.")
        best_matches.append((matches["equivalent"][0], 10))
    elif len(matches["indirect"]) > 0: # Other matches
        # TODO select best match from multiple
        for m in matches["indirect"]:
            certainty = 5
            # Lower likelihood of *1 because other matches are more likely
            # TODO is this valid?   
            if matches_core(m) == "CYP2D6*1":
                certainty -= 1
            best_matches.append((m, certainty))
    else: # Return default
        # QUESTION: is this valid
        best_matches.append(("CYP2D6*1", 0))

    # Also add core alleles
    extended_best_matches = best_matches[:]
    for match, certainty in best_matches:
        core = matches_core(match)
        # Less certain because of extra variants
        certainty -= 1
        if core not in [m[0] for m in extended_best_matches]:
            extended_best_matches.append((core, certainty)) 
    # TODO check for equally good matches
    return extended_best_matches

def star_allele_calling(sample, nodes, edges):
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

    matches = {"equivalent": [], "contained": [], "variants": [], "indirect": []}
    # STEP 1: Trivial matching
    if sample in eq_graph.nodes():
        matches["equivalent"] = [match for match in eq_graph[sample]]
    # STEP 2: Matching by containment
    if sample in cont_graph.nodes():
        # TODO split personal variants
        for match, _ in cont_graph.in_edges(sample):
            if sort_types(match) == 3:
                matches["variants"].append(match)
            elif sort_types(match) in (1, 2):
                matches["contained"].append(match)
            else:
                raise Exception(f"Unexpected match type: {match}")
        # STEP 3: matching by indirect containment (more detail)    
        find_contained_alleles(sample, cont_graph, matches["indirect"])
        # TODO also find indirectly contained cores (from suballeles)
    # Filter and return
    return find_best_match(matches) 

def matches_core(match):
    """Print the core allele from a match."""
    if sort_types(match) == 1:
        return match
    elif sort_types(match) == 2:
        return match[:-4] # QUESTION is this the best way to get the core?
    else:
        raise Exception(f"Unexpected match type: {match}")
        # QUESTION needed to also handle variants?

def print_classification(classifications, detail_level=1):
    """Print the classification of samples.
    
    Different detail levels are available.
    0: Only print best match based on certainty
    1: Simplify to core matches
    2: Print all matches
    """
    # TODO move matches core to find best match function
    for sample, classification in classifications.items():
        for key, _ in classification.items():
            staralleles = classification[key]
            staralleles.sort(key=lambda c: c[1], reverse=True)
            print(f"{sample}{key}:", end='\n')
            for starallele, certainty in staralleles:
                if detail_level <= 1: # Simplify to core matches
                    if sort_types(starallele) == 2:
                        continue
                elif detail_level == 2: # Print all matches
                    pass
                print(f"{starallele} ({certainty})", end='\n')
                if detail_level == 0: # Only print best match based on certainty
                    break
        print()