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
    """Find the best match from a list of matches.
    
    Also return a measure of certainty.
    """
    # TODO how to express certainty here (depend on variants?)
    if len(matches["equivalent"]) > 0: # Best match
        if len(matches["equivalent"]) > 1:
            raise Exception("This should not happen, multiple equivalent matches found.")
        return matches["equivalent"], 1
    elif len(matches["indirect"]) > 0: # Other matches
        # TODO is there a preference for direct matches or are indirect ones equally important?
        # TODO select best match
        return matches["indirect"], 0.5
    else: # Return default
        # QUESTION: is this valid
        return ["CYP2D6*1"], 0

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
    # Filter and return
    return find_best_match(matches) 

def print_classification(classifications):
    # TODO add detail level
    # TODO print certainty
    # TODO add simplification
    for sample, classification in classifications.items():
        print(f"{sample}:", end=' ')
        for key, _ in classification.items():
            classes, certainty = classification[key]
            print(",".join(classes), end=' ')
        print()