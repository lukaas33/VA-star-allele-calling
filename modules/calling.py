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

def find_contained_cores(start, cont_graph, matches):
    """Find the first contained cores from a given start node."""
    for node, _ in cont_graph.in_edges(start):
        if sort_types(node) == 1: # Found core, don't go deeper (don't want core alleles contained in core alleles)
            matches.add(node)
        else: # Go deeper to find
            find_contained_cores(node, cont_graph, matches)

def find_least_specific(matches, cont_graph):
    """Find least specific match (core) from a list of matches."""
    if len(matches) == 0:
        return None
    core_matches = [match for match in matches if sort_types(match) == 1]
    if "CYP2D6*1" in core_matches and len(core_matches) > 1: # More than *1 
        core_matches.remove("CYP2D6*1") # Can be ignored # TODO reason?
    if len(core_matches) > 1: # TODO handle this case
        raise Exception(f"Multiple core matches found after reduction: {core_matches}")
    if len(core_matches) == 0: # Found no cores but may be able to find them from the suballeles
        # TODO always check deeply for corealleles? 
        indirect_matches = set()
        for match in matches:
            find_contained_cores(match, cont_graph, indirect_matches)
        return find_least_specific(list(indirect_matches), cont_graph) # Try again for indirect matches
    return core_matches[0] # Only match


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

    matches = set()
    # STEP 1: Trivial matching
    for match in eq_graph[sample]:
        matches.add(match)
    # QUESTION should equivalence be more important than containment?
    if len(matches) > 0: 
        return find_least_specific(matches, cont_graph)
    # STEP 2: Matching by containment
    # Traverse graph to find directly or indirectly contained cores
    find_contained_cores(sample, cont_graph, matches) 
    return find_least_specific(matches, cont_graph)

def print_classification(classifications):
    classified = 0 
    for sample, classification in classifications.items():
        if classification['A'] is None and classification['B'] is None:
            continue
        # TODO how to represent uncertainty vs *1? (handle by classification?)
        print(f"{sample}:", end=' ')
        for key, value in classification.items():
            if value is None:
                classification[key] = '?'
            else:
                classified += 1
            print(classification[key], end=' ')
        print()
    total = len(classifications) * 2
    print(f"Classified {classified} of {total}")