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

def find_least_specific(matches, cont_graph):
    """Find least specific match (core) from a list of matches."""
    if len(matches) == 0:
        return None
    core_matches = [match for match in matches if sort_types(match) == 1]
    if len(core_matches) > 1: # TODO handle this case
        raise Exception(f"Multiple core matches found after reduction: {core_matches}")
    if len(core_matches) == 0: # Found no cores but may be able to find them from the suballeles
        # TODO find most specific classification (return core of sub?)
        indirect_matches = []
        for match in matches:
            for indirect_match, _ in cont_graph.in_edges(match): # Nodes contained in this one
                indirect_matches.append(indirect_match)
        return find_least_specific(indirect_matches, cont_graph) # Try again for indirect matches
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

    matches = list()
    # STEP 1: Trivial matching
    for match in eq_graph[sample]:
        matches.append(match)
    if len(matches) > 0:
        return find_least_specific(matches, cont_graph)
    return None

    # TODO do by traversing pruned graph
    # specific_relations = [] # TODO use graph ds
    # for left, right, relation in relations:
    #     if left == sample or right == sample:
    #         specific_relations.append((left, right, relation))

    # # Try to find a core/suballele or variant that is equivalent to the sample
    # matches = set() # TODO avoid duplicates in a different way
    # for left, right, relation in specific_relations:
    #     if relation == va.Relation.EQUIVALENT:
    #         matches.add(left if left != sample else right)
    # if len(matches) > 0:
    #     return find_most_specific(matches, relations)
    # # STEP 2: Matching by containment
    # for left, right, relation in specific_relations:
    #     if relation == va.Relation.IS_CONTAINED and right == sample:
    #         matches.add(left)
    # if len(matches) > 0:
    #     return find_most_specific(matches, relations) 
    # return None

def print_classification(classifications):
    classified = 0 
    for sample, classification in classifications.items():
        if classification['A'] == classification['B'] == None:
            continue
        # TODO how to represent uncertainty vs *1? (handle by classification?)
        classified += 1
        print(f"{sample}: {classification['A']} {classification['B']}")
    total = len(classifications)
    print(f"Classified {classified} of {total}")