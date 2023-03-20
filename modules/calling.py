import algebra as va


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

def find_most_specific(matches, relations):
    """Find most specific match from a list of matches."""
    if len(matches) == 0:
        return None
    # TODO find most specific classification (return core of sub?)
    core_matches = [match for match in matches if sort_types(match) == 1]
    if len(core_matches) > 1:
        reduced_matches = core_matches[:]
        for match in core_matches: # TODO use reduction methods here
            for left, right, relation in relations: # Find if this match is contained in other matches
                if left != match:
                    continue
                if right not in core_matches:
                    continue
                if relation != va.Relation.IS_CONTAINED:
                    continue
                if left in reduced_matches:
                    reduced_matches.remove(left)
        if len(reduced_matches) > 1:
            raise Exception(f"Multiple core matches found after reduction: {reduced_matches}")
        if len(reduced_matches) == 0:
            raise Exception("No core matches found after reduction")
        return reduced_matches[0]
    if len(core_matches) == 0:
        raise Exception("No core matches found")
    return core_matches[0]

def star_allele_calling(sample, relations):
    """Determine star allele calling for a sample based on va relations."""
    # TODO do by traversing pruned graph
    specific_relations = [] # TODO use graph ds
    for left, right, relation in relations:
        if left == sample or right == sample:
            specific_relations.append((left, right, relation))

    # QUESTION: is it needed to look at suballeles for calling?
    # QUESTION: is it needed to look at individual variants for calling?
    # STEP 1: Trivial matching
    # Try to find a core/suballele or variant that is equivalent to the sample
    matches = set() # TODO avoid duplicates in a different way
    for left, right, relation in specific_relations:
        if relation == va.Relation.EQUIVALENT:
            matches.add(left if left != sample else right)
    if len(matches) > 0:
        return find_most_specific(matches, relations)
    # STEP 2: Matching by containment
    for left, right, relation in specific_relations:
        if relation == va.Relation.IS_CONTAINED and right == sample:
            matches.add(left)
    if len(matches) > 0:
        return find_most_specific(matches, relations) 
    return None

def print_classification(classifications):
    classified = 0 
    for sample, classification in classifications.items():
        if classification['A'] == classification['B'] == None:
            continue
        # TODO print None = *1?
        # TODO how to represent uncertainty vs *1?
        classified += 1
        print(f"{sample}: {classification['A']} {classification['B']}")
    total = len(classifications)
    print(f"Classified {classified} of {total}")