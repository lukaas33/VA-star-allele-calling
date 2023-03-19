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
        return
    matches.sort(key=sort_types) 
    # TODO find most specific classification (return core of sub?)
    # TODO do on pruned to avoid multiple
    core_matches = [match for match in matches if sort_types(match) == 1]
    reduced_matches = core_matches[:]
    for match in core_matches: # TODO use reduction methods here
        for _, right, relation in relations: # Find if this match is contained in other matches
            if relation != va.Relation.IS_CONTAINED:
                continue
            if right not in core_matches:
                continue
            reduced_matches.remove(match)
            break
    if len(core_matches) > 1:
        print(f"WARNING: multiple core matches: {core_matches}")
    return core_matches[0]

def star_allele_calling(sample, relations):
    """Determine star allele calling for a sample based on va relations."""
    specific_relations = [] # TODO use graph ds
    for relation in relations:
        if relation[0] == sample:
            specific_relations.append(relation)
        elif relation[1] == sample:
            specific_relations.append((relation[1], relation[0], relation[2]))
    # QUESTION: is it needed to look at suballeles for calling?
    # QUESTION: is it needed to look at individual variants for calling?
    # STEP 1: Trivial matching
    # Try to find a core/suballele or variant that is equivalent to the sample
    matches = []
    for _, allele, relation in specific_relations:
        if relation == va.Relation.EQUIVALENT:
            matches.append(allele)
    if len(matches) > 0:
        return find_most_specific(matches, specific_relations)
    # STEP 2: Matching by containment
    for _, allele, relation in specific_relations:
        if relation == va.Relation.CONTAINS:
            matches.append(allele)
    if len(matches) > 0:
        return find_most_specific(matches, specific_relations) 
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