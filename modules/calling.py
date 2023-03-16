import algebra as va

def star_allele_calling(sample, relations):
    """Determine star allele calling for a sample based on va relations."""
    specific_relations = [relation for relation in relations if relation[0] == sample] # TODO use graph ds
    # QUESTION: is it needed to look at suballeles for calling?
    # QUESTION: is it needed to look at individual variants for calling?
    # STEP 1: Trivial matching
    # Try to find a core/suballele or variant that is equivalent to the sample
    matches = []
    for _, allele, relation in specific_relations:
        if relation == va.Relation.EQUIVALENT:
            matches.append(allele)
    # STEP _: Find the most specific classification (core, sub, variant)
    if len(matches) == 0:
        return
    matches.sort(key=lambda x: len(x)) # Can use length because of notation TODO watch out for alternative variant notation; find a better way
    # TODO find most specific classification (return core of sub?)
    return matches[0]

def print_classification(classifications):
    classified = 0 
    for sample, classification in classifications.items():
        if classification['A'] == classification['B'] == None:
            continue
        # TODO print None = *1?
        classified += 1
        print(f"{sample}: {classification['A']} {classification['B']}")
    total = len(classifications)
    print(f"Classified {classified} of {total}")