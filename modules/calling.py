import algebra as va

def star_allele_calling(sample, relations):
    """Determine star allele calling for a sample based on va relations."""
    specific_relations = [relation for relation in relations if relation[0] == sample] # TODO use graph ds
    # STEP 1: Trivial matching
    # Try to find a core allele that is equivalent to the sample
    for _, allele, relation in specific_relations:
        if relation == va.Relation.EQUIVALENT:
            print(f"{sample} is equivalent to {allele}")