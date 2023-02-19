import algebra as va

def count_relations(relations):
    counts = {relationType: 0 for relationType in va.Relation}
    for relation in relations:
        counts[relation[2]] += 1
    return counts

def count_arity(nodes, relations):
    arities = {node: {relationType: 0 for relationType in va.Relation} for node in nodes}
    for relation in relations:
        arities[relation[0]][relation[2]] += 1
    return arities
