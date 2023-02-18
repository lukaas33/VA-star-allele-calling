import algebra as va
from .graph import prune_relations

def count_relations(relations):
    print(f"{len(relations)} relations:")
    counts = {relationType: 0 for relationType in va.Relation}
    for relation in relations:
        counts[relation[2]] += 1
    for relationType, count in counts.items():
        print(f"  {count} {relationType}")

def relation_statistics(nodes, relations):
    count_relations(relations)
    print("After pruning:")
    count_relations(prune_relations(nodes, relations))