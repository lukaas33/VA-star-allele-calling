def image_reduction_transitive(relations):
    img_sel = ("CYP2D6*39", "CYP2D6*10", "CYP2D6*147")
    positions = ((0, 0), (0, 100), (100, 100))
    edges = set()
    for edge in relations:
        if edge[0] == edge[1]:
            continue
        if edge[0] not in img_sel or edge[1] not in img_sel:
            continue
        edges.add(edge)
    return img_sel, edges, positions

def image_reduction_symmetric(relations):
    img_sel = ("CYP2D6*115", "CYP2D6*109")
    positions = ((0, 0), (100, 0))
    edges = set()
    for edge in relations:
        if edge[0] == edge[1]:
            continue
        if edge[0] not in img_sel or edge[1] not in img_sel:
            continue
        edges.add(edge)
    return img_sel, edges, positions

def image_reduction_reflexive(relations):
    img_sel = ("CYP2D6*109",)
    positions = ((0, 0),)
    nodes = img_sel
    edges = set()
    for edge in relations:
        if edge[0] not in img_sel or edge[1] not in img_sel:
            continue
        edges.add(edge)
        break
    return nodes, edges, positions

def image_reduction_most_specific(relations):
    img_sel = ("CYP2D6*17", "CYP2D6*82", "CYP2D6*58")
    positions = ((0, 0), (100, 100), (100, 0))
    edges = set()
    for edge in relations:
        if edge[0] == edge[1]:
            continue
        if (edge[1], edge[0], edge[2]) in edges:
            continue
        if edge[0] not in img_sel or edge[1] not in img_sel:
            continue
        edges.add(edge)
    return img_sel, edges, positions

def image_reduction_common_ancestor(relations):
    img_sel = ("CYP2D6*9", "CYP2D6*109", "CYP2D6*115")
    positions = ((0, 0), (100, 100), (0, 100))
    edges = set()
    for edge in relations:
        if edge[0] == edge[1]:
            continue
        if (edge[1], edge[0], edge[2]) in edges:
            continue
        if edge[0] not in img_sel or edge[1] not in img_sel:
            continue
        edges.add(edge)
    return img_sel, edges, positions

def image_reduction_equivalence(relations):
    img_sel = ("CYP2D6*44", "CYP2D6*22", "CYP2D6*44.001")
    positions = ((0, 0), (100, 100), (100, 0))
    edges = set()
    for edge in relations:
        if edge[0] == edge[1]:
            continue
        if (edge[1], edge[0], edge[2]) in edges:
            continue
        if edge[0] not in img_sel or edge[1] not in img_sel:
            continue
        edges.add(edge)
    return img_sel, edges, positions