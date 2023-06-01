# TODO move to scripts
import algebra as va 

image_configs = {}

image_configs["transitive"] = {
    "selection": ("CYP2D6*39", "CYP2D6*10", "CYP2D6*147"),
    "positions": ((0, 0), (0, 100), (100, 100)),
    "edges": [
        ("CYP2D6*39", "CYP2D6*10", va.Relation.IS_CONTAINED),
        ("CYP2D6*39", "CYP2D6*147", va.Relation.IS_CONTAINED),
        ("CYP2D6*10", "CYP2D6*147", va.Relation.IS_CONTAINED),    
    ],
    "color": False
}

image_configs["symmetric"] = {
    "selection": ("CYP2D6*115", "CYP2D6*109"),
    "positions": ((0, 0), (100, 0)),
    "edges": [
        ("CYP2D6*115", "CYP2D6*109", va.Relation.OVERLAP),
        ("CYP2D6*109", "CYP2D6*115", va.Relation.OVERLAP)
    ],
    "color": False
}

image_configs["reflexive"] = {
    "selection": ("CYP2D6*109",),
    "positions": ((0, 0),),
    "edges": [
        ("CYP2D6*109", "CYP2D6*109", va.Relation.EQUIVALENT)
    ],
    "color": False
}

image_configs["most-specific"] = {
    "selection": ("CYP2D6*17", "CYP2D6*82", "CYP2D6*58"),
    "positions": ((0, 0), (0, 100), (100, 100)),
    "edges": [
        ("CYP2D6*17", "CYP2D6*82", va.Relation.IS_CONTAINED),
        ("CYP2D6*82", "CYP2D6*58", va.Relation.OVERLAP),
        ("CYP2D6*17", "CYP2D6*58", va.Relation.OVERLAP)
    ],
    "color": False
}

image_configs["common-ancestor"] = {
    "selection": ("CYP2D6*9", "CYP2D6*109", "CYP2D6*115"),
    "positions": ((0, 0), (100, 100), (0, 100)),
    "edges": [
        ("CYP2D6*9", "CYP2D6*109", va.Relation.IS_CONTAINED),
        ("CYP2D6*9", "CYP2D6*115", va.Relation.IS_CONTAINED),
        ("CYP2D6*109", "CYP2D6*115", va.Relation.OVERLAP)
    ],
    "color": False
}

image_configs["edge-contraction"] = {
    "selection": ("CYP2D6*44", "CYP2D6*22", "CYP2D6*44.001"),
    "positions": ((0, 0), (100, 100), (0, 100)),
    "edges": [
        ("CYP2D6*44", "CYP2D6*44.001", va.Relation.EQUIVALENT),
        ("CYP2D6*44", "CYP2D6*22", va.Relation.IS_CONTAINED),
        ("CYP2D6*44.001", "CYP2D6*22", va.Relation.IS_CONTAINED),
    ],
    "color": False
}

image_configs["node-coreallele"] = {
    "selection": ("CYP2D6*2",),
    "edges": [],
    "color": False
}
image_configs["node-suballele"] = {
    "selection": ("CYP2D6*2.001",),
    "edges": [],
    "color": False
} 
image_configs["node-variant"] = {
    "selection": ("NC_000022.11:g.42126611C>G",),
    "edges": [],
    "color": False
} 
image_configs["node-input"] = {
    "selection": ("HG00373_A",),
    "edges": [],
    "color": False
} 
image_configs["node-personal-variant"] = {
    "selection": ("42125924>G",),
    "edges": [],
    "color": False
}