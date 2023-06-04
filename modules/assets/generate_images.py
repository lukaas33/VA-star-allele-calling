# TODO move to scripts?
import algebra as va 

image_configs = {}

image_configs["NA07348_B-simple"] = {
    "selection": ("NA07348_B", "CYP2D6*169", "NC_000022.11:g.42126962C>G", "NC_000022.11:g.42126963C>T"),
    "edges": set([
        ("NA07348_B", "NC_000022.11:g.42126963C>T", va.Relation.EQUIVALENT),
        ("CYP2D6*169", "NC_000022.11:g.42126962C>G", va.Relation.EQUIVALENT),
        ("CYP2D6*169", "NA07348_B", va.Relation.OVERLAP),
    ]),
    "color": True,
    "layout": "dagre"
}

image_configs["NA07348_B-extended"] = {
    "selection": ("NA07348_B", "CYP2D6*169", "CYP2D6*169.001", "NC_000022.11:g.42126962C>G", "NC_000022.11:g.42126963C>T", "CYP2D6*1.002"),
    "edges": set([
        ("NA07348_B", "CYP2D6*1.002", va.Relation.EQUIVALENT),
        ("CYP2D6*1.002", "NC_000022.11:g.42126963C>T", va.Relation.EQUIVALENT),
        ("CYP2D6*1.002", "CYP2D6*169", va.Relation.OVERLAP),
        ("CYP2D6*169", "NC_000022.11:g.42126962C>G", va.Relation.EQUIVALENT),
        ("CYP2D6*169", "CYP2D6*169.001", va.Relation.EQUIVALENT),
    ]),
    "color": True,
    "layout": "dagre"
}

image_configs["transitive"] = {
    "selection": ("CYP2D6*39", "CYP2D6*10", "CYP2D6*147"),
    "positions": ((0, 0), (0, 100), (100, 100)),
    "edges": set([
        ("CYP2D6*39", "CYP2D6*10", va.Relation.IS_CONTAINED),
        ("CYP2D6*39", "CYP2D6*147", va.Relation.IS_CONTAINED),
        ("CYP2D6*10", "CYP2D6*147", va.Relation.IS_CONTAINED),    
    ]),
    "color": False
}

image_configs["symmetric"] = {
    "selection": ("CYP2D6*115", "CYP2D6*109"),
    "positions": ((0, 0), (100, 0)),
    "edges": set([
        ("CYP2D6*115", "CYP2D6*109", va.Relation.OVERLAP),
        ("CYP2D6*109", "CYP2D6*115", va.Relation.OVERLAP)
    ]),
    "color": False
}

image_configs["reflexive"] = {
    "selection": ("CYP2D6*109",),
    "positions": ((0, 0),),
    "edges": set([
        ("CYP2D6*109", "CYP2D6*109", va.Relation.EQUIVALENT)
    ]),
    "color": False
}

image_configs["most-specific"] = {
    "selection": ("CYP2D6*17", "CYP2D6*82", "CYP2D6*58"),
    "positions": ((0, 0), (0, 100), (100, 100)),
    "edges": set([
        ("CYP2D6*17", "CYP2D6*82", va.Relation.IS_CONTAINED),
        ("CYP2D6*82", "CYP2D6*58", va.Relation.OVERLAP),
        ("CYP2D6*17", "CYP2D6*58", va.Relation.OVERLAP)
    ]),
    "color": False
}

image_configs["common-ancestor"] = {
    "selection": ("CYP2D6*9", "CYP2D6*109", "CYP2D6*115"),
    "positions": ((0, 0), (100, 100), (0, 100)),
    "edges": set([
        ("CYP2D6*9", "CYP2D6*109", va.Relation.IS_CONTAINED),
        ("CYP2D6*9", "CYP2D6*115", va.Relation.IS_CONTAINED),
        ("CYP2D6*109", "CYP2D6*115", va.Relation.OVERLAP)
    ]),
    "color": False
}

image_configs["edge-contraction"] = {
    "selection": ("CYP2D6*44", "CYP2D6*22", "CYP2D6*44.001"),
    "positions": ((0, 0), (100, 100), (0, 100)),
    "edges": set([
        ("CYP2D6*44", "CYP2D6*44.001", va.Relation.EQUIVALENT),
        ("CYP2D6*44", "CYP2D6*22", va.Relation.IS_CONTAINED),
        ("CYP2D6*44.001", "CYP2D6*22", va.Relation.IS_CONTAINED),
    ]),
    "color": False
}

image_configs["node-coreallele"] = {
    "selection": ("CYP2D6*2",),
    "edges": set([]),
    "color": False
}
image_configs["node-suballele"] = {
    "selection": ("CYP2D6*2.001",),
    "edges": set([]),
    "color": False
} 
image_configs["node-variant"] = {
    "selection": ("NC_000022.11:g.42126611C>G",),
    "edges": set([]),
    "color": False
} 
image_configs["node-input"] = {
    "selection": ("HG00373_A",),
    "edges": set([]),
    "color": False
} 
image_configs["node-personal-variant"] = {
    "selection": ("42125924>G",),
    "edges": set([]),
    "color": False
}

image_configs["function-normal"] = {
    "selection": ("CYP2D6*2",),
    "edges": set([]),
    "color": True
}
image_configs["function-no"] = {
    "selection": ("CYP2D6*4",),
    "edges": set([]),
    "color": True
}
image_configs["function-decreased"] = {
    "selection": ("CYP2D6*10",),
    "edges": set([]),
    "color": True
}
image_configs["function-uncertain"] = {
    "selection": ("CYP2D6*22",),
    "edges": set([]),
    "color": True
}
image_configs["function-unknown"] = {
    "selection": ("CYP2D6*58",),
    "edges": set([]),
    "color": True
}
image_configs["function-uncertain"] = {
    "selection": ("CYP2D6*22",),
    "edges": set([]),
    "color": True
}
image_configs["function-na"] = {
    "selection": ("CYP2D6*162",),
    "edges": set([]),
    "color": True
}
image_configs["impact-3"] = {
    "selection": ("NC_000022.11:g.42128945C>T",),
    "edges": set([]),
    "color": True
}

image_configs["impact-2"] = {
    "selection": ("NC_000022.11:g.42130692G>A",),
    "edges": set([]),
    "color": True
}

image_configs["impact-1"] = {
    "selection": ("NC_000022.11:g.42131469C>T",),
    "edges": set([]),
    "color": True
}

image_configs["impact-0"] = {
    "selection": ("42132027delinsCT",),
    "edges": set([]),
    "color": True
}