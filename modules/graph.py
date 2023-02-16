import algebra as va
from dash import Dash, html
import dash_cytoscape as cyto
from .data import cache_get, cache_set
from .parse import parse_multi_hgvs

def find_relations(corealleles, reference_sequence):
    """Find the relation between all corealleles.

    Relations are cached since they take a long time to generate.

    Returns nodes and edges.
    """
    # Find relation for each pair of corealleles directionally
    coreallele_names = list(corealleles.keys())
    name = "relations"
    try:
        relations = cache_get(name)
    except:
        relations = {
            coreallele: {
                coreallele2: None 
                for coreallele2 in coreallele_names
            } 
            for coreallele in coreallele_names
        }
        for i, left_coreallele in enumerate(coreallele_names):
            # TODO change loop to remove all none values
            for right_coreallele in coreallele_names[i+1:]: # Only check each pair once
                # Parse allele
                left_hgvs = [variant["hgvs"] for variant in corealleles[left_coreallele]["variants"]]
                left_variants = parse_multi_hgvs(left_hgvs, reference_sequence, left_coreallele)
                right_hgvs = [variant["hgvs"] for variant in corealleles[right_coreallele]["variants"]]
                right_variants = parse_multi_hgvs(right_hgvs, reference_sequence, right_coreallele)
                # Store relation
                relation = va.compare(reference_sequence, left_variants, right_variants)
                relations[left_coreallele][right_coreallele] = relation
                # Store inverse relation (the same for most relations)
                if relation == va.Relation.CONTAINS:
                    inv_relation = va.Relation.IS_CONTAINED
                elif relation == va.Relation.IS_CONTAINED:
                    inv_relation = va.Relation.CONTAINS
                else:
                    inv_relation = relation
                relations[right_coreallele][left_coreallele] = inv_relation
        cache_set(relations, name) # TODO cache later
    pruned = prune_relations(coreallele_names, relations)
    return pruned

def prune_relations(allele_names, relations):
    """Prune relations which are redundant.

    Useful for displaying as a graph.
    Symmetric relations should not be displayed twice (disjoint, equivalent, overlap).
    Reflexive relations should not be displayed (equivalence to self).
    Transitive relations can be reduced to a tree structure (equivalence, containment).
    Disjoint relation can be left out.

    returns list of nodes and edges.
    """
    nodes = allele_names[:]
    edges = []
    check_symmetric = set() 
    check_transitivity = {r: [] for r in va.Relation}
    for node in allele_names:
        for other in relations[node].keys():
            relation = relations[node][other]
            # Skip trivial self equivalence (reflexivity)
            if node == other: 
                continue
            if relation is None:
                raise ValueError("Relation data is incomplete")
            # Don't display disjointedness explicitly
            if relation == va.Relation.DISJOINT: 
                continue
            # Don't display symmetric relations twice
            if relation in (va.Relation.EQUIVALENT, va.Relation.OVERLAP):
                pair = (node, other)
                inv_pair = (other, node)
                if pair in check_symmetric or inv_pair in check_symmetric: 
                    continue
                check_symmetric.add(pair)
                check_symmetric.add(inv_pair)
            # Store transitive relations to prune and add later
            if relation in (va.Relation.EQUIVALENT, va.Relation.CONTAINS, va.Relation.IS_CONTAINED):
                check_transitivity[relation].append((node, other))
                continue
            edges.append((node, other, relation))

    print(len(nodes), len(edges))
    return nodes, edges


def display_graph(nodes, edges):
    """Display relations as a graph

    Uses dash Cytoscape which creates a localhost website.
    The underlying framework is Cytoscape.js, a standard tool in biological network visualization.
    https://dash.plotly.com/cytoscape
    """
    # Convert to proper format
    elements = []
    for node in nodes:
        elements.append({            
            "data": {
                "id": node, 
                "label": node.split("CYP2D6")[1]
            }
        })
    for node, other, relation in edges:
        elements.append({
            "data": {
                "source": node,
                "target": other
            },
            "classes": str(relation).split('.')[1]
        })
    # Start webpage
    cyto.load_extra_layouts()
    app = Dash(__name__)
    app.layout = html.Div([
        cyto.Cytoscape(
            # TODO highlight connected edges for selected node
            # TODO add filters 
            # TODO filter out hubs?
            id='graph',
            # TODO find layout concentric, breadthfirst
            layout={'name': 'concentric'}, 
            # TODO make full screen
            style = {
                "width": "100%",
                "height": "400px"
            },
            # TODO make stylesheet external
            # TODO use colors/symbols and add legend 
            stylesheet = [
                {
                    'selector': 'node',
                    'style': {
                        'content': 'data(label)'
                    }
                }, {
                    'selector': 'edge',
                    'style': {
                        'curve-style': 'bezier'
                    }
                }, {
                    'selector': '.EQUIVALENT',
                    'style': {
                        'line-color': 'black'
                    }
                }, {
                    'selector': '.CONTAINS',
                    'style': {
                        'target-arrow-shape': 'triangle',
                        'target-arrow-color': 'grey',
                        'line-color': 'grey'
                    }
                }, {
                    'selector': '.IS_CONTAINED',
                    'style': {
                        'target-arrow-shape': 'triangle',
                        'target-arrow-color': 'grey',
                        'line-color': 'grey'
                    }
                }, {
                    'selector': '.OVERLAP',
                    'style': {
                        'line-color': 'lightgrey'
                    }
                }
            ],
            elements=elements
        )
    ])
    app.run_server(debug=True)