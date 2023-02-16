import algebra as va
from dash import Dash, html
import dash_cytoscape as cyto
from modules.data import cache_get, cache_set
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
    Reflexive relations (equivalence to self) should not be displayed.
    Transitive relations can be reduced by using a tree structure (equivalence, containment).
    Disjoint relation can be left out.

    returns list of nodes and edges.
    """
    nodes = []
    for allele_name in allele_names:
        nodes.append({"data": {
            "id": allele_name, 
            "label": allele_name
        }})
    edges = []
    check_symmetric = set() # Used to check if a reflexive relation was present already
    for node in allele_names:
        for other in relations[node].keys():
            relation = relations[node][other]
            # Skip trivial self equivalence (reflexivity)
            if node == other: 
                continue
            if relation is None:
                raise ValueError("Relation data is incomplete")
            # Don't display disjointness explicitly
            if relation == va.Relation.DISJOINT: 
                continue
            # Don't display symmetric relations twice
            if relation == va.Relation.EQUIVALENT or relation == va.Relation.OVERLAP:
                pair = (node, other)
                inv_pair = (other, node)
                if pair in check_symmetric or inv_pair in check_symmetric: 
                    continue
                check_symmetric.add(pair)
                check_symmetric.add(inv_pair)
            # TODO prune transitive
            edges.append({"data": {
                "source": node,
                "target": other,
                "label": str(relation),
                "classes": str(relation)
            }})
    print(len(nodes), len(edges))
    return nodes, edges


def display_graph(nodes, edges):
    """Display relations as a graph

    Uses dash Cytoscape which creates a localhost website.
    The underlying framework is Cytoscape.js, a standard tool in biological network visualization.
    https://dash.plotly.com/cytoscape
    """
    app = Dash(__name__)
    app.layout = html.Div([
        cyto.Cytoscape(
            # TODO highlight connected edges for selected node
            id='graph',
            layout={'name': 'cose'},
            # TODO make fullscreen
            style={},
            # TODO make stylesheet external
            # TODO use colours and add legend 
            stylesheet = [
                {
                    'selector': 'node',
                    'style': {
                        'content': 'data(label)'
                    }
                },
                {
                    'selector': 'edge',
                    'style': {
                        # 'content': 'data(label)',
                        'curve-style': 'bezier',
                        'target-arrow-shape': 'none',
                        'source-arrow-shape': 'none',
                    }
                },
                {
                    'selector': '.Relation.EQUIVALENT',
                    'style': {
                    }
                },
                {
                    'selector': '.Relation.CONTAINS',
                    'style': {
                        'target-arrow-shape': 'triangle'

                    }
                },
                {
                    'selector': '.Relation.IS_CONTAINED',
                    'style': {
                        'source-arrow-shape': 'triangle'
                    }
                },
                {
                    'selector': '.Relation.OVERLAP',
                    'style': {
                    }
                }
            ],
            elements=(nodes + edges)
        )
    ])
    app.run_server(debug=True)