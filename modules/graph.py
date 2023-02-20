import algebra as va
import json
from dash import Dash, html, dcc, Input, Output, ctx
import dash_cytoscape as cyto
import plotly.express as px
from pandas import DataFrame as df
from .data import cache_get, cache_set
from .parse import parse_multi_hgvs
from .analyse import count_arity, count_relations
from .assets.graph_styles import default_stylesheet, selection_stylesheet

# TODO use OOP graph class for everything?
# TODO split into more files

def find_relations(corealleles, reference_sequence):
    """Find the relation between all corealleles.

    Relations are cached since they take a long time to generate.

    Returns list of edges.
    """
    # Find relation for each pair of corealleles directionally
    coreallele_names = list(corealleles.keys())
    name = "relations"
    try:
        edges = cache_get(name)
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
            for right_coreallele in coreallele_names[i:]: # Only check each pair once
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
        edges = []
        for l_node in corealleles:
            for r_node in corealleles:
                edges.append((l_node, r_node, relations[l_node][r_node]))
        cache_set(edges, name) # Cache
    return edges

class Graph():
    def __init__(self, nodes):
        self.nodes = nodes
        self.graph = {node: [] for node in nodes}
 
    def addEdge(self, node, edge):
        self.graph[node].append(edge)

    def removeEdge(self, node, edge):
        self.graph[node].remove(edge)
 
    def isCyclicUtil(self, v, visited, stack):
        # Mark current node as visited and
        # adds to recursion stack
        visited[v] = True
        stack[v] = True
 
        # Recur for all neighbours
        # if any neighbour is visited and in
        # stack then graph is cyclic
        for neighbour in self.graph[v]:
            if visited[neighbour] == False:
                if self.isCyclicUtil(neighbour, visited, stack) == True:
                    return True
            elif stack[neighbour] == True:
                return True
 
        # The node needs to be popped from
        # recursion stack before function ends
        stack[v] = False
        return False
 
    # Returns true if graph is cyclic else false
    def isCyclic(self):
        visited = {node: False for node in self.nodes}
        stack = {node: False for node in self.nodes}
        for node in self.nodes:
            if visited[node] == False:
                if self.isCyclicUtil(node, visited, stack) == True:
                    return True
        return False

def spanning_tree(nodes, edges):
    """Find spanning tree given a graph.
    
    Can be used to minimize graph structure with transitive relations.
    """
    # TODO delete, NOT a correct approach
    # TODO use existing algorithms
    min_edges = []
    tree = Graph(nodes)
    for edge in edges:
        tree.addEdge(edge[0], edge[1])
        if tree.isCyclic():
            tree.removeEdge(edge[0], edge[1])
            continue
        min_edges.append(edge)
    return min_edges

def prune_relations(nodes, relations):
    """Prune relations which are redundant.

    Symmetric relations should not be displayed twice (disjoint, equivalent, overlap).
    Reflexive relations should not be displayed (equivalence to self).
    Transitive relations can be reduced to a tree structure (equivalence, containment).
    Disjoint relation can be left out.

    returns list of edges.
    """
    transitive = (va.Relation.EQUIVALENT, va.Relation.CONTAINS, va.Relation.IS_CONTAINED)
    symmetric = (va.Relation.EQUIVALENT, va.Relation.OVERLAP, va.Relation.DISJOINT)
    edges = []
    check_symmetric = set() 
    check_transitivity = {r: [] for r in transitive}
    for node, other, relation in relations:
        if relation is None:
            raise ValueError("Relation data is incomplete")
        # Skip trivial self equivalence (reflexivity)
        if node == other: 
            continue
        # Don't represent disjointedness explicitly, 
        # it can be seen as 'no relation'
        if relation == va.Relation.DISJOINT: 
            continue
        # Don't represent symmetric relations twice
        if relation in symmetric:
            pair = (node, other)
            inv_pair = (other, node)
            if pair in check_symmetric or inv_pair in check_symmetric: 
                continue
            check_symmetric.add(pair)
            check_symmetric.add(inv_pair)
        # Store transitive relations to prune and add later
        if relation in transitive:
            check_transitivity[relation].append((node, other, relation))
            continue
        edges.append((node, other, relation))
    # Reduce transitive graph to tree and add
    for relation, subset_edges in check_transitivity.items():
        edges += spanning_tree(nodes, subset_edges)
    return edges

def plot_arity(nodes, relations):
    """Create a plot of the arity values"""
    arities = count_arity(nodes, relations)
    # Covert data to pandas dataframe
    data = {"allele": [], "relation": [], "arity": []}
    for rel in va.Relation:
        for node in nodes:
            data["allele"].append(node)
            data["relation"].append(str(rel).split('.')[1])
            data["arity"].append(arities[node][rel])
    data = df(data)
    # Display
    figure = px.bar(df(data), x="allele", y="arity", color="relation", title="Pruned relationships")
    figure.update_layout(barmode='stack', xaxis={'categoryorder': 'total ascending'})
    return figure


def display_graph(nodes, relations, data):
    """Display relations as a graph

    Uses dash Cytoscape which creates a localhost website.
    The underlying framework is Cytoscape.js, a standard tool in biological network visualization.
    https://dash.plotly.com/cytoscape
    """
    # print(f"All {len(relations)} relationships:")
    count_all = count_relations(relations)
    # for relationType, count in count_all.items():
    #     print(f"  {count} {relationType}")
    # TODO switch to js for a more extensive app (filtering, searching, expanding, etc.)
    # Prune redundant relations
    edges = prune_relations(nodes, relations)
    # print(f"Pruned {len(edges)} relationships for graph:")
    count_pruned = count_relations(edges)
    # for relationType, count in count_pruned.items():
    #     print(f"  {count} {relationType}")
    # Convert to proper format
    elements = []
    for node in nodes:
        elements.append({            
            "data": {
                "id": node, 
                "label": node.split("CYP2D6")[1],
                "data": data[node]
            }
        })
    for node, other, relation in edges:
        elements.append({
            "data": {
                "source": node,
                "target": other,
            },
            "classes": str(relation).split('.')[1]
        })
    # Setup graph webpage
    cyto.load_extra_layouts()
    app = Dash(__name__)
    default_layout = 'cose-bilkent' 
    app.layout = html.Div([
        dcc.Location(id='url', refresh=False), # Page load
        dcc.Tabs([
            dcc.Tab(
                label="Network",
                children=[
                    cyto.Cytoscape(
                        id='graph',
                        style = {
                            "width": "100%",
                            "height": "75vh"
                        },
                        stylesheet = default_stylesheet,
                        elements = elements
                    ),
                    html.Button('Reset view', id='reset-view'),
                    html.Button('Reset selection', id='reset-selection'),
                    html.Button("Export image", id="image-svg"),
                    dcc.Dropdown(
                        id='change-layout',
                        value=default_layout,
                        clearable=False,
                        options=[
                            # Layouts which load efficiently enough
                            {'label': name.capitalize(), 'value': name}
                            for name in ['grid', 'random', 'circle', 'cose', 'concentric', 'cola', 'spread', 'breadthfirst', 'cose-bilkent']
                        ]
                    ),
                    html.Pre(id='data'),
                ]
            ),
            dcc.Tab(
                label="Statistics",
                children=[
                    dcc.Graph(id="plot-arity")
                ]
            )
        ])        
    ])
    # Add interactive component callbacks 
    # Change layout
    @app.callback(Output('graph', 'layout'), Input('change-layout', 'value'))
    def update_layout(layout):
        settings = {"name": layout}
        settings["nodeDimensionsIncludeLabels"] = True
        if layout == 'cose-bilkent' or layout == 'cose':
            settings["idealEdgeLength"] = 400
            settings["tile"] = False
            settings["animate"] = False
        return settings
    # Display information about selection
    @app.callback(Output('data', 'children'), Input('graph', 'tapNodeData'))
    def displayTapNodeData(data):
        return json.dumps(data, indent=2)
    # Display connections of selected
    @app.callback([Output('graph', 'stylesheet'), Output('reset-selection', 'n_clicks')], [Input('graph', 'tapNode'), Input('reset-selection', 'n_clicks')])
    def generate_stylesheet(node, n_clicks):
        if not node or n_clicks == 1: # No input or resetting
            return [default_stylesheet, None]
        return [selection_stylesheet(node), None]
    # Reset view
    @app.callback([Output('graph', 'zoom'), Output('graph', 'elements')], [Input('reset-view', 'n_clicks')])
    def reset_layout(n_clicks):
        return [1, elements]
    # Show plots
    @app.callback(Output('plot-arity', "figure"), Input('url', 'pathname'))
    def show_arity_plot(_):
        return plot_arity(nodes, edges)
    # Export image
    @app.callback(
        Output("graph", "generateImage"),
        [
            Input("image-svg", "n_clicks"),
        ])
    def get_image(n_clicks):
        if n_clicks is None:
            return { # TODO why needed?
                'type': 'png',
                'action': 'store'
            }
        return {
            'type': 'svg',
            'action': 'download'
            }
    # Start webpage
    app.run_server(debug=True)