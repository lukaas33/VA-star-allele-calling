import algebra as va
import json
from dash import Dash, html, dcc, Input, Output
import dash_cytoscape as cyto
from .data import cache_get, cache_set
from .parse import parse_multi_hgvs
from .assets.graph_styles import default_stylesheet, selection_stylesheet

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
    min_edges = []
    tree = Graph(nodes)
    for edge in edges:
        tree.addEdge(edge[0], edge[1])
        if tree.isCyclic():
            tree.removeEdge(edge[0], edge[1])
            continue
        min_edges.append(edge)
    return min_edges

def prune_relations(allele_names, relations):
    """Prune relations which are redundant.

    Useful for displaying as a graph.
    Symmetric relations should not be displayed twice (disjoint, equivalent, overlap).
    Reflexive relations should not be displayed (equivalence to self).
    Transitive relations can be reduced to a tree structure (equivalence, containment).
    Disjoint relation can be left out.

    returns list of nodes and edges.
    """
    transitive = (va.Relation.EQUIVALENT, va.Relation.CONTAINS, va.Relation.IS_CONTAINED)
    symmetric = (va.Relation.EQUIVALENT, va.Relation.OVERLAP, va.Relation.DISJOINT)
    nodes = allele_names[:]
    edges = []
    check_symmetric = set() 
    check_transitivity = {r: [] for r in transitive}
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
            # Only show one (arbitrary) direction of containment 
            if relation == va.Relation.IS_CONTAINED:
                continue
            # Don't display symmetric relations twice
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
    print(len(relations) - len(edges))
    return nodes, edges


def display_graph(nodes, edges):
    """Display relations as a graph

    Uses dash Cytoscape which creates a localhost website.
    The underlying framework is Cytoscape.js, a standard tool in biological network visualization.
    https://dash.plotly.com/cytoscape
    """
    # TODO switch to js for a more extensive app
    # Convert to proper format
    # TODO add more info?
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
    # Setup graph webpage
    cyto.load_extra_layouts()
    app = Dash(__name__)
    default_layout = 'concentric' # TODO find best
    app.layout = html.Div([
        html.Button('Reset view', id='reset-view'),
        html.Button('Reset selection', id='reset-selection'),
        cyto.Cytoscape(
            id='graph',
            style = {
                "width": "100%",
                "height": "75vh"
            },
            layout = {"name": default_layout},
            stylesheet = default_stylesheet,
            elements = elements
        ),
        dcc.Dropdown(
            id='change-layout',
            value=default_layout,
            clearable=False,
            options=[
                # Layouts which load efficiently enough
                {'label': name.capitalize(), 'value': name}
                for name in ['grid', 'random', 'circle', 'cose', 'concentric', 'cola', 'spread', 'breadthfirst']
            ]
        ),
        html.Pre(id='data'),
    ])
    # TODO add filters 
    #    TODO filter out hubs?
    # TODO make expanding?
    # Add interactive component callbacks 
    # Change layout
    @app.callback(Output('graph', 'layout'), Input('change-layout', 'value'))
    def update_layout(layout):
        return {'name': layout}
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
    # Start webpage
    app.run_server(debug=True)