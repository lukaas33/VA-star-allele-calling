import algebra as va
import json
from dash import Dash, html, dcc, Input, Output, no_update
import dash_cytoscape as cyto
import plotly.express as px
import pandas
from .utils import count_arity, count_relations
from .assets.graph_styles import default_stylesheet, selection_stylesheet
from .relations import prune_relations, find_context
from .other_sources import severity_GO, severity_pharmvar
from .calling import sort_types

def plot_arity(nodes, relations):
    """Create a plot of the arity values"""
    arity = count_arity(nodes, relations)
    # Covert data to pandas data frame
    data = {"allele": [], "relation": [], "arity": []}
    for node in nodes:
        if "*" not in node:
            continue
        for rel in va.Relation:
            data["allele"].append(node)
            data["relation"].append(rel.name)
            data["arity"].append(arity[node][rel.name])
    data = pandas.DataFrame(data)
    # Display
    figure = px.bar(pandas.DataFrame(data), x="allele", y="arity", color="relation", title="Pruned relationships between all variants and alleles, per allele")
    figure.update_layout(barmode='stack', xaxis={'categoryorder': 'total ascending'})
    return figure

def plot_relations(relations, pruned=False):
    data = {"relation": [], "count": []}
    counts = count_relations(relations)
    for relation in counts.keys():
        data["relation"].append(relation)
        data['count'].append(counts[relation])
    title = "Relationships between all variants and alleles"
    if pruned: title = "Pruned " + title.lower()
    figure = px.bar(pandas.DataFrame(data), x="relation", y="count", title=title, log_y=(not pruned))
    return figure

def plot_counts(elements):
    data = {"category": ["core", "sub", "variant"], "count": [0, 0, 0]}
    for element in elements:
        if element["classes"] in data["category"]:
            data["count"][data["category"].index(element["classes"])] += 1
    figure = px.bar(pandas.DataFrame(data), x="category", y="count", title="Amount of alleles and variants")
    return figure

def layout_graph(elements, nodes, edges, default_layout='cose-bilkent'):
    """Returns the layout for the Dash graph"""
    layouts = list(set(['grid', 'random', 'circle', 'concentric', 'cola', 'spread', 'breadthfirst', 'cose-bilkent', "dagre", "euler", "klay", default_layout]))
    layouts.sort()
    return dcc.Tabs([
        dcc.Tab(
            label="Network",
            children=[
                html.Div(
                    id='settings',
                    style = {
                        "display": "flex", 
                        "flexWrap": "wrap",
                        "flex-direction": "column",
                        "width": "10vw"
                    },
                    children = [
                        dcc.Input(
                            id="filter",
                            placeholder="CYP2D6*4, ...",
                            type="text",
                            debounce=True
                        ),
                        dcc.Dropdown(
                            id='change-layout',
                            value=default_layout,
                            clearable=False,
                            options=[{'label': name.capitalize(), 'value': name} for name in layouts]
                        ),
                        dcc.Dropdown(
                            id='image-type',
                            value='svg',
                            clearable=False,
                            options=[
                                {'label': 'SVG', 'value': 'svg'},
                                {'label': 'PNG', 'value': 'png'},
                                {'label': 'JPEG', 'value': 'jpeg'},
                            ]
                        ),
                        html.Button("Export image", id="image"),
                        html.Button('Subgraph selection', id='subgraph'),
                    ]
                ),
                cyto.Cytoscape(
                    id='graph',
                    style = {
                        "height": "75vh",
                        "width": "90vw"
                    },
                    layout={
                        "name": default_layout,
                        "fit": default_layout != 'preset' 
                    },
                    stylesheet = default_stylesheet,
                    elements = elements
                ),
            ]
        ),
        dcc.Tab(
            label="Statistics",
            children=[
                dcc.Graph(
                    id="plot-arity",
                    figure=plot_arity(nodes, edges)
                ),
                # TODO fix?
                # dcc.Graph(
                #     id="relation-counts",
                #     figure=plot_relations(edges)
                # ),
                dcc.Graph(
                    id="pruned-relation-counts",
                    figure=plot_relations(edges, pruned=True)
                ),                    
                dcc.Graph(
                    id="counts",
                    figure=plot_counts(elements)
                )
            ]
        ),
        dcc.Tab(
            label="Selection data",
            children=[
                html.Pre(id='data')
            ]
        )
    ])

def interactive_graph(app, original_elements, edges):
    """Add interactive components to graph"""
    # Change layout
    @app.callback(
        Output('graph', 'layout'), 
        [Input('graph', 'layout'), Input('change-layout', 'value')])
    def update_layout(current_layout, new_layout):
        if current_layout is None:
            return no_update
        if current_layout["name"] == new_layout:
            return no_update
        settings = {}
        settings["name"] = new_layout
        if new_layout == 'cose-bilkent' or new_layout == 'cose':
            settings["nodeDimensionsIncludeLabels"] = True
            settings["tile"] = False
        settings["animate"] = False
        return settings
    # Display information about selection
    @app.callback(
        Output('data', 'children'), 
        Input('graph', 'selectedNodeData'))
    def displayTapNodeData(data):
        if not data:
            return no_update
        return json.dumps(data, indent=2)
    # Display connections of selected
    @app.callback(
        Output('graph', 'stylesheet'), 
        [Input('graph', 'selectedNodeData'), Input('graph', 'layout')])
    def generate_stylesheet(nodes, layout):
        if not nodes: # No input or resetting
            return default_stylesheet
        # SHow both incoming and outgoing containment (different from subgraph)
        # TODO use same neighbourhood definition?
        context = find_context([node["id"] for node in nodes], edges)
        return selection_stylesheet(context, layout["name"])
    @app.callback(
        Output("graph", "generateImage"),
        [Input("image", "n_clicks"), Input("image-type", "value")])
    def get_image(n_clicks, type):
        # TODO fix downloading on change
        if n_clicks is None:
            return no_update
        return {
            'type': type,
            'action': 'download'
            }
    # Subgraph selection
    @app.callback(
            [Output('graph', 'elements'), Output('subgraph', 'n_clicks')], 
            [Input('graph', 'selectedNodeData'), Input('subgraph', 'n_clicks')])
    def subgraph(nodes, n_clicks):
        if n_clicks is None: # On initial load
            return [no_update, None]
        if n_clicks >= 1 and nodes == []: # Reset
            # TODO fix resetting from filter
            return [original_elements, None]
        # Only shows incoming containments in subgraph
        context = find_context([node["id"] for node in nodes], edges, directional=True)
        # Translate to elements
        selected_elements = []
        for element in original_elements:
            if "id" in element["data"].keys() and element["data"]["id"] in context:
                selected_elements.append(element)
            elif "source" in element["data"].keys() and element["data"]["source"] in context and \
                    "target" in element["data"].keys() and element["data"]["target"] in context:
                selected_elements.append(element)
        return [selected_elements, 1]
    # Filter
    @app.callback(
        Output('graph', 'selectedNodeData'),
        Input('filter', 'value'))
    def filter(filterValues):
        # TODO inexact selection
        # TODO fix filter
        if not filterValues:
            return no_update
        filterValues = [filterValue.strip() for filterValue in filterValues.split(",")]
        selection = []
        for element in original_elements:
            if "id" not in element["data"].keys():
                continue
            if any([filterValue == element["data"]["id"] for filterValue in filterValues]):
                selection.append(element["data"])
        return selection

def display_graph(nodes, edges, data, functions, positions=None, default_layout="cose-bilkent"):
    """Display relations as a graph

    Uses dash Cytoscape which creates a localhost website.
    The underlying framework is Cytoscape.js, a standard tool in biological network visualization.
    https://dash.plotly.com/cytoscape
    This framework should not be used to make a complete visualization since it is limited in functionality.
    """
    # TODO create js page
    # TODO why does it call this function twice
    # TODO add link between phased samples
    # TODO add display filters
    print("Displaying graph")
    # Convert to proper format for cytoscape
    elements = []
    for i, node in enumerate(nodes):
        # TODO use enum or classes
        function, impact, severity = None, None, None
        if sort_types(node) == 1: # Core allele
            category = "core"
            label = "*" + node.split("*")[1]
            function = functions[node]
        elif sort_types(node) == 2: # Suballele
            category = "sub"
            label = "*" + node.split("*")[1]
            function = functions[node]
        elif sort_types(node) == 3: # PharmVar variant
            category = "variant"
            label = node.split(':')[1].split('.')[1]
            impact = functions[node]
            severity = severity_pharmvar(functions[node])
        elif sort_types(node) == 4: # Sample
            category = "sample"
            label = node
        elif sort_types(node) == 5: # Personal variant
            category = "observed variant"
            label = node
            impact = "; ".join(functions[node])
            severity = severity_GO(functions[node])
        element = {            
            "data": {
                "id": node, 
                "label": label,
                "data": data[node] if node in data else None,
                "function": function,
                "impact": impact,
                "severity": severity
            },
            "classes": category,
        }
        if positions is not None: element["position"] = positions[i]
        elements.append(element)
    for node, other, relation in edges:
        element = {
            "data": {
                "source": node,
                "target": other,
            },
            "classes": relation.name
        }
        elements.append(element)
    # Setup graph webpage
    cyto.load_extra_layouts()
    app = Dash(__name__)
    # Show graph
    app.layout = layout_graph(elements, nodes, edges, default_layout=default_layout) 
    # Add interactive component callbacks 
    interactive_graph(app, elements, edges)
    # Start webpage
    app.run_server(debug=True)