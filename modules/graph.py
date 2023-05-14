import algebra as va
import json
from dash import Dash, html, dcc, Input, Output, no_update
import dash_cytoscape as cyto
import plotly.express as px
import pandas
import requests
import os
import signal
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
    return html.Div(
        style = {
            "display": "flex",
        },
        children = [
            html.Div(
                style = {
                    "display": "flex", 
                    "flex-direction": "column",
                    "width": "20vw",
                    "height": "100vh",
                    "margin": 0,
                    "box-sizing": "border-box",
                },
                children = [
                    html.P("Search"),
                    dcc.Input(
                        id="filter",
                        placeholder="CYP2D6*4, NC_000022.11:g.42128945C>T, ...",
                        type="text",
                        debounce=True
                    ),
                    html.P("Layout"),
                    dcc.Dropdown(
                        id='change-layout',
                        value=default_layout,
                        clearable=False,
                        options=[{'label': name.capitalize(), 'value': name} for name in layouts]
                    ),
                    html.P("Download"),
                    dcc.Dropdown(
                        id='image-type',
                        value='svg',
                        clearable=False,
                        options=[{'label': name.upper(), 'value': name} for name in ['svg', 'png', 'jpeg']]
                    ),
                    html.Button(
                        "Export image", 
                        id="image"
                    ),
                    html.P("Filter nodes"),
                    html.Button(
                        'Create Subgraph', 
                        id='subgraph',
                    ),
                    html.Button(
                        id='shutdown',
                        style = {
                            "display": "none"
                        }
                    ),
                    html.Button( # TODO remove
                        id='dummy',
                        style = {
                            "display": "none"
                        }
                    ), 
                    html.P("Data"),
                    html.Pre(
                        id='data',
                        style = {
                            "width": "100%",
                            "height": "100%",
                            "overflow-y": "scroll",
                            "overflow-x": "hidden"
                        }
                    )
                ]
            ),
            cyto.Cytoscape(
                id='graph',
                style = {
                    "height": "100vh",
                    "width": "80vw"
                },
                layout={
                    "name": default_layout,
                    "fit": default_layout != 'preset',
                    "nodeDimensionsIncludeLabels": True,
                    "tile": False,
                    "animate": False,
                },
                stylesheet = default_stylesheet,
                elements = elements,
                zoom = 1,
                pan = {"x": 0, "y": 0},
                minZoom = 0.1,
                maxZoom = 10
            ),
        ]
    )

def interactive_graph(app, original_elements, edges, auto_download):
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
        current_layout["name"] = new_layout
        return current_layout
    # Display information about selection
    @app.callback(
        Output('data', 'children'), 
        Input('graph', 'selectedNodeData'))
    def displayTapNodeData(data):
        if not data or "data" not in data:
            return no_update
        return json.dumps(data["data"], indent=2)
    # Display connections of selected
    @app.callback(
        Output('graph', 'stylesheet'), 
        [Input('graph', 'selectedNodeData'), Input('graph', 'layout')])
    def generate_stylesheet(nodes, layout):
        if not nodes: # No input or resetting
            return default_stylesheet
        # show both incoming and outgoing containment (different from subgraph neighbourhood)
        # TODO use same neighbourhood definition?
        context, _ = find_context(set([node["id"] for node in nodes]), edges)
        return selection_stylesheet(context, layout["name"])
    # Download image
    @app.callback(
        [Output("graph", "generateImage"), Output("image", "n_clicks"), Output("shutdown", "n_clicks")],
        [Input("image", "n_clicks"), Input("image-type", "value")])
    def get_image(n_clicks, type):
        if n_clicks is None and not auto_download:
            return no_update
        name = auto_download if auto_download else "image" # TODO use selection
        shutdown = 1 if auto_download else None # Shutdown after download
        return [{'type': type, 'action': 'download', 'filename': name}, None, shutdown]
    # Shutdown server
    @app.callback(
        Output("dummy", "n_clicks"), # No output
        Input("shutdown", "n_clicks"))
    def shutdown(n_clicks):
        if n_clicks is None:
            return no_update
        print('Server shutting down...')
        os.kill(os.getpid(), signal.SIGTERM)
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
        context, _ = find_context(set([node["id"] for node in nodes]), edges, directional=True)
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

def display_graph(nodes, edges, data, functions, positions=None, default_layout="cose-bilkent", relevance=None, auto_download=None, marked_calling=None, group_variants=None):
    """Display relations as a graph

    Uses dash Cytoscape which creates a localhost website.
    The underlying framework is Cytoscape.js, a standard tool in biological network visualization.
    https://dash.plotly.com/cytoscape
    This framework should not be used to make a complete visualization since it is limited in functionality.
    """
    # TODO create js page
    # TODO why does it call this function twice
    # TODO allow for multiple samples
    # TODO add display filters
    # TODO make auto download faster?
    # TODO make calling graph separate function
    print("Displaying graph")
    # Convert to proper format for cytoscape
    elements = []
    # Add nodes
    for i, node in enumerate(nodes):
        function, impact, severity, parent = None, None, None, None
        relevant = True
        if sort_types(node) in (1, 2):
            function = functions[node]
            label = "*" + node.split("*")[1]
            if sort_types(node) == 1: # Core allele
                category = "core"
            elif sort_types(node) == 2: # Suballele
                category = "sub"
            if marked_calling and node in marked_calling:
                category += " called"
        elif sort_types(node) in (3, 5): # Variant
            category = "variant"
            # Show relevance to calling
            relevant = relevance[node] if relevance is not None and node in relevance else True
            # Group variants
            parent = 'variants' if group_variants is not None and node in group_variants else None 
            if sort_types(node) == 3: # PharmVar variant
                label = node.split(':')[1].split('.')[1]
                impact = functions[node]
                severity = severity_pharmvar(functions[node])
            elif sort_types(node) == 5: # Personal variant
                category += " observed"
                label = node
                impact = "; ".join(functions[node])
                severity = severity_GO(functions[node])
        elif sort_types(node) == 4: # Sample
            category = "sample"
            label = node
        element = {            
            # TODO don't store all fields of data
            "data": {
                "id": node, 
                "label": label,
                "function": function,
                "impact": impact,
                "severity": severity,
                "relevant": relevant,
                "parent": parent,
                "data": data[node] if node in data else None,
            },
            "classes": category,
        }
        if positions is not None: element["position"] = positions[i]
        elements.append(element)
    # Treat group as single node
    if group_variants is not None and len(group_variants) > 0: 
        elements.append({"data": {"id": "variants"}})
        # Attach edges to some in group to group
        # TODO only attach edges to all?
        for source, target, relation in list(edges):
            if source in group_variants: 
                edges.remove((source, target, relation))
                edges.add(("variants", target, relation))
            if target in group_variants:
                edges.remove((source, target, relation))
                edges.add((source, "variants", relation))
    # Add edges
    for source, target, relation in edges:
        element = {
            "data": {
                "source": source,
                "target": target,
            },
            "classes": relation.name
        }
        if element in elements:
            continue
        elements.append(element)
    # Setup graph webpage
    cyto.load_extra_layouts()
    app = Dash(__name__)
    # Show graph
    app.layout = layout_graph(elements, nodes, edges, default_layout=default_layout) 
    # Add interactive component callbacks 
    interactive_graph(app, elements, edges, auto_download=auto_download)
    # Start webpage
    app.run_server(debug=True)