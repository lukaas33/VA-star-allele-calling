import algebra as va
import json
from dash import Dash, html, dcc, Input, Output
import dash_cytoscape as cyto
import plotly.express as px
import pandas
from .va_tools import count_arity, count_relations
from .assets.graph_styles import default_stylesheet, selection_stylesheet
from .relations import prune_relations

# TODO add search function
# TODO add subgraph function

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

def layout_graph(elements, nodes, edges, relations):
    """Returns the layout for the Dash graph"""
    default_layout = 'cose-bilkent' 
    return html.Div([
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
                    dcc.Graph(
                        id="plot-arity",
                        figure=plot_arity(nodes, edges)
                    ),
                    dcc.Graph(
                        id="relation-counts",
                        figure=plot_relations(relations)
                    ),
                    dcc.Graph(
                        id="pruned-relation-counts",
                        figure=plot_relations(edges, pruned=True)
                    ),                    
                    dcc.Graph(
                        id="counts",
                        figure=plot_counts(elements)
                    )
                ]
            )
        ])        
    ])

def interactive_graph(app, original_elements):
    """Add interactive components to graph"""
    @app.callback(Output('graph', 'layout'), Input('change-layout', 'value'))
    def update_layout(layout):
        # Change layout
        settings = {"name": layout}
        settings["nodeDimensionsIncludeLabels"] = True
        if layout == 'cose-bilkent' or layout == 'cose':
            settings["idealEdgeLength"] = 250
            settings["tile"] = False
            settings["animate"] = False
        return settings
    @app.callback(Output('data', 'children'), Input('graph', 'tapNodeData'))
    def displayTapNodeData(data):
        # Display information about selection
        return json.dumps(data, indent=2)
    @app.callback([Output('graph', 'stylesheet'), Output('reset-selection', 'n_clicks')], [Input('graph', 'tapNode'), Input('reset-selection', 'n_clicks')])
    def generate_stylesheet(node, n_clicks):
        # Display connections of selected
        if not node or n_clicks == 1: # No input or resetting
            return [default_stylesheet, None]
        return [selection_stylesheet(node), None]
    @app.callback([Output('graph', 'zoom'), Output('graph', 'elements')], [Input('reset-view', 'n_clicks')])
    def reset_layout(n_clicks):
        # Reset view
        return [1, original_elements]
    @app.callback(
        Output("graph", "generateImage"),
        [
            Input("image-svg", "n_clicks"),
        ])
    def get_image(n_clicks):
        # Export image
        if n_clicks is None:
            return { # TODO why needed?
                'type': 'png',
                'action': 'store'
            }
        return {
            'type': 'svg',
            'action': 'download'
            }

def display_graph(relations, data):
    """Display relations as a graph

    Uses dash Cytoscape which creates a localhost website.
    The underlying framework is Cytoscape.js, a standard tool in biological network visualization.
    https://dash.plotly.com/cytoscape
    This framework should not be used to make a complete visualization since it is limited in functionality.
    """
    # TODO why does it call this function twice
    nodes, edges = prune_relations(relations)
    # Convert to proper format for cytoscape
    elements = []
    for node in nodes:
        if "*" in node:
            label = "*" + node.split("*")[1]
            info = data[node]
            if "." in node:
                category = "sub"
            else:
                category = "core"
        else:
            label = node.split(':')[1].split('.')[1]
            info = None # TODO fix data for this type
            category = "variant"
        elements.append({            
            "data": {
                "id": node, 
                "label": label,
                "data": info
            },
            "classes": category
        })

    for node, other, relation in edges:
        elements.append({
            "data": {
                "source": node,
                "target": other,
            },
            "classes": relation.name
        })
    # Setup graph webpage
    cyto.load_extra_layouts()
    app = Dash(__name__)
    # Show graph
    app.layout = layout_graph(elements, nodes, edges, relations)
    # Add interactive component callbacks 
    interactive_graph(app, elements)
    # Start webpage
    app.run_server(debug=True)