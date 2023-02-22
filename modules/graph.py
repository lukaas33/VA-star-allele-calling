import algebra as va
import json
from dash import Dash, html, dcc, Input, Output
import dash_cytoscape as cyto
import plotly.express as px
from pandas import DataFrame as df
from .va_tools import count_arity
from .assets.graph_styles import default_stylesheet, selection_stylesheet
from .relations import prune_relations

def plot_arity(nodes, relations):
    """Create a plot of the arity values"""
    arity = count_arity(nodes, relations)
    # Covert data to pandas data frame
    data = {"allele": [], "relation": [], "arity": []}
    for rel in va.Relation:
        for node in nodes:
            data["allele"].append(node)
            data["relation"].append(str(rel).split('.')[1])
            data["arity"].append(arity[node][rel])
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
    This framework should not be used to make a complete visualization since it is limited in functionality.
    """
    edges = prune_relations(nodes, relations)
    return
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
            "classes": relation
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