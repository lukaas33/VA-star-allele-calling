# Define default style
selection_color = "#6a1b9a" # Purple 800
adj_color = '#B71C1C' # Red 900


default_stylesheet = [
    {
        'selector': 'node',
        'style': {
            'content': 'data(label)',
            'background-color': '#455A64',# Blue gray 700
            "text-valign": "center",
            "text-halign": "center",
            "color": "white",
        }
    }, {
        'selector': 'edge',
        'style': {
            'curve-style': 'bezier'
        }
    }, {
        "selector": "node:selected",
        "style": {
            "background-color": selection_color,
        }
    }, {
        "selector": ".core",
        "style": {
            "font-size": "15px",
            "width": "40px",
            "height": "40px",
            "shape": "ellipse"
        }
    }, {
        "selector": ".sub",
        "style": {
            "font-size": "5px",
            "width": "25px",
            "height": "25px",
            "shape": "ellipse"
        }
    }, {
        "selector": ".variant",
        "style": {
            "font-size": "5px",
            "width": "65px", 
            "height": "15px",
            "shape": "round-rectangle",
        }   
    }, {
        'selector': '.EQUIVALENT',
        'style': {
            'line-color': '#37474F', # Blue gray 800
            'width': '3'
        }
    }, {
        'selector': '.IS_CONTAINED',
        'style': {
            'target-arrow-shape': 'triangle',
            'target-arrow-color': '#78909C', # Blue gray 400
            'line-color': '#78909C',
            'width': '1'
        }
    }, {
        'selector': '.CONTAINS',  
        'style': {
            # Only display one arbitrary direction of containment, 
            # the other direction follows
            'display': 'none'
        }     
    }, {
        'selector': '.OVERLAP',
        'style': {
            'line-color': '#90A4AE', # Blue gray 300
            'width': '1',
            'line-style': 'dashed'
        }
    }, {
        'selector': '.DISJOINT',
        'style': {
            'display': 'none'
        }
    }
]

# Selection style
def selection_stylesheet(nodes):
    stylesheet = default_stylesheet.copy()
    stylesheet += [{
        "selector": 'node:unselected',
        'style': {
            'opacity': 0.3
        }
    }, {
        'selector': 'edge',
        'style': {
            'opacity': 0.2
            }
    }]
    for node in nodes:
        stylesheet += [{
            "selector": f"node[id = '{node}']:unselected",
            "style": {
                'background-color': adj_color,
                "opacity": 1,
                'z-index': 9999
            }
        }]
    for source in nodes:
        for target in nodes:
            stylesheet += [{
                "selector": f"edge[source = '{source}'][target = '{target}']",
                "style": {
                    "line-color": adj_color,
                    'target-arrow-color': adj_color,
                    'opacity': 1,
                    'z-index': 5000
                }
            }]
    return stylesheet