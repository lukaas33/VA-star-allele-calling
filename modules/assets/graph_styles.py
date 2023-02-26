# Define styles
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


def connected_styles(edge, direction):
    connected_color = '#D32F2F' # Red 700
    edge_color = '#E57373' # Red 300
    return [
        {
            "selector": f"node[id = '{edge[direction]}']",
            "style": {
                'background-color': connected_color,
                'opacity': 1
            }
        }, {
            "selector": f"edge[id = '{edge['id']}']",
            "style": {
                "line-color": edge_color,
                'target-arrow-color': edge_color,
                'opacity': 1,
                'z-index': 5000
            }
        }
    ]

def selection_stylesheet(node):
    selection_color = '#B71C1C' # Red 900
    stylesheet = default_stylesheet.copy()
    stylesheet += [
        {
            "selector": 'node',
            'style': {
                'opacity': 0.3
            }
        }, {
            'selector': 'edge',
            'style': {
                'opacity': 0.2
            }
        }, {
            "selector": f"node[id = '{node['data']['id']}']",
            "style": {
                'background-color': selection_color,
                "border-color": selection_color,
                "opacity": 1,
                'z-index': 9999
            }
        }
    ]

    for edge in node['edgesData']:
        if edge['source'] == node['data']['id']:
            stylesheet += connected_styles(edge, 'target')
        if edge['target'] == node['data']['id']:
            stylesheet += connected_styles(edge, 'source')
    return stylesheet