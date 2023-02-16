# Define styles
# TODO use colors/symbols and add legend 
default_stylesheet = [
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
]


def selection_stylesheet(node):
    selection_color = 'purple'
    connected_color = 'red'
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
                "color": selection_color,
                'z-index': 9999
            }
        }
    ]
    for edge in node['edgesData']:
        if edge['source'] == node['data']['id']:
            stylesheet.append({
                "selector": f"node[id = '{edge['target']}']",
                "style": {
                    'background-color': connected_color,
                    "color": connected_color,
                    'opacity': 0.9
                }
            })
            stylesheet.append({
                "selector": f"edge[id = '{edge['id']}']",
                "style": {
                    "line-color": connected_color,
                    'opacity': 0.9,
                    'z-index': 5000
                }
            })
        if edge['target'] == node['data']['id']:
            stylesheet.append({
                "selector": f"node[id = '{edge['source']}']",
                "style": {
                    'background-color': connected_color,
                    "color": connected_color,
                    'opacity': 0.9,
                    'z-index': 9999
                }
            })
            stylesheet.append({
                "selector": f"edge[id = '{edge['id']}']",
                "style": {
                    "line-color": connected_color,
                    'opacity': 1,
                    'z-index': 5000
                }
            })
    return stylesheet