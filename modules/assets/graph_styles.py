# Define default style
# Colour codes from https://m2.material.io/design/color/the-color-system.html#tools-for-picking-colors
selection_color = '#1E88E5' # Light blue 600
adj_color = "#4FC3F7" # Light blue 300


default_stylesheet = [
    # General style of types
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
            "font-size": "3px",
            "width": "30px", 
            "height": "10px",
            "shape": "round-rectangle",
            "text-wrap": "ellipsis",
            "text-max-width": "25px",
        }   
    }, {
        "selector": ".variant.observed",
        "style": {
            "shape": "hexagon",
        }
    }, {
        "selector": ".sample",
        "style": {
            "font-size": "10px",
            "width": "75px", 
            "height": "75px",
            "shape": "hexagon",
        }  
    # Style for different relations
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

# Style for different functional annotations 
function_colours = (
    ('no function', '#F44336'), # Red 500
    ('decreased function', '#FFEB3B'), # Yellow 500
    ('normal function', '#4CAF50'), # Green 500
    ('function not assigned', '#E0E0E0'), # Gray 300
    ('unknown function', '#BDBDBD'), # Gray 400
    ('uncertain function', '#9E9E9E'), # Gray 500
)
for function, colour in function_colours:
    default_stylesheet.append({
        "selector": f"node[function = '{function}']",
        "style": {
            "border-width": 2,
            "border-color": colour,
        }
    })

# Style for impact
impact_colours = (
    (0, '#BDBDBD'), # Gray 400
    (1, '#4CAF50'), # Green 500
    (2, '#FFEB3B'), # Yellow 500
    (3, '#F44336'), # Red 500
)
for severity, colour in impact_colours:
    default_stylesheet.append({
        "selector": f"node[severity = {severity}]",
        "style": {
            "border-width": 1,
            "border-color": colour,
        }
    })

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
                    "curve-style": "taxi",
                    "taxi-direction": "vertical",
                    "line-color": adj_color,
                    'target-arrow-color': adj_color,
                    'opacity': 1,
                    'z-index': 5000
                }
            }]
    return stylesheet