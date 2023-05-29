# Define default style
# Colour codes from https://m2.material.io/design/color/the-color-system.html#tools-for-picking-colors
default_color = '#455A64' # Blue gray 700
calling_color = "#03a9f4" # Blue 500

default_stylesheet = []
# General style of types
default_stylesheet.append({
    'selector': 'node',
    'style': {
        'content': 'data(label)',
        'background-color': default_color, 
        "text-valign": "center",
        "text-halign": "center",
        "color": "white",
        "font-family": "times new roman",
        "border-width": "2px",
        "border-color": default_color
    }
})
default_stylesheet.append({
    'selector': 'edge',
    'style': {
        'curve-style': 'bezier' 
    }
})
default_stylesheet.append({
    "selector": "node[?homozygous]",
    "style": {
        "border-style": "double",
        "border-width": "6px",
    }
})
default_stylesheet.append({
    "selector": ".core",
    "style": {
        "font-size": "15px",
        "width": "45px",
        "height": "45px",
        "shape": "ellipse",
    }
 })
default_stylesheet.append({
    "selector": ".sub",
    "style": {
        "font-size": "10px",
        "width": "35px",
        "height": "35px",
        "shape": "ellipse"
    }
})
default_stylesheet.append({
    "selector": ".called",
    "style": {
        "overlay-color": calling_color,
        "overlay-padding": "10px",
        "overlay-opacity": "0.15",
    }
})
default_stylesheet.append({
    "selector": ".variant",
    "style": {
        "font-size": "10px",
        "width": "80px", 
        "height": "20px",
        "shape": "rectangle",
        "text-wrap": "ellipsis",
        "text-max-width": "80px"
    }   
})
default_stylesheet.append({
    "selector": ".variant.personal",
    "style": {
        "shape": "tag",
    }
})
default_stylesheet.append({
    "selector": ".group",
    "style": {
        "font-size": "10px",
        "color": "black",
        "text-wrap": "wrap",
        "text-justification": "left",
        "background-color": "#b0bec5",
        "border-width": 0,
        "text-max-width": "35px",
        "height": "15px",
        "text-valign": "bottom",
        "text-halign": "right",
        "text-margin-y": "-10px",
    }
})
# TODO change height?
default_stylesheet.append({
    "selector": ".group.var",
    "style": {
        "shape": "square",
        "width": "75px",
        "text-margin-x": "-70px",
    }
})
default_stylesheet.append({
    "selector": ".group.p_var",
    "style": {
        "shape": "tag",
        "width": "120px",
        "text-margin-x": "-115px",
    }
})
default_stylesheet.append({
    "selector": ":compound",
    "style": {
        "width": "10px",
        "background-color": "white",
    }
})
default_stylesheet.append({
    "selector": ".sample",
    "style": {
        "font-size": "10px",
        "width": "75px", 
        "height": "35px",
        "shape": "hexagon",
    }  
})
# Style for different relations
default_stylesheet.append({
    'selector': '.EQUIVALENT',
    'style': {
        'line-color': '#37474F', # Blue gray 800
        'width': '2',
    }
})
default_stylesheet.append({
    'selector': '.IS_CONTAINED',
    'style': {
        'target-arrow-shape': 'triangle',
        'target-arrow-color': '#78909C', # Blue gray 400
        'line-color': '#78909C',
        'width': '1'
    }
})
default_stylesheet.append({
    'selector': '.CONTAINS',  
    'style': {
        # Only display one arbitrary direction of containment, 
        # the other direction follows
        'display': 'none'
    }     
})
default_stylesheet.append({
    'selector': '.OVERLAP',
    'style': {
        'line-color': '#90A4AE', # Blue gray 300
        'width': '1',
        'line-style': 'dashed'
    }
})
default_stylesheet.append({
    'selector': '.DISJOINT',
    'style': {
        'display': 'none'
    }
})


# Style for different functional annotations 
function_colours = (
    ('no function', '#F44336'), # Red 500
    ('decreased function', '#ffc107'), # Orange 500
    ('normal function', '#4CAF50'), # Green 500
    ('function not assigned', '#E0E0E0'), # Gray 300
    ('unknown function', '#BDBDBD'), # Gray 400
    ('uncertain function', '#9E9E9E'), # Gray 500
)
for function, colour in function_colours:
    default_stylesheet.append({
        "selector": f"node[function = '{function}']",
        "style": {
            "border-color": colour,
        }
    })

# Style for impact
impact_colours = (
    (0, '#BDBDBD'), # Gray 400
    (1, '#4CAF50'), # Green 500
    (2, '#ffc107'), # Orange 500
    (3, '#F44336'), # Red 500
)
for severity, colour in impact_colours:
    default_stylesheet.append({
        "selector": f"node[severity = {severity}]",
        "style": {
            "border-color": colour,
        }
    })

# Style for relevance
default_stylesheet.append({
    # TODO how to display this
    # TODO don't display?
    "selector": "node[!relevant]",
    "style": {
        "opacity": 0.8, 
    }
})

# Selection style
def selection_stylesheet(nodes, layout):
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
                "opacity": 0.9,
                'z-index': 9999
            }
        }]
    for source in nodes:
        for target in nodes:
            style = {
                "selector": f"edge[source = '{source}'][target = '{target}']",
                "style": {
                    'opacity': 1,
                    'z-index': 5000
                }
            }
            stylesheet.append(style)
    return stylesheet