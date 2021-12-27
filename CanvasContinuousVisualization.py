from collections import defaultdict
from mesa.visualization.ModularVisualization import VisualizationElement


class CanvasContinuous(VisualizationElement):
    """Creation of Custom VisualizationElement to render ContinuousSpace (CanvasContinuous)

        Mesa Package did not come with a VisualizationElement preloaded for its ContinuousSpace model, so I designed a custom one. Included are three modified JavaScript files that are
        designed to support this VisualizationElement and render the model onto the ModularServer. These JS Files are necessary to run the model.

        Space_Width and Space_Height (similar to grid_width and grid_height) are used to scale agents in the ContinuousSpace to proportionally visualize on the HTML5 canvas element.
        """

    package_includes = ["ContinuousDraw.js", "CanvasContinuousModule.js", "ContinuousInteractionHandler.js"]

    def __init__(self, portrayal_method, space_width, space_height, canvas_width=500, canvas_height=500):
        self.portrayal_method = portrayal_method
        self.space_width = space_width
        self.space_height = space_height
        self.canvas_width = canvas_width
        self.canvas_height = canvas_height
        new_element = "new ContinuousCanvasModule({}, {}, {}, {})".format(self.canvas_width, self.canvas_height, self.space_width, self.space_height)
        self.js_code = "elements.push(" + new_element + ");"




    def render(self , model):
        space_state = defaultdict(list)
        for agent in model.space._agent_to_index:
            portrayal = self.portrayal_method(agent)
            portrayal["x"] = agent.pos[0]
            portrayal["y"] = agent.pos[1]
            space_state[portrayal["Layer"]].append(portrayal)
        return space_state

