from collections import defaultdict
from mesa.visualization.ModularVisualization import VisualizationElement

class CanvasContinuous(VisualizationElement):



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

