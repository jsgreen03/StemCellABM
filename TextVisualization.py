"""
Text Module
============

Module for drawing live-updating text.

"""
from mesa.visualization.ModularVisualization import VisualizationElement
from mesa.visualization.TextVisualization import TextData



class TextElement(VisualizationElement):
    package_includes = ["TextModule.js"]
    js_code = "elements.push(new TextModule());"


    def __init__(self , attribute):
        self.js_code = "elements.push(new TextModule());"
        self.attribute = attribute

    def render(self, model):
        return TextData(model , self.attribute).render()


        
        





