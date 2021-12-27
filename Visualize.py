from mesa.agent import Agent
from mesa.visualization.modules.CanvasContinuousVisualization import CanvasContinuous
from mesa.visualization.ModularVisualization import ModularServer
from StemCellABM import StemCell , Morphogen , Nodal , Lefty , ABM
import Constants


def agent_portrayal(agent : Agent):
    """Definition of Agent Portrayals

        Can be adjusted to make images of agents appear differently in the visualization. 
        Currently, I am using small circles for all with different colors to distinguish between agents

        Colors:
            StemCell = Black
            Morphogen = Orange
            Nodal = Blue
            Lefty = Purple"""

    if agent.__class__ == StemCell:
        portrayal = {
            "Shape" : "circle",
            "Filled" : "true",
            "r" : Constants.STEMCELL_R,
            "Color" : "black",
            "Layer" : 1
        }
    elif agent.__class__ == Morphogen:
        portrayal = {
            "Shape" : "circle",
            "Filled" : "true",
            "r" : Constants.MORPHOGEN_R,
            "Color" : "orange",
            "Layer" : 0
        }
    elif agent.__class__ == Nodal:
        portrayal = {
            "Shape" : "circle",
            "Filled" : "true",
            "r" : Constants.NODAL_R,
            "Color" : "blue",
            "Layer" : 2
        }
    elif agent.__class__ == Lefty:
        portrayal = {
            "Shape" : "circle",
            "Filled" : "true",
            "r" : Constants.LEFTY_R,
            "Color" : "purple",
            "Layer" : 3
        }
    return portrayal

#Create the ContinuousSpace
space = CanvasContinuous(agent_portrayal , Constants.MAX_X , Constants.MAX_Y , 500 , 500)
#Create the ModularServer
server = ModularServer(ABM , [space] , "Stem Cell ABM" , {"num_stem_cells": Constants.NUM_STEM_CELLS , "num_morph": Constants.NUM_MORPH , "sauce": Constants.SAUCE , "num_nodals": Constants.NUM_NODALS , "num_leftys": Constants.NUM_LEFTYS , "spawn_freq": Constants.SPAWN_FREQ , "diff_timer": Constants.DIFF_TIMER , "endo_min": Constants.ENDO_MIN , "ecto_max": Constants.ECTO_MAX , "max_x": Constants.MAX_X , "max_y": Constants.MAX_Y})
server.port = 8521
#Start Server
server.launch()