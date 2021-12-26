from mesa.agent import Agent
from mesa.visualization.modules.CanvasContinuousVisualization import CanvasContinuous
from mesa.visualization.ModularVisualization import ModularServer
from StemCellABM import StemCell , Morphogen , Nodal , Lefty , ABM

def agent_portrayal(agent : Agent):
    if agent.__class__ == StemCell:
        portrayal = {
            "Shape" : "circle",
            "Filled" : "true",
            "r" : 0.1,
            "Color" : "black",
            "Layer" : 1
        }
    elif agent.__class__ == Morphogen:
        portrayal = {
            "Shape" : "circle",
            "Filled" : "true",
            "r" : 0.05,
            "Color" : "orange",
            "Layer" : 0
        }
    elif agent.__class__ == Nodal:
        portrayal = {
            "Shape" : "circle",
            "Filled" : "true",
            "r" : 0.01,
            "Color" : "blue",
            "Layer" : 2
        }
    elif agent.__class__ == Lefty:
        portrayal = {
            "Shape" : "circle",
            "Filled" : "true",
            "r" : 0.01,
            "Color" : "purple",
            "Layer" : 3
        }
    return portrayal

space = CanvasContinuous(agent_portrayal , 20 , 20 , 500 , 500)
server = ModularServer(ABM , [space] , "Stem Cell ABM" , {"num_stem_cells": 10 , "num_morph": 10 , "sauce": False , "num_nodals": 10 , "num_leftys": 10 , "spawn_freq": 10 , "diff_timer": 10 , "endo_min": 10 , "ecto_max": 10})
server.port = 8521
server.launch()