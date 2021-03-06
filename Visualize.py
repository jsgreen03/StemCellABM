from mesa.agent import Agent
from mesa.visualization.modules.CanvasContinuousVisualization import CanvasContinuous
from mesa.visualization.ModularVisualization import ModularServer
from mesa.visualization.modules.TextVisualization import TextElement
from StemCellABM import StemCell, BMP4 , NOG , ABM
import Constants


def agent_portrayal(agent : Agent):
    """Definition of Agent Portrayals

        Can be adjusted to make images of agents appear differently in the visualization. 
        Currently, I am using small circles for all with different colors to distinguish between agents

        Colors:
            StemCell = Black* 
            BMP4 = Blue
            NOG = Purple
            
            *Differentiated StemCells have color dependent on their classification:
                endo: Blue
                meso: Green
                ecto: Red 
            """

    if agent.__class__ == StemCell:
        dumb = 1
        if agent.differentiated == "virgin":
            color = "black"
        elif agent.differentiated == "endo":
            color = "blue"
        elif agent.differentiated == "meso":
            color = "green"
        elif agent.differentiated == "ecto":
            color = "red"
        elif agent.differentiated == "dumb":
            dumb = 0
            color = "white"
        portrayal = {
            "Shape" : "circle",
            "Filled" : "true",
            "r" : Constants.STEMCELL_R,
            "Color" : color,
            "Layer" : dumb
        }
        
    elif agent.__class__ == BMP4:
        portrayal = {
            "Shape" : "circle",
            "Filled" : "true",
            "r" : Constants.BMP4_R,
            "Color" : "blue",
            "Layer" : 2
        }
    elif agent.__class__ == NOG:
        if agent.absorbed:
                portrayal = {
                "Shape" : "circle",
                "Filled" : "true",
                "r" : Constants.NOG_R,
                "Color" : "white",
                "Layer" : 0
            }
        else:
            portrayal = {
                "Shape" : "circle",
                "Filled" : "true",
                "r" : Constants.NOG_R,
                "Color" : "purple",
                "Layer" : 3
            }
    return portrayal

#Create the ContinuousSpace
space = CanvasContinuous(agent_portrayal , Constants.MAX_X , Constants.MAX_Y , 500 , 500)

#Create TextElements for Model Paramaters
displayNumStemCells = TextElement("num_stem_cells")
displaySauce = TextElement("sauce")
displayNumBMP4 = TextElement("num_BMP4")
displayNumNOG = TextElement("num_NOG")
displaySpawFreq = TextElement("spawn_freq")
displayDiffTimer = TextElement("diff_timer")
displayEndoMin = TextElement("endo_min")
displayEctoMax = TextElement("ecto_max")
displayOneCellsContact = TextElement("one_cells_contact")
displayStartDiff = TextElement("start_diff")
displayStemCellEx = TextElement("stem_cell_ex")
displayStemCellExDiff = TextElement("stem_cell_ex_diff")
displayAvgX = TextElement("avg_x")
displayAvgY = TextElement("avg_y")
displayAvgRadius = TextElement("avg_radius")
displayMConcX = TextElement("mConcX")
displayMConcY = TextElement("mConcY")
displayMR = TextElement("mR")




#Create the ModularServer
server = ModularServer( ABM , 
                        [space , displayAvgX , displayAvgY , displayAvgRadius , displayEndoMin , displayEctoMax , displayNumStemCells , displaySauce , displayNumBMP4 , displayNumNOG , displaySpawFreq , 
                        displayDiffTimer ,  displayOneCellsContact , displayStartDiff , displayStemCellEx , 
                        displayStemCellExDiff , displayMConcX , displayMConcY , displayMR] ,
                         "Stem Cell ABM" , {"num_stem_cells": Constants.NUM_STEM_CELLS , "sauce": Constants.SAUCE , 
                         "num_BMP4": Constants.NUM_BMP4 , "num_NOG": Constants.NUM_NOG , "spawn_freq": Constants.SPAWN_FREQ , "diff_timer": Constants.DIFF_TIMER , 
                         "endo_min": Constants.ENDO_MIN , "ecto_max": Constants.ECTO_MAX , "max_x": Constants.MAX_X , "max_y": Constants.MAX_Y})
server.port = 8521
#Start Server
server.launch()