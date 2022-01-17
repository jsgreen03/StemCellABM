from mesa import Model, Agent
from mesa.time import BaseScheduler, StagedActivation
from mesa.datacollection import DataCollector
from mesa.space import ContinuousSpace
from mesa.visualization.TextVisualization import TextData
import math
import copy
import Constants

class ABM(Model):
    '''Creation of Agent-Based-Model Class
        Stages of Model Per Tick: "movement" , "reaction_regulation" , "differentiation_tick" , "tracking_update"
        Global Variables:
            num_stem_cells : Int
            num_morph : Int
            sauce : Boolean : True if model will measure cell differentiation after a certain diff_timer elapses
            num_nodals : Int
            num_leftys : Int
            spawn_freq : Int : How much energy a given StemCell needs to reproduce
            diff_timer : Int
            endo_min : Int : How much concentration of Nodal a Stem Cell needs to come into contact to at minimum differntiate into an endoderm cell
            ecto_max : Int : How little concentration of Nodal a Stem Cell needs to come into contact to at maximum differntiate into an ectoderm cell
            one_cells_contact : Int : This attribute is used to test, it is designed to take an arbitraty Stem Cell each step and access its chemical contact
            start_diff : Boolean : When this is True, start positional differentiation if other Stem Cells have already differentiated nearby
            stem_cell_ex : StemCell : This attribute is used to latch onto a StemCell for purposes of tracking
            stem_cell_ex_diff : Boolean : This attribute indicates the example StemCell's differentiation status
            avg_x : Float : Indicates the Average X Value of all Stem Cells
            avg_y : Float : Indicates the Average Y Value of all Stem Cells
            avg_radius : Float : Indicates the Average Distance of all Stem Cells to the Centroid
            schedule : StagedActivation : Staged Activation Schedule that Follows the Stages listed above
            running : True : Batch will continually run this model's steps indefinitely'''

    def __init__(self, num_stem_cells: int , num_morph: int , sauce: bool , num_nodals: int , num_leftys: int , spawn_freq: int , diff_timer: int , endo_min: int , ecto_max: int , max_x:int=20 , max_y:int=20) -> None:
        self.num_stem_cells = num_stem_cells
        self.num_morph = num_morph
        self.sauce = sauce
        self.num_nodals = num_nodals
        self.num_leftys = num_leftys
        self.spawn_freq = spawn_freq
        self.diff_timer = diff_timer
        self.endo_min = endo_min
        self.ecto_max = ecto_max
        self.one_cells_contact = 0
        self.start_diff = False
        self.stem_cell_ex = None
        self.stem_cell_ex_diff = None
        self.avg_x = 0
        self.avg_y = 0
        self.avg_radius = 0
        self.schedule = BaseScheduler(self)
        self.running = True
        self.space = ContinuousSpace(max_x , max_y , False , 0 , 0)
        self.center_pos = self.space.center
        self.currentIDNum = 0
        self.hasCells = True#False
        self.setup()
        
        

    def setup(self):
        #Add StemCells to the Space
        if self.hasCells == True:
            for i in range(self.num_stem_cells):
                self.currentIDNum += 1
                c = StemCell(self.currentIDNum , self)
                self.schedule.add(c)
                r = self.random.random() * 3
                theta = self.random.random() * 2 * math.pi
                x = r * math.cos(theta) + self.center_pos[0]
                y = r * math.sin(theta) + self.center_pos[1]
                self.space.place_agent(c , (x , y))
            self.stem_cell_ex = self.space._index_to_agent[self.random.randrange(0 , len(self.space._agent_points))]
            self.stem_cell_ex_diff = self.stem_cell_ex.differentiated


            

        #Add Morphogens to the Space
        for i in range(self.num_morph):
            self.currentIDNum += 1
            m = Morphogen(self.currentIDNum, self)
            self.schedule.add(m)

            r = self.random.random() * 3
            theta = self.random.random() * 2 * math.pi
            x = r * math.cos(theta) + self.center_pos[0]
            y = r * math.sin(theta) + self.center_pos[1]
            
            self.space.place_agent(m , (x , y))


        #Add Nodals to the Space
        for i in range(self.num_nodals):
            self.currentIDNum += 1
            n = Nodal(self.currentIDNum, self)
            self.schedule.add(n)
            r = self.random.random()
            theta = self.random.random() * 2 * math.pi
            x = r * math.cos(theta) + self.center_pos[0]
            y = r * math.sin(theta) + self.center_pos[1]
            self.space.place_agent(n , (x , y))

        #Add Leftys to the Space
        for i in range(self.num_leftys):
            self.currentIDNum += 1
            l = Lefty(self.currentIDNum, self)
            self.schedule.add(l)
            r = self.random.random()
            theta = self.random.random() * 2 * math.pi
            x = r * math.cos(theta) + self.center_pos[0]
            y = r * math.sin(theta) + self.center_pos[1]
            self.space.place_agent(l , (x , y))

        self.calcAvgs()


    def calcAvgs(self):
        x = 0
        y = 0
        r = 0
        for agent in self.space._agent_to_index:
            if type(agent) == StemCell:
                x += agent.pos[0]
                y += agent.pos[1]
        x = x / self.num_stem_cells
        y = y / self.num_stem_cells
        centroid = (x , y)
        for agent in self.space._agent_to_index:
            if type(agent) == StemCell:
                r += self.space.get_distance(centroid , agent.pos)
        r = r / self.num_stem_cells
        self.avg_x = x
        self.avg_y = y
        self.avg_radius = r


    def spawnCells(self):
        agentDict = copy.deepcopy(self.space._agent_to_index)
        for agent in agentDict:
            if type(agent) == StemCell:
                if agent.energy >= self.spawn_freq:
                    self.currentIDNum += 1
                    newCell = StemCell(self.currentIDNum , self)
                    self.schedule.add(newCell)
                    self.space.place_agent(newCell , agent.pos)
                    agent.energy = agent.energy // 2
                    self.num_stem_cells += 1
    

    def cascade(self):
        for agent in self.space._agent_to_index:
            if type(agent) == StemCell:
                neighbors = self.space.get_neighbors(agent.pos , agent.internalR , True)
                for neighbor in neighbors:
                    if type(neighbor) == StemCell and neighbor.differentiated != "virgin":
                        agent.differentiated = neighbor.differentiated
                        agent.time_for_diff = 0
                 #       if agent.chemical_contact < agent.model.endo_min and agent.chemical_contact >= agent.model.ecto_max:
                  #          agent.differentiated = "meso"
                   #     if agent.chemical_contact < agent.model.ecto_max:
                    #        agent.differentiated = "ecto"


    def updateParams(self):
        agents = self.space._agent_to_index.keys()
        StemCells = [x for x in agents if type(x) == StemCell]
        c = StemCells[self.random.randrange(0 , len(StemCells))]
        self.one_cells_contact = c.chemical_contact
        self.stem_cell_ex_diff = self.stem_cell_ex.differentiated










    def step(self):
        self.spawnCells()
        self.calcAvgs()
        self.updateParams()
        if self.start_diff == True:
            self.cascade()
        self.schedule.step()

    






class StemCell(Agent):
    '''Creation of Stem Cell Agent
            Attributes: 
                Differentiated: Boolean
                Chemical Contact: Int
                Energy: Int
                Time For Diff: Int'''

    def __init__(self, unique_id: int, model: Model) -> None:
        super().__init__(unique_id , model)
        self.differentiated = "virgin"
        self.chemical_contact = 0
        self.energy = 0
        self.time_for_diff = self.random.randrange(Constants.TIME_FOR_DIFF_UPPER - 10 , Constants.TIME_FOR_DIFF_UPPER + 1)
        self.internalR = Constants.STEMCELL_R


    def step(self):
        self.movement()
        self.differentiation_tick()



    def movement(self):
        self.energy += self.random.randrange(0 , 3)
        hComp = 0
        vComp = 0
        for agent in self.model.space._agent_to_index:
            if type(agent) == Morphogen:
                h = (agent.pos[0] - self.pos[0])
                v = (agent.pos[1] - self.pos[1])
                n = (h ** 2 + v ** 2) ** (0.5)
                hComp += h/n
                vComp += v/n
        norm = (hComp ** 2 + vComp ** 2) ** (0.5)
        newPos = (self.pos[0] + (hComp / norm) , self.pos[1] + (vComp / norm))

        head =  self.model.space.get_heading((self.model.avg_x , self.model.avg_y) , newPos)
        headMag = (head[0] ** 2 + head[1] ** 2) ** (0.5)
        newPos = (newPos[0] + head[0]/(1.5*headMag) , newPos[1] + head[1]/(1.5*headMag)) 



        self.model.space.move_agent(self , newPos)


    def isTouchingOtherCells(self , point):
        d = self.model.space.get_distance(self.pos , point)
        l = d - (2*self.internalR)
        if l <= 0:
            return True 
        return False
        




    def isTouching(self , other:Agent):
        d = self.model.space.get_distance(self.pos , other.pos)
        l = d - self.internalR - other.internalR
        if l <= 0:
            return True 
        return False







    def reaction_regulation(self):
        return

    def differentiation_tick(self):
        if self.time_for_diff > 0:
            self.time_for_diff -= 1
            touchingNodal = 0
            for agent in self.model.space._agent_to_index:
                if type(agent) == Nodal and agent.active == True:
                    if self.isTouching(agent):
                        touchingNodal += 1
            self.chemical_contact += touchingNodal
        else:
            if self.differentiated == "virgin":
                if self.model.start_diff == False:
                    self.model.start_diff = True
                if self.chemical_contact >= self.model.endo_min:
                    self.differentiated = "endo"
                if self.chemical_contact < self.model.endo_min and self.chemical_contact >= self.model.ecto_max:
                    self.differentiated = "meso"
                if self.chemical_contact < self.model.ecto_max:
                    self.differentiated = "ecto"
                

        return

    def tracking_update(self):    
        return    







class Morphogen(Agent):
    '''Creation of Morphogen Agent'''
    
    def __init__(self, unique_id: int, model: Model) -> None:
        super().__init__(unique_id, model)
        self.internalR = Constants.MORPHOGEN_R

    def step(self):
        return
    



class Nodal(Agent):
    '''Creation of Nodal Agent'''

    def __init__(self, unique_id: int, model: Model) -> None:
        super().__init__(unique_id, model)
        self.immobilized = False
        self.immobilized_timer = 0
        self.active = True
        self.active_timer = 0
        self.internalR = Constants.NODAL_R

    def step(self):
        self.movement()
        self.reaction_regulation()


    def movement(self):
        if self.immobilized == False:
            heading = self.model.space.get_heading(self.pos , self.model.center_pos)
            if heading[0] > 0 and heading[1] > 0:
                if self.random.randrange(0 , 2) < 1: #Vertical Shift

                    x = heading[0]
                    y = self.random.uniform(-heading[0] , heading[1])

                else: #Horizontal Shift

                    x = self.random.uniform(-heading[1] , heading[0])
                    y = heading[1]

            elif heading[0] < 0 and heading[1] < 0:

                if self.random.randrange(0 , 2) < 1: #Vertical Shift

                    x = heading[0]
                    y = self.random.uniform(heading[1] , -heading[0])

                else: #Horizontal Shift

                    x = self.random.uniform(heading[0] , -heading[1])
                    y = heading[1]

            elif heading[0] < 0 and heading[1] > 0:

                if self.random.randrange(0 , 2) < 1: #Vertical Shift

                    x = heading[0]
                    y = self.random.uniform(heading[0] , heading[1])

                else: #Horizontal Shift

                    x = self.random.uniform(heading[0] , heading[1])
                    y = heading[1]

            elif heading[0] > 0 and heading[1] < 0:

                if self.random.randrange(0 , 2) < 1: #Vertical Shift

                    x = heading[0]
                    y = self.random.uniform(heading[1] , heading[0])

                else: #Horizontal Shift

                    x = self.random.uniform(heading[1] , heading[0])
                    y = heading[1]

            elif heading == (0,0):
                
                x = 1
                y = 1


            norm = (x ** 2 + y ** 2) ** 0.5
            xDisplacement = x / (norm * 3)
            yDisplacement = y / (norm * 3)
            
            self.model.space.move_agent(self , (self.pos[0] + xDisplacement , self.pos[1] + yDisplacement))
            for agent in self.model.space._agent_to_index:
                if type(agent) == StemCell:
                    if self.isTouching(agent):
                        xDisplacement = -xDisplacement / 2
                        yDisplacement = -yDisplacement / 2
            self.model.space.move_agent(self , (self.pos[0] + xDisplacement , self.pos[1] + yDisplacement))
            if self.random.randrange(0 , 100) < 49:
                self.immobilized = True
                self.immobilized_timer = self.random.randrange(0 , 11)
        else:
            self.immobilized_timer -= 1
            if self.immobilized_timer == 0:
                self.immobilized = False
            
            


    def reaction_regulation(self):
        if self.active_timer == 0:
            for agent in self.model.space._agent_to_index:
                if type(agent) == Lefty:
                    if self.isTouching(agent):
                        self.active = False
                        self.active_timer = self.random.randrange(0 , 11)
        else:
            self.active_timer -= 1
            if self.active_timer == 0:
                self.active = True


 

    def isTouching(self , other:Agent):
        d = self.model.space.get_distance(self.pos , other.pos)
        l = d - self.internalR - other.internalR
        if l <= 0:
            return True 
        return False




class Lefty(Agent):
    '''Creation of Lefty Agent'''

    def __init__(self, unique_id: int, model: Model) -> None:
        super().__init__(unique_id, model)
        self.internalR = Constants.LEFTY_R


    def step(self):
        self.movement()


    def movement(self):
            heading = self.model.space.get_heading(self.pos , self.model.center_pos)
            if heading[0] > 0 and heading[1] > 0:
                if self.random.randrange(0 , 2) < 1: #Vertical Shift

                    x = heading[0]
                    y = self.random.uniform(-heading[0] , heading[1])

                else: #Horizontal Shift

                    x = self.random.uniform(-heading[1] , heading[0])
                    y = heading[1]

            elif heading[0] < 0 and heading[1] < 0:

                if self.random.randrange(0 , 2) < 1: #Vertical Shift

                    x = heading[0]
                    y = self.random.uniform(heading[1] , -heading[0])

                else: #Horizontal Shift

                    x = self.random.uniform(heading[0] , -heading[1])
                    y = heading[1]

            elif heading[0] < 0 and heading[1] > 0:

                if self.random.randrange(0 , 2) < 1: #Vertical Shift

                    x = heading[0]
                    y = self.random.uniform(heading[0] , heading[1])

                else: #Horizontal Shift

                    x = self.random.uniform(heading[0] , heading[1])
                    y = heading[1]

            elif heading[0] > 0 and heading[1] < 0:

                if self.random.randrange(0 , 2) < 1: #Vertical Shift

                    x = heading[0]
                    y = self.random.uniform(heading[1] , heading[0])

                else: #Horizontal Shift

                    x = self.random.uniform(heading[1] , heading[0])
                    y = heading[1]

            elif heading == (0,0):

                x = 1
                y = 1


            norm = (x ** 2 + y ** 2) ** 0.5
            xDisplacement = x / (norm * 0.5)
            yDisplacement = y / (norm * 0.5)
            
            self.model.space.move_agent(self , (self.pos[0] + xDisplacement , self.pos[1] + yDisplacement))


    def isTouching(self , other:Agent):
        d = self.model.space.get_distance(self.pos , other.pos)
        l = d - self.internalR - other.internalR
        if l <= 0:
            return True 
        return False







