from mesa import Model, Agent
from mesa.time import StagedActivation
from mesa.datacollection import DataCollector
from mesa.space import ContinuousSpace
import math

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
class ABM(Model):

    def __init__(self, num_stem_cells: int , num_morph: int , sauce: bool , num_nodals: int , num_leftys: int , spawn_freq: int , diff_timer: int , endo_min: int , ecto_max: int) -> None:
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
        self.stem_cell_ex_diff = False
        self.avg_x = 0
        self.avg_y = 0
        self.avg_radius = 0
        self.schedule = StagedActivation(self , ["movement" , "reaction_regulation" , "differentiation_tick" , "tracking_update"])
        self.running = True
        self.space = ContinuousSpace(20 , 20 , False , 0 , 0)
        self.setup()
        

    def setup(self):
        #Add StemCells to the Space
        for i in range(self.num_stem_cells):
            c = StemCell(i , self)
            self.schedule.add(c)
            x = self.random.randrange(self.space.x_min , self.space.x_max)
            y = self.random.randrange(self.space.y_min , self.space.y_max)
            self.space.place_agent(c , (x , y))

        #Add Morphogens to the Space
        for i in range(self.num_morph):
            m = Morphogen(i + self.num_stem_cells, self)
            self.schedule.add(m)
            x = self.random.randrange(self.space.x_min , self.space.x_max)
            y = self.random.randrange(self.space.y_min , self.space.y_max)
            self.space.place_agent(m , (x , y))


        #Add Nodals to the Space
        for i in range(self.num_nodals):
            n = Nodal(i + self.num_stem_cells + self.num_morph, self)
            self.schedule.add(n)
            x = self.random.randrange(self.space.x_min , self.space.x_max)
            y = self.random.randrange(self.space.y_min , self.space.y_max)
            self.space.place_agent(n , (x , y))

        #Add Leftys to the Space
        for i in range(self.num_leftys):
            l = Lefty(i + self.num_stem_cells + self.num_morph + self.num_nodals, self)
            self.schedule.add(l)
            x = self.random.randrange(self.space.x_min , self.space.x_max)
            y = self.random.randrange(self.space.y_min , self.space.y_max)
            self.space.place_agent(l , (x , y))



    def step(self):
        self.schedule.step()


'''Creation of Stem Cell Agent
Attributes: 
   Differentiated: Boolean
   Chemical Contact: Int
   Energy: Int
   Time For Diff: Int'''
class StemCell(Agent):

    def __init__(self, unique_id: int, model: Model) -> None:
        super().__init__(unique_id , model)
        self.differentiated = False
        self.chemical_contact = 0
        self.energy = 0
        self.time_for_diff = 0
        self.theta = 0

    def movement(self):
        r = 1
        x = self.pos[0]
        y = self.pos[1]
        self.theta = self.random.random() * 2 * math.pi
        x += r * math.cos(self.theta)
        y += r * math.sin(self.theta)
        if self.model.space.out_of_bounds((x,y)):
            if self.model.space.out_of_bounds((-x,-y)):
                if self.model.space.out_of_bounds((x,-y)):
                    if self.model.space.out_of_bounds((-x,y)):
                        pass
                    else:
                        self.model.space.move_agent(self , (-x , y))
                else:
                    self.model.space.move_agent(self , (x , -y))
            else:
                self.model.space.move_agent(self , (-x , -y))
        else:
            self.model.space.move_agent(self , (x , y))


    def reaction_regulation(self):
        return

    def differentiation_tick(self):
        return

    def tracking_update(self):    
        return    



'''Creation of Morphogen Agent'''
class Morphogen(Agent):
    
    def __init__(self, unique_id: int, model: Model) -> None:
        super().__init__(unique_id, model)

    def movement(self):
        return

    def reaction_regulation(self):
        return

    def differentiation_tick(self):
        return

    def tracking_update(self):    
        return    
    


'''Creation of Nodal Agent'''
class Nodal(Agent):

    def __init__(self, unique_id: int, model: Model) -> None:
        super().__init__(unique_id, model)
        self.immobilized = False
        self.immobilized_timer = 0
        self.active = True
        self.active_timer = 0

    def movement(self):
        r = 1
        x = self.pos[0]
        y = self.pos[1]
        self.theta = self.random.random() * 2 * math.pi
        x += r * math.cos(self.theta)
        y += r * math.sin(self.theta)
        if self.model.space.out_of_bounds((x,y)):
            if self.model.space.out_of_bounds((-x,-y)):
                if self.model.space.out_of_bounds((x,-y)):
                    if self.model.space.out_of_bounds((-x,y)):
                        pass
                    else:
                        self.model.space.move_agent(self , (-x , y))
                else:
                    self.model.space.move_agent(self , (x , -y))
            else:
                self.model.space.move_agent(self , (-x , -y))
        else:
            self.model.space.move_agent(self , (x , y))


    def reaction_regulation(self):
        return

    def differentiation_tick(self):
        return

    def tracking_update(self):    
        return    



'''Creation of Lefty Agent'''
class Lefty(Agent):

    def __init__(self, unique_id: int, model: Model) -> None:
        super().__init__(unique_id, model)

    def movement(self):
        r = 1
        x = self.pos[0]
        y = self.pos[1]
        self.theta = self.random.random() * 2 * math.pi
        x += r * math.cos(self.theta)
        y += r * math.sin(self.theta)
        if self.model.space.out_of_bounds((x,y)):
            if self.model.space.out_of_bounds((-x,-y)):
                if self.model.space.out_of_bounds((x,-y)):
                    if self.model.space.out_of_bounds((-x,y)):
                        pass
                    else:
                        self.model.space.move_agent(self , (-x , y))
                else:
                    self.model.space.move_agent(self , (x , -y))
            else:
                self.model.space.move_agent(self , (-x , -y))
        else:
            self.model.space.move_agent(self , (x , y))


    def reaction_regulation(self):
        return

    def differentiation_tick(self):
        return

    def tracking_update(self):    
        return    









