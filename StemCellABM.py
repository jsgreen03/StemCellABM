
from typing import DefaultDict
from matplotlib.pyplot import sca
from mesa import Model, Agent
from mesa.time import BaseScheduler, StagedActivation
from mesa.datacollection import DataCollector
from mesa.space import ContinuousSpace
from mesa.visualization.TextVisualization import TextData
import math
import Brandon
import copy
import Constants

class ABM(Model):
    '''Creation of Agent-Based-Model Class
        Stages of Model Per Tick: "movement" , "reaction_regulation" , "differentiation_tick" , "tracking_update"
        Global Variables:
            num_stem_cells : Int
            sauce : Boolean : True if model will measure cell differentiation after a certain diff_timer elapses
            num_BMP4 : Int
            num_NOG : Int
            spawn_freq : Int : How much energy a given StemCell needs to reproduce
            diff_timer : Int
            endo_min : Int : How much concentration of BMP4 a Stem Cell needs to come into contact to at minimum differntiate into an endoderm cell
            ecto_max : Int : How little concentration of BMP4 a Stem Cell needs to come into contact to at maximum differntiate into an ectoderm cell
            one_cells_contact : Int : This attribute is used to test, it is designed to take an arbitraty Stem Cell each step and access its chemical contact
            start_diff : Boolean : When this is True, start positional differentiation if other Stem Cells have already differentiated nearby
            stem_cell_ex : StemCell : This attribute is used to latch onto a StemCell for purposes of tracking
            stem_cell_ex_diff : Boolean : This attribute indicates the example StemCell's differentiation status
            avg_x : Float : Indicates the Average X Value of all Stem Cells
            avg_y : Float : Indicates the Average Y Value of all Stem Cells
            avg_radius : Float : Indicates the Average Distance of all Stem Cells to the Centroid
            schedule : StagedActivation : Staged Activation Schedule that Follows the Stages listed above
            running : True : Batch will continually run this model's steps indefinitely'''

    def __init__(self, num_stem_cells: int , sauce: bool , num_BMP4: int , num_NOG: int , spawn_freq: int , diff_timer: int , endo_min: int , ecto_max: int , max_x:int=20 , max_y:int=20) -> None:
        self.num_stem_cells = num_stem_cells
        self.sauce = sauce
        self.num_BMP4 = num_BMP4
        self.num_NOG = num_NOG
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
        self.hasCells = True
        self.cellTouchingDict = {}
        for i in range(1 , Constants.NUM_STEM_CELLS + 1):
            self.cellTouchingDict[i] = []
        self.end_time = 0
        self.mConcX = 0
        self.mConcY = 0
        self.mR = 0
        self.cells = []
        self.NOG = []
        self.BMP4 = []
        self.BMP4vector = Brandon.unot
        self.setup()
        
        

    def setup(self):
        #Add StemCells to the Space
        if self.hasCells == True:
            for i in range(self.num_stem_cells):
                self.currentIDNum += 1
                c = StemCell(self.currentIDNum , self)
                self.schedule.add(c)
                r = self.random.random() * 1
                theta = self.random.random() * 2 * math.pi
                x = r * math.cos(theta) + self.center_pos[0]
                y = r * math.sin(theta) + self.center_pos[1]
                self.space.place_agent(c , (x , y))
                self.cells.append(c)
            self.stem_cell_ex = self.space._index_to_agent[self.random.randrange(0 , len(self.space._agent_points))]
            self.stem_cell_ex_diff = self.stem_cell_ex.differentiated



        #Add BMP4 to the Space
        for i in range(self.num_BMP4):
            self.currentIDNum += 1
            n = BMP4(self.currentIDNum, self)
            self.schedule.add(n)
            r = self.random.random()
            theta = self.random.random() * 2 * math.pi
            x = r * math.cos(theta) + self.center_pos[0]
            y = r * math.sin(theta) + self.center_pos[1]
            self.space.place_agent(n , (x , y))

            self.BMP4.append(n)

        #Add NOG to the Space
        for i in range(self.num_NOG):
            self.currentIDNum += 1
            l = NOG(self.currentIDNum, self)
            self.schedule.add(l)
            r = self.random.random()
            theta = self.random.random() * 2 * math.pi
            x = r * math.cos(theta) + self.center_pos[0]
            y = r * math.sin(theta) + self.center_pos[1]
            self.space.place_agent(l , (x , y))

            self.NOG.append(l)

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



   
    

    def cascade(self):
        for agent in self.space._agent_to_index:
            if type(agent) == StemCell:
                neighbors = self.space.get_neighbors(agent.pos , agent.internalR , True)
                for neighbor in neighbors:
                    if type(neighbor) == StemCell and neighbor.differentiated == "virgin":
                        agent.differentiated = neighbor.differentiated
                        agent.time_for_diff = 0



    def updateParams(self):
        agents = self.space._agent_to_index.keys()
        StemCells = [x for x in agents if type(x) == StemCell]
        c = StemCells[self.random.randrange(0 , len(StemCells))]
        self.one_cells_contact = c.chemical_contact
        self.stem_cell_ex_diff = self.stem_cell_ex.differentiated



    def updateTouchingDict(self):
        #self.cellTouchingDict
        for key in self.cellTouchingDict.keys():
            self.cellTouchingDict[key] = []

        for agent in self.space._agent_to_index:

                if type(agent) == StemCell:
             
                    neighbors = self.space.get_neighbors(agent.pos , agent.internalR -.01 , include_center=False)
                    for neighbor in neighbors:
                        if type(neighbor) == StemCell:
                            if neighbor not in self.cellTouchingDict[agent.unique_id]:
                                self.cellTouchingDict[agent.unique_id].append(neighbor)
                            if agent not in self.cellTouchingDict[neighbor.unique_id]:
                                self.cellTouchingDict[neighbor.unique_id].append(agent)
                      

    def updateBMP4(self):
        self.BMP4vector = Brandon.reaction(self.BMP4)


    def step(self):
        self.calcAvgs()
        if self.hasCells:
            self.updateParams()
            self.updateTouchingDict()
            self.updateBMP4()
        if self.start_diff == True:
            self.cascade()
            if self.end_time == -2:
                self.end_time = 3
            self.end_time -= 1
        self.schedule.step()
        if self.end_time == -1 and self.start_diff == True:
            self.running = False

    






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
        self.absorbedNOG = []


    def step(self):
        self.spawnCells()
        self.movement2()
        if self.model.sauce == True:
            self.differentiation_tick()

    def movement2(self):
        scaleFactor = 5
        self.energy += self.random.randrange(0 , 3)
        if self.differentiated == "virgin":
            neighbors = self.model.cellTouchingDict[self.unique_id]
            neighborPoints = []
            for neighbor in neighbors:
                points = self.intersectingPoints(neighbor)
                if points[0] != self.pos:
                    neighborPoints.append(points)
            polarNeighborPoints = []
            for intersection in neighborPoints:
                newInter = self.convertIntersectingPointsPolar(intersection)
                polarNeighborPoints.append(newInter)
            nSets = []
            for intersection in polarNeighborPoints:
                nSets.append(SetRange(intersection[0][1] , intersection[1][1] , True , intersection[1]))
            neighborSet = SetRangeUnion(nSets)
            r = neighborSet.getRange()
            thetas0 = neighborSet.getPoints()
            thetas = []
            for pTuple in thetas0:
                thetas.append(pTuple[0])
                thetas.append(pTuple[1])
            centerDir = self.model.space.get_heading(self.pos , (self.model.avg_x , self.model.avg_y))
            magCenterDir = (centerDir[0] ** 2 + centerDir[1] ** 2)**0.5
            zpAngle = math.atan(centerDir[1] / centerDir[0])
            if centerDir[0] < 0 and centerDir[1] > 0:
                zpAngle += math.pi
            if centerDir[0] < 0 and centerDir[1] < 0:
                zpAngle += math.pi
            if centerDir[0] > 0 and centerDir[1] < 0:
                zpAngle += 2*math.pi
            zpAngle += self.random.random()*(math.pi)/magCenterDir
            while zpAngle > math.pi * 2:
                zpAngle -= math.pi*2
            if zpAngle < (math.pi / 2) or zpAngle > (3 * math.pi / 2):
                l = True
            else:
                l = False
            zpRange = SetRange(zpAngle + math.pi / 2 , zpAngle - math.pi / 2 , True , l)
            noThetas = True
            for theta in thetas:
                if zpRange.numInRange(theta):
                    noThetas = False
            if not noThetas:
                if neighborSet.numInRange(zpAngle):
                    zpAngle = self.getMinDistanceAlongCircumference(zpAngle , thetas)
            x = self.internalR * math.cos(zpAngle)
            y = self.internalR * math.sin(zpAngle)

            head1 = self.model.space.get_heading(self.pos , (self.pos[0] + x, self.pos[1] + y))
            head1 = (head1[0] *(((2*math.pi)-r)/(2*math.pi)) , head1[1]*(((2*math.pi)-r)/(2*math.pi)))
            h = head1[0]
            v = head1[1]
            norm = scaleFactor * (h ** 2 + v ** 2) ** (0.5)
            if h == 0 and v == 0:
                norm = 1
            newPos = (self.pos[0] + (h/norm) , self.pos[1] + (v/norm))              
            self.model.space.move_agent(self , newPos)



    def spawnCells(self):
        if self.energy >= self.model.spawn_freq:
            self.model.currentIDNum += 1
            newCell = StemCell(self.model.currentIDNum , self.model)
            self.model.schedule.add(newCell)
            self.model.space.place_agent(newCell , self.pos)
            self.energy = self.energy // 2
            self.model.num_stem_cells += 1
            self.model.cellTouchingDict[self.model.currentIDNum] = []
            self.model.cells.append(newCell)


#Does not work do not use 
#(No more generic morphogens are in this model)
    def BROKEN(self):
        scaleFactor = 0.9
        self.energy += self.random.randrange(0 , 3)
        if self.differentiated == "virgin":
            neighbors = self.model.cellTouchingDict[self.unique_id]
            neighborPoints = []
            for neighbor in neighbors:
                points = self.intersectingPoints(neighbor)
                if points[0] != self.pos:
                    neighborPoints.append(points)
            polarNeighborPoints = []
            for intersection in neighborPoints:
                newInter = self.convertIntersectingPointsPolar(intersection)
                polarNeighborPoints.append(newInter)
            nSets = []
            for intersection in polarNeighborPoints:
                nSets.append(SetRange(intersection[0][1] , intersection[1][1] , True , intersection[1]))
            neighborSet = SetRangeUnion(nSets)
            r = neighborSet.getRange()
            thetas0 = neighborSet.getPoints()
            thetas = []
            for pTuple in thetas0:
                thetas.append(pTuple[0])
                thetas.append(pTuple[1])
            hComp = 0
            vComp = 0
            for agent in self.model.morphs:
                if agent.dead == False:
                    head0 = self.model.space.get_heading(self.pos , agent.pos)
                    h = head0[0]
                    v = head0[1]
                    n = (h ** 2 + v ** 2) ** 0.5
                    zpChange = (h * self.internalR / n , v * self.internalR / n)
                    zeroPoint = (self.pos[0] + (zpChange[0]) , self.pos[1] + (zpChange[1]))
                    zpAngle = self.determineAngle(zeroPoint)
                    if zpAngle < (math.pi / 2) or zpAngle > (3 * math.pi / 2):
                        l = True
                    else:
                        l = False
                    zpRange = SetRange(zpAngle + math.pi / 2 , zpAngle - math.pi / 2 , True , l)
                    x = 0
                    y = 0
                    noThetas = True
                    for theta in thetas:
                        if zpRange.numInRange(theta):
                            noThetas = False
                    if not noThetas:
                        if neighborSet.numInRange(zpAngle):
                            zpAngle = self.getMinDistanceAlongCircumference(zpAngle , thetas)
                        x = self.internalR * math.cos(zpAngle)
                        y = self.internalR * math.sin(zpAngle)
                    head1 = self.model.space.get_heading(self.pos , (self.pos[0] + x, self.pos[1] + y))
                    head1 = (head1[0] *(((2*math.pi)-r)/(2*math.pi)) , head1[1]*(((2*math.pi)-r)/(2*math.pi)))
                    hComp += head1[0]
                    vComp += head1[1]
            norm = scaleFactor * (hComp ** 2 + vComp ** 2) ** (0.5)
            if hComp == 0 and vComp == 0:
                norm = 1
            newPos = (self.pos[0] + (hComp/norm) , self.pos[1] + (vComp/norm))              
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

    def spawnBMP4(self, num):
        for i in range(num):
            self.model.currentIDNum += 1
            n = NOG(self.model.currentIDNum, self.model)
            self.model.schedule.add(n)
            r = self.internalR + 0.001
            theta = self.model.random.random() * 2 * math.pi
            x = r * math.cos(theta) + self.pos[0]
            y = r * math.sin(theta) + self.pos[1]
            self.model.space.place_agent(n , (x , y))
        self.absorbedNOG = []


    def differentiation_tick(self):
        if self.time_for_diff > 0:
            self.time_for_diff -= 1
        else:
            matrixIndex = (self.pos[0] // Brandon.dx, self.pos[1] // Brandon.dx)
            vectorIndex = points * matrixIndex[0] + matrixIndex[1]
            BMP4conc = self.model.BMP4vector[vectorIndex]


            if self.differentiated == "virgin":
                if self.model.start_diff == False:
                    self.model.start_diff = True
                if self.BMP4conc >= self.model.endo_min:
                    self.differentiated = "endo"
                if self.BMP4conc < self.model.endo_min and self.chemical_contact >= self.model.ecto_max:
                    self.differentiated = "meso"
                if self.BMP4conc < self.model.ecto_max:
                    self.differentiated = "ecto"
                


    def tracking_update(self):    
        return    

    def intersectingPoints(self , other):
        a = self.pos[0]
        b = self.pos[1]
        c = other.pos[0]
        d = other.pos[1]
        r = self.internalR
        D = ((a-c)**2 + (b-d)**2)**0.5
        if d == b:
            x1 = min(a,c) + r * ( (2)**0.5 / 2)
            x2 = x1
            y1 = ((r)**2 - (r * ( (2)**0.5 / 2))**2) ** 0.5
            y2 = -1 * y1
        else:
            if D > 2*self.internalR:
                return [(self.pos[0] , self.pos[1])]
            E = (r**2 - (D/2)**2)**0.5
            M = (c-a)/(b-d)
            B = (b**2 - d**2)/(c**2 - a**2)
            theta = math.atan2(c-a , b-d)
            midX = (a+c)/2
            midY = (b+d)/2
            dist = E * math.cos(theta)
            x2 = midX + dist
            y2 = M*x2 + B
            x1 = midX - dist
            y1 = M*x1 + B
        val = False
        head = self.model.space.get_heading(self.pos , other.pos)
        if head[0] > 0 and b > min(y2,y1) and b < max(y2,y1):
            val = True            
        return [(x1,y1) , (x2,y2) , val]

    def convertIntersectingPointsPolar(self, points):
        polarCoords = []
       
        p = (points[0] , points[1])
        polar0 = [self.internalR]
        polar1 = [self.internalR]
        polar0.append(self.determineAngle(p[0]))
        polar1.append(self.determineAngle(p[1]))

        polarCoords.append(tuple(polar0))
        polarCoords.append(tuple(polar1))

        polarCoords.append(points[2])
        return polarCoords

    def determineAngle(self, point):
        radius = self.internalR
        center = self.pos
        if point == (center[0] , center[1] + radius):
            angle = math.pi / 2
        elif point == (center[0] + radius , center[1]):
            angle = 0
        elif point == (center[0] - radius , center[1]):
            angle = math.pi
        elif point == (center[0], center[1] - radius):
            angle = math.pi * 3 / 2

        else:
            x = point[0]
            y = point[1]
            theta = math.atan(y/x)
            if point[0] > center[0] and point[1] > center[1]:
                angle = theta
            elif point[0] < center[0] and point[1] > center[1]:
                angle = math.pi - theta
            elif point[0] > center[0] and point[1] < center[1]:
                angle = (2*math.pi) - theta
            elif point[0] < center[0] and point[1] < center[1]:
                angle = theta+(math.pi)
        return angle




    def getMinDistanceAlongCircumference(self, zeroPoint , thetas):
        distances = []
        if thetas == []:
            return zeroPoint
        else:
            for theta in thetas:
                distances.append(abs(theta - zeroPoint))
            angle = min(distances)
            corrIndex = distances.index(angle)
            return thetas[corrIndex]


        
        



class SetRange:

    def __init__(self , num1 , num2 , closed , lapped):
        self.lower = min(num1 , num2)
        self.upper = max(num1 , num2)
        self.closed = closed
        self.lapped = lapped

    def __lt__(self , other):
        return self.lower < other.lower


    def numInRange(self , num):
        if self.closed:
            if num == self.lower or num == self.upper:
                return True
        if not self.lapped:
            return (num > self.lower) and (num < self.upper)
        else:
            return (num < self.lower) or (num > self.upper)
 

    def getPoints(self):
        return (self.upper , self.lower)

    def getRange(self):
        if self.lapped:
            return (2*math.pi) - (self.upper - self.lower)
        return self.upper - self.lower

    def __eq__(self , other):
        return (self.upper == other.upper) and (self.lower == self.lower) and (self.closed == other.closed)

    def __str__(self) -> str:
        if self.closed:
            return "[{} , {}]".format(self.lower , self.upper)
        return "({} , {})".format(self.lower , self.upper)


class SetRangeUnion:
    
    def __init__(self , sets):
        setCopy = []
        for set in sets:
            if set.lapped == True:
                setCopy.append(SetRange(0 , set.lower , True , False))
                setCopy.append(SetRange(set.upper , math.pi*2 , True , False))
            else: 
                setCopy.append(set)

        self.sets = sorted(setCopy)
        if len(self.sets) > 1:
            i = 1
            new = []
            currSet = self.sets[0]
            while i <= len(self.sets) - 1:
                nextSet = self.sets[i]
                if currSet.upper > nextSet.lower:
                    if nextSet.upper > currSet.upper:
                        currSet.upper = nextSet.upper
                else: 
                    new.append(currSet)
                    currSet = nextSet
                i += 1
            new.append(currSet)
            self.sets = new


    def numInRange(self , num):
        for set in self.sets:
            val = set.numInRange(num)
            if val:
                return True
        return False

    def getRange(self):
        r = 0
        for s in self.sets:
            r += s.getRange()
        return r


    def getPoints(self):
        points = []
        for s in self.sets:
            points.append(s.getPoints())
        return points


class BMP4(Agent):
    '''Creation of BMP4 Agent'''

    def __init__(self, unique_id: int, model: Model) -> None:
        super().__init__(unique_id, model)
        self.immobilized = False
        self.immobilized_timer = 0
        self.active = True
        self.active_timer = 0
        self.internalR = Constants.BMP4_R

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
                if type(agent) == NOG:
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




#class NOG(Agent):
 #   '''Creation of NOG Agent'''

  #  def __init__(self, unique_id: int, model: Model) -> None:
   #     super().__init__(unique_id, model)
    #    self.internalR = Constants.NOG_R
     #   self.absorbed = False


    #def step(self):
     #   if not self.absorbed:
      #      self.movement()


    #def movement(self):
     #       heading = self.model.space.get_heading(self.pos , self.model.center_pos)
      #      if heading[0] > 0 and heading[1] > 0:
       #         if self.random.randrange(0 , 2) < 1: #Vertical Shift

        #            x = heading[0]
         #           y = self.random.uniform(-heading[0] , heading[1])

          #      else: #Horizontal Shift

                    #x = self.random.uniform(-heading[1] , heading[0])
           #         y = heading[1]

         #   elif heading[0] < 0 and heading[1] < 0:

          #      if self.random.randrange(0 , 2) < 1: #Vertical Shift

           #         x = heading[0]
            #        y = self.random.uniform(heading[1] , -heading[0])

             #   else: #Horizontal Shift

              #      x = self.random.uniform(heading[0] , -heading[1])
               #     y = heading[1]

            #elif heading[0] < 0 and heading[1] > 0:

             #   if self.random.randrange(0 , 2) < 1: #Vertical Shift

              #      x = heading[0]
               #     y = self.random.uniform(heading[0] , heading[1])

                #else: #Horizontal Shift

                 #   x = self.random.uniform(heading[0] , heading[1])
                  #  y = heading[1]

           # elif heading[0] > 0 and heading[1] < 0:

            #    if self.random.randrange(0 , 2) < 1: #Vertical Shift

             #       x = heading[0]
              #      y = self.random.uniform(heading[1] , heading[0])

               # else: #Horizontal Shift

                #    x = self.random.uniform(heading[1] , heading[0])
                 #   y = heading[1]

            #elif heading == (0,0):

              #  x = 1
             #   y = 1


            #norm = (x ** 2 + y ** 2) ** 0.5
           # xDisplacement = x / (norm * 0.5)
          #  yDisplacement = y / (norm * 0.5)
            
            #self.model.space.move_agent(self , (self.pos[0] + xDisplacement , self.pos[1] + yDisplacement))


   # def isTouching(self , other:Agent):
    #    d = self.model.space.get_distance(self.pos , other.pos)
     #   l = d - self.internalR - other.internalR
      #  if l <= 0:
       #     return True 
        #return False





    

