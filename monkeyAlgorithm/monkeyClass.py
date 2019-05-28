# -*- coding: utf-8 -*-
"""
Created on Mon Jun 12 13:36:14 2017

@author: rafael.torrese
"""
import numpy as np
import pandas as pd
from pprint import pprint
from itertools import product,cycle,combinations
from collections import Counter
import random
from copy import deepcopy

from containerClass import Container
from circleClass import Circle



class Monkey(Container):
    'Class for generates Monkey object'
    def __init__(self,radii=[2.4,5.3,1.5],length=50, width=50, numPointsLength=20, numPointsWidth=20,p=0.001):
        super().__init__(length,width,numPointsLength,numPointsWidth)
        self.radii = radii # list
        self.probability = p

# ==============================================================================
    @property
    def initialPosition(self):
        '''
        Generate DataFrame of shape (len(radii)*len(points), M). This data frame only generates monkeys with points inside the container, points in the frontier are excluded.

        radii:list
        --------------------
           List of Radii Objects

        points:list
        --------------------
        List of tuples Points inside container.
        [(x1,y1),(x2,y2),(x3,y3),...,(xi,yi),...,(xn,yn)]

        M:int
        --------------------
            Number of monkeys
        positions:list
            List of

        p:float
        --------------------
            Probability for get ones
        '''

        positionsIndex = self.positionsLabels
        lengthVector = len(self.radii)*len(self.points())
        vec = np.random.binomial(1,self.probability,size=lengthVector)
        return pd.DataFrame(vec,index=positionsIndex)

# ==============================================================================
    @property
    def positionsLabels(self):
        'This functions gets positions as index for a dataFrame pandas object'
        return [(x,y,r)  for x,y in self.points() for r in self.radii ] # list of positions, tuple(k,x,y)

# ==============================================================================
    def getSolutionsDF(self):
        '''
        Get list of tuples that represent solutions from a pandas data frame. Solutions are represented by a monkey
        return:
            List of tuples(xCoord,yCoord,radius). This list of tuples is a postion for monkey
        '''
        df = self.initialPosition
        # if every component in the solutions equals to one then save in a list
        return df[df[0] == 1].index.tolist()

# ==============================================================================
    @property
    def objectiveValue(self):
        'Get values of Objective function'
        tuplesolution = self.getSolutionsDF()
        listOfTuples = list(zip(*tuplesolution)) # yield -> [(x,x,x),(y,y,y),(r,r,r)]
        objSolution = listOfTuples[2] # only radii
        nObjects = Counter(objSolution)
        return sum((np.pi)*r*r*n for r,n in nObjects.items())/(self.areaContainer)
 # ==============================================================================
    @property
    def listObjects(self):
        'Creates a list of objects from a list of tuple solution [(x,x,x,...,x,x),(y,y,y,...,y,y),(r,r,r,...,r,r)]'
        list_TupleSolution = self.getSolutionsDF()
        return [Circle(*t) for t in list_TupleSolution]




# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
class Monkey2:
    'Class for Monkey Algorithm'
    def __init__(self,grid=None,radii=None,objects=None):
        self.grid = grid # Container Class
        if radii == None:
            self.radii = [] # Radii List
        else:
            self.radii = list(radii)
        if objects == None:
            self.objects = [] # List of objects
        else:
            self.objects = list(objects) # List of objects

    @property
    def initialGrid(self):
        'Initial grid for selecting random initial points'
        coords_xInside = self.grid.xcoords[1:-1] # only coordinates inside container, for this reason remove the first and last element
        coords_yInside = self.grid.ycoords[1:-1] # only coordinates inside container
        Radius = min(self.radii) # choose a minimum radii
        x,y = coords_xInside[0],coords_yInside[0] # choose the very first point (corner point)
        circle = Circle(x,y,Radius) # set a circle in a choosen point
        while circle.is_outside(self.grid): # test if circle is outside
            coords_yInside = coords_yInside[1:-1] # update list
            coords_xInside = coords_xInside[1:-1] # update list
            x,y = coords_xInside[0],coords_yInside[0] # choose the first point, this is a corner point
            circle.center = x,y  # set a circle in the new point
        return list(product(coords_xInside,coords_yInside))



    def initialPosition(self,low=2,up=12):
        "Returns List of Objects. This list is a position, i.e, a pattern of circles inside a container"
        initGrid = self.initialGrid # initial grid, only inside points
        #objects = []
        numObj = random.randint(low,up)
        while (len(self.objects) < numObj) and initGrid:
            p = random.choice(initGrid) # select a point from grid, point is a tuple (x,y)
            idxp = initGrid.index(p) # index of selected point in the list initGrid
            r = random.choice(self.radii) # Choose randomly an object
            circle = Circle(*p,radius=r) # Set circle in selected point p
            #if len(self.objects): # if there are at least two cirlcles in the list
                #flag = True
            while any(circle.overlapping(obj) for obj in self.objects):
                #flag = any(circle.overlapping(obj) for obj in self.objects)
                #if flag:
                p = random.choice(initGrid) # select point
                circle.center = p # assign another center
                idxp = initGrid.index(p) # index of selected point
                initGrid.pop(idxp) # remove from init list
            self.removeOutside()
            #if circle.is_outside(self.grid): # test if circle is outside
            #    circle.radius = min(self.radii) # try with minimum radius for verify. It is a point in the inside frontier
            self.objects.append(circle)
            #pts.append(initGrid.pop(idxp)) # remove from init list and put in list of points
            initGrid = set(initGrid) - set(self.grid.coverage(circle) ) # forbidden points for next iteration
            initGrid = list(initGrid) # update grid  points
        for obj,c in zip(self.objects,cycle(Circle.colors)):
            obj.color = c
        return self.objects


    @property
    def position(self):
        'List of objects that represents a position'
        return self.objects


    @position.setter
    def position(self,objs):
        'Set List of objects that represents a position'
        self.objects = objs

    @property
    def objectiveValue(self):
        'Gets Objective Value'
        return sum(obj.area for obj in self.objects)

    @property
    def positionTuple(self):
        'Converts Solution Objects to Solution Tuple (x,y,r)'
        return [(obj.x,obj.y,obj.radius) for obj in self.objects]

    #@property
    def showPosition(self,nombre,flag=False):
        'Draw objects '
        self.grid.drawObjects2(self.position,name=nombre,saveHere=flag)

    def __add__(self,other):
        'Sum two Monkeys'
        sumObjects = self.objects + other.objects
        setObjects = {(obj.x,obj.y,obj.radius) for obj in sumObjects} # remove duplicates
        sumObjects = (Circle(x,y,r,color=c) for (x,y,r),c in zip(setObjects,cycle(Circle.colors)))
        return Monkey2(self.grid,self.radii,sumObjects)


#==============================================================================
#     def removeOverlappings(self):
#         'Detect overlapings from object list'
#         #numberOverlapingsList = []
#         for i,circle_1 in enumerate(self.objects):
#             for circle_2 in self.objects[ i+1: ]:
#                 if circle_2.overlapping(circle_1):
#                     if circle_1.radius <= circle_2.radius:
#                         self.objects.remove(circle_1)
#                     else:
#                         self.objects.remove(circle_2)
#                     #numberOverlapingsList.extend([circle_1,circle_2])
#         #print("There was(were) {} overlapping(s)".format(len(numberOverlapingsList)//2))
#         #numberObjectsOverlapped = list(set(numberOverlapingsList))
#         #return numberOverlapingsList, numberObjectsOverlapped # Number overlappings is different from number OBJECTS overlapped
#         #return numberOverlapingsList
#==============================================================================
    def removeOutside(self):
        'If object is outside of the grid, then remove it'
        self.radii.sort() # sort in place in asecending order
        for c in self.objects:
            while c.is_outside(self.grid): # While object is outside
                idxRadius = self.radii.index(c.radius) # Position of radius in radii list
                c.radius = self.radii[idxRadius - 1] # Select the next small radius


    def removeOverlappings(self,save='n'):
        'If object is overlapping with another object, then remove it'
        toRemove = []
        objects = deepcopy(self.objects)
        for circle1,circle2 in combinations(objects, 2):
            if circle1.overlapping(circle2): # If overlaps
                if circle1.radius < circle2.radius:
#                    objects.remove(circle1)
                    toRemove.append(circle1)
                else:
#                    self.objects.remove(circle2)
                    toRemove.append(circle2) # save circle 2 in list toRemove
        if toRemove: # if exists overlapping objects, otherwise do nothing
            for c in toRemove:
               if c in objects:
                   objects.remove(c)
            self.objects = objects
        if save != 'n': # If you wish to save overlapping list
            return toRemove


    def combiner(self,other,p=0.6):
        'For Watch Jump Process'
        for c in self.objects:
            keep = np.random.binomial(1,p) # to keep an alement
            if c in other.objects and not keep:
                other.objects.remove(c)
            elif c not in other.objects and keep:
                other.objects.append(c)
        return Monkey2(self.grid,self.radii,other.objects)


    def cooperates(self,other,p=0.6):
        'Combines Monkey and Best Monkey with probability p. Probability p is a measure of accordance with best monkey'
        M = deepcopy(self.objects)
        Mbest = deepcopy(other.objects)
        Y = []
        for m in M:
            if m in Mbest:
                Y.append(m)
                Mbest.remove(m)
                M.remove(m)
        # Apply probability of change or reject
        for m in M:
            rejectChange = np.random.binomial(1,1-p) # If elemen reject best monkey element
            if rejectChange: # This means that original element no change
                Y.append(m)
        for m in Mbest: # if Monkey Accept element of best monkey
            acceptChange = np.random.binomial(1,p)
            if acceptChange:
                Y.append(m)
        return Monkey2(self.grid,self.radii, objects= Y)


    def somersault(self,other,a=-1,b=1,threshold=0.3):
        'Monkey Samersault in along the direction pointing to the pivot. Other is Pivot'
        M = deepcopy(self.objects)
        Mp = deepcopy(other.objects) # monkey pivot
        Y = []
        for m in Mp:
            if m in M:
                Y.append(m)
                M.remove(m)
        # Apply probability of change or reject
        for m in M:
            x = random.uniform(a,b)
            p = (x-a)/(b-a) # cumulative uniform distribution
            if p > threshold: # This means that original element no change
                Y.append(m)
        for m in Mp: # if Monkey Accept element of best monkey
            x = random.uniform(a,b)
            p = (x-a)/(b-a) # cumulative uniform distribution
            if p < threshold:
                Y.append(m)
        return Monkey2(self.grid,self.radii, objects= Y)


    def outsideExists(self):
        'Returns True if an object is outside of the grid'
        for c in self.objects:
            if c.is_outside(self.grid):
                return True


    def overlappingExists(self):
        'Check in overlapping exists for at least one object'
        for circle1,circle2 in combinations(self.objects, 2):
            if circle1.overlapping(circle2):
                return True

    @property
    def occupiedPoints(self):
        'Returns Points that have been occupied by objects'
        occupied = []
        for c in self.objects:
            occupied += c.pointsCoveraged(self.grid)
        occupied = set(occupied)
        return list(occupied)

    @property
    def freePoints(self):
        'Return Points that have not been yet occupied'
        whole = set(self.initialGrid)
        occupied = set(self.occupiedPoints)
        free = whole - occupied
        return list(free)


    @property
    def intensityOverlap(self):
        'Measure of the intensity of overlappings (scalar)'
        measure = 0
        for circle1,circle2 in combinations(self.objects, 2):
            overlap = circle1.overlapping(circle2)
            if overlap:
                measure += overlap
        return measure

    @property
    def intensityOut(self):
        'Measure of the intensity of domain violations (scalar)'
        measure = 0
        for obj in self.objects:
            out = obj.is_outside(self.grid)
            if out:
                measure += out # accumulator
        return measure


    def value(self,v='overall'):
        'Value of Monkey'
        v = v.lower()
        if v == 'overall': # All
            return self.objectiveValue - self.intensityOut - self.intensityOverlap
        elif v == 'objective': # Only objective Value
            return self.objectiveValue
        elif v == 'overlap': # Only overlapping
            return self.intensityOverlap
        elif v == 'out': # Only domain violations
            return self.intensityOut
        else:
            print("Don't recognize this option. Value is setting to zero")
            return 0


    def accommodate(self,trials=3):
        'Remove overlappings and place objects in free points while available trials'
        objs = self.removeOverlappings('saveList') # list of objects to accommodate [C1,C2,...,Cn]
        while objs and trials:# Maybe it is imposible to accomodate all objcets, so, we need to specify a number o trials
            p = self.freePoints # save free Points list
            objs.sort(key=lambda c: c.radius) # sort by radius from smaller to bigger
            for obj in objs:
                if p: # if exists points in list
                    point = random.choice(p) # choice random point from free points
                    obj.center = point # Assign new center to object
                    p.remove(point)
            self.objects += objs # update Monkey position list (objects list)
            objs = self.removeOverlappings('saveList') # list of objects to accomadte [C1,C2,...,Cn]
            trials -= 1 # decrease number of trials and test condition

#    def accommodate(self,trials=3):
#            'Remove overlappings and place objects in free points while available trials'
#            objs = self.removeOverlappings('saveList') # list of objects to accommodate [C1,C2,...,Cn]
#            while trials:# Maybe it is imposible to accomodate all objcets, so, we need to specify a number o trials
#                if objs:
#                    p = deepcopy(self.freePoints) # save free Points list
#                    objs.sort(key=lambda c: c.radius) # sort by radius from smaller to bigger
#                    for obj in objs:
#                        if p: # if exists points in list
#                            point = random.choice(p) # choice random point from free points
#                            obj.center = point # Assign new center to object
#                            p.remove(point)
#                    self.objects += objs # update Monkey position list (objects list)
#                    objs = self.removeOverlappings('saveList') # list of objects to accomadte [C1,C2,...,Cn]
#                    trials -= 1 # decrease number of trials and test condition
#                else:
#                    break

#==============================================================================
#     def accommodate(self,trials=3):
#         'Remove overlappings and place objects in free points while available trials'
#         for _ in np.arange(trials):
#             objs = self.removeOverlappings('saveList') # list of objects to accommodate [C1,C2,...,Cn]
#             if objs:
#                 p = self.freePoints # save free Points list
#                 #objs.sort(key=lambda c: c.radius) # sort by radius from smaller to bigger
#                 for obj in objs:
#                     if p: # if exists points in list
#                         point = random.choice(p) # choice random point from free points
#                         obj.center = point # Assign new center to object
#                         p.remove(point)
#                 self.objects += objs # update Monkey position list (objects list)
#             else:
#                 break
#==============================================================================

    def countObjects(self):
        'Count different elements in objects list '
        self.objects.sort(key=lambda s: s.radius)
        idxList = []
        for obj1,obj2 in zip(self.objects, self.objects[1:]):
            if obj1.radius != obj2.radius:
                idxList.append(self.objects.index(obj2))
        counts = {r:0 for r in self.radii} # initilize counter
        for c in self.objects:
            for k in counts:
                if c.radius == k:
                  counts[k] += 1 # count
        return counts # returns a Dict

    @property
    def totalAreaObjects(self):
        return np.pi*sum(r*r*v for r,v in self.countObjects().items())

    @property
    def density(self):
        return self.totalAreaObjects / self.grid.areaContainer


class Population:
    def __init__(self,individuals:list):
        self.individuals = individuals # List of objects

    def best(self):
        'Returns the best individual in Population'
        bestIndex, bestIndividual = max(enumerate(self.individuals), key=lambda M: M[1].value())
        return bestIndividual


    def showIndividual(self,number=0):
        #numberIndividuals = len(self.individuals)
        self.individuals[number].showPosition




# ==============================================================================
# Functions
# ==============================================================================
def detectOverlapings(objects:list):
    'Detect overlapings from object list'
    #numberOverlapingsList = []
    for i,circle_1 in enumerate(objects):
        for circle_2 in objects[ i+1: ]:
            if circle_1.overlapping(circle_2):
                if circle_1.radius <= circle_2.radius:
                    print("Removing object {}".format(circle_1))
                    objects.remove(circle_1)
                else:
                    print("Removing object {}".format(circle_2))
                    objects.remove(circle_2)
                #numberOverlapingsList.extend([circle_1,circle_2])
    #print("There was(were) {} overlapping(s)".format(len(numberOverlapingsList)//2))
    #numberObjectsOverlapped = list(set(numberOverlapingsList))
    #return numberOverlapingsList, numberObjectsOverlapped # Number overlappings is different from number OBJECTS overlapped
    #return numberOverlapingsList

# ==============================================================================
def detectObjectsOutside(objects:list,container):
    'Detect objects that lies outside of container from object list'
    objectsOutsideList = []
    for obj in objects:
        if obj.is_outside(container):
            objectsOutsideList.append(obj)
    return objectsOutsideList

def getPositions(radii,points):
    '''
    This functions gets positions as index for a dataFrame object pandas
    ----------------------------------------
    radii:
       List of radii objects for packing problem
    ------------------------------------------------------------
    points:list
    ------------------------------------------------------------
       List of tuples Points inside container.points:
    ------------------------------------------------------------
    return list of tuples (k,x,y)
    ------------------------------------------------------------
    '''
    return [(x,y,r)  for x,y in points for r in radii] # list of positions, tuple(x,y,r)
# ==============================================================================


# ==============================================================================

if __name__ == "__main__":
    grid = Container(60,50,20,20)
    #radiiList = [3.5,7.1,4,5.6,1.5,2.3,8.9]
    radiiList = [5.0,3.8,  2.9,  5.6,   2.2,  4.6,6.3,3.2]
    random.seed(845)
#    monkey2 = Monkey2(grid,radiiList)
#    monkey2.position = monkey2.initialPosition()


    numonk = 1
    # Initial Population
    X0 = [Monkey2(grid,radiiList) for _ in range(numonk)]
    delta_1 = [Monkey2(grid,radiiList) for _ in range(numonk)]
    delta_2 = [Monkey2(grid,radiiList) for _ in range(numonk)]



    # Initialize positions
#==============================================================================
#     for m in range(numonk):
#         X0[m].position = X0[m].initialPosition()
#         delta_1[m].position = delta_1[m].initialPosition()
#         delta_2[m].position = delta_2[m].initialPosition()
#==============================================================================

    # Initialize positions
    for x0,d1,d2 in zip(X0,delta_1,delta_2):
        x0.position = x0.initialPosition(2,6)
        d1.position = d1.initialPosition(2,6)
        d2.position = d2.initialPosition(2,8)


#==============================================================================
    # ==================== CLIMB PROCESS ====================
#==============================================================================
#==============================================================================
#     X1 = [X0[m] + delta_1[m] for  m in range(numonk)] # positions
#     X2 = [X0[m] + delta_2[m] for  m in range(numonk)]
#==============================================================================

    X1 = [x0 + d1 for x0,d1 in zip(X0,delta_1)] # positions
    X2 = [x0 + d2 for x0,d2 in zip(X0,delta_2)]




    for m in range(numonk):
        X0[m].showPosition
        delta_1[m].showPosition
        delta_2[m].showPosition
        X1[m].showPosition
        X2[m].showPosition
        #X1[m].removeOverlappings()
        #X2[m].removeOverlappings()
        #X1[m].showPosition




    for m in range(numonk):
        if X1[m].objectiveValue > X2[m].objectiveValue and X1[m].objectiveValue > X0[m].objectiveValue:
            X0[m] = X1[m]
            print("X1")
        elif X2[m].objectiveValue > X1[m].objectiveValue and X2[m].objectiveValue > X0[m].objectiveValue:
            X0[m] = X2[m]
            print("X2")
        else:
            print("X0")

    for m in range(numonk):
        X1[m].showPosition
        X2[m].showPosition
        X0[m].showPosition

    for m in X0:
        #m.showPosition
        m.removeOverlappings()
        #m.showPosition
    #X0[0].showPosition
#==============================================================================
#     for m in range(numonk):
#         print(X0[m].objectiveValue,X1[m].objectiveValue,X2[m].objectiveValue)
#         print(len(X0[m].objects),len(X1[m].objects),len(X2[m].objects))
#==============================================================================
#==============================================================================
    # ==================== WATCH JUMP PROCESS ====================
#==============================================================================
    y = [Monkey2(grid,radiiList) for _ in range(numonk)]
    # Initialize positions
    for m in range(numonk):
        y[m].position = y[m].initialPosition(3,5)


#    idx = 1
#    X0[idx].showPosition
#    print(len(X0[idx].objects))
#    y[idx].showPosition
#    print(len(y[idx].objects))
    #Y = [X0[m] + y[m] for m in range(numonk)]
    Y = [y[m].combiner(X0[m]) for m in range(numonk)] # combine X0 with WatchJump
    for m in range(numonk):
        if Y[m].objectiveValue > X0[m].objectiveValue:
            X0[m] = Y[m]


    #X0[1].showPosition
    monkeyPopulation = Population(X0)
    bestMonkey = monkeyPopulation.best()


    #bestMonkey.showPosition
    #X0[0].showPosition
    #monkeyPopulation.showIndividual(number=1)
    newMonkey = X0[0].cooperates(bestMonkey,p=0.7)
    #newMonkey.showPosition
#    Y[idx].showPosition
#    print(len(Y[idx].objects))
    #print(X0[1].objects[0] in delta_1[1].objects)


    #X0[1].showPosition
    pivotMonkey = random.choice(X0) # select random monkey

    #X0[0].showPosition
    #pivotMonkey.showPosition
    monkeySomersault = X0[0].somersault(pivotMonkey,threshold=0.8)
    #monkeySomersault.showPosition
    print(len(monkeySomersault.position))
#    monkeySomersault.removeOverlappings()
    monkeySomersault.accommodate()
    #monkeySomersault.showPosition
    print(len(monkeySomersault.position))
    monkeySomersault.removeOutside()
    #monkeySomersault.showPosition
    monkeySomersault.countObjects()
    for obj in monkeySomersault.position:
        print(obj)


#    grid.markPoints(monkeySomersault.freePoints)
#==============================================================================
#     for m in range(numonk):
#         X0[m].showPosition
#         y[m].showPosition
#         Y[m].showPosition
#==============================================================================


#==============================================================================
#     for m in monkeys:
#         print("{:.1%}".format(m.objectiveValue(m.initialPosition())/grid.areaContainer))
#==============================================================================





