# This file contains code to extract the inequalities
# from a given STN
# Authors: Viva Ojha and Joon Lee
# Contact Info: vmojha@g.hmc.edu and joolee@g.hmc.edu

import math
import random
from pulp import *
import numpy as np
from stntools.stn import STN

# \fn getInequal(STN)
#  \brief gets the inequality constraints from a Robot Brunch STN; converts
#         it to the dictionary format we use


def getInequal(RBS):
    # get a list of all the stn edges
    edgeList = RBS.getAllEdges()

    # create a dictionary of inequalities
    inequalDict = {}

    # convert edge into an inequality and store in
    # dictionary
    for edge in edgeList:
        start = edge.i
        end = edge.j

        # Add to dictionary
        if start not in inequalDict:
            inequalDict[start] = []
        if end not in inequalDict:
            inequalDict[end] = []

        if edge.Cij != float('inf'):
            inequalDict[start].append([start, end, edge.Cij])
            inequalDict[end].append([start, end, edge.Cij])
        if edge.Cji != float('inf'):
            inequalDict[start].append([end, start, edge.Cji])
            inequalDict[end].append([end, start, edge.Cji])

    return inequalDict

# calculates the chebyshev center and radius; based on the code from
# Polybots '17 but works with our STN representation


def chebyshev(inequalDict):
    events = [LpVariable("t" + str(k)) for k in sorted(inequalDict.keys())]
    n = len(events)
    prob = LpProblem("LP", LpMaximize)

    # maximize the radius r
    r = LpVariable("r", lowBound=0)
    prob += r

    # add constraints to our lp
    for l in inequalDict.values():
        for constraint in l:
            # print constraint
            i = constraint[0]
            j = constraint[1]
            coeff1 = 0
            coeff2 = 0

            if i > 0:
                coeff1 = 1
            if j > 0:
                coeff2 = 1

            norm = math.sqrt(coeff1 + coeff2)
            prob += events[j] - events[i] + r * norm <= constraint[2]

    status = prob.solve()
    center = [value(x) for x in events]
    radius = value(r)
    return center, radius

# chebyshev function that works with robot brunch STNs


def chebyshevRBS(RBS):
    inequalDict = getInequal(RBS)
    return chebyshev(inequalDict)

# Takes in an STN and a point, and perturbs one dimension of the point
# by one unit, repeating the process until the point


def perturb(inequalDict, p, vectWalk, posTime, ray):
    if vectWalk:
        return vectorPerturb(inequalDict, p, posTime)
    if ray:
        return rayPerturb(inequalDict, p, posTime)
    p = [coord - p[0] for coord in p[1:]]
    random.seed()
    count = 0
    p = [0] + p
    while True:
        # randomly choose dimension to perturb and which direction
        dim = random.randint(1, len(p)-1)
        coeff = random.randint(0, 1)

        unit = 1

        if coeff == 0 and (not posTime):
            unit *= -1

        p[dim] += unit

        l = inequalDict[dim]

        # check that all constraints involving the perturbed dimension
        # still hold
        for ineq in l:
            start = ineq[0]
            end = ineq[1]
            constraint = ineq[2]

            if p[end] - p[start] > constraint:
                return count

        count += 1

# Takes in an STN in Robot Brunch's format and uses the perturb function
def perturbRBS(RBS, p):
    inequalDict = getInequal(RBS)
    return perturb(inequalDict, p)

# Function that perturbs a given schedule, but in a normalized vector
# that only goes in positive direction (supposed to account for time
# by assuming that everything happens later and can't happen earlier
# than scheduled)
def vectorPerturb(inequalDict, p, posTime):
    count = 0
    p = [coord - p[0] for coord in p[1:]]
    p = [0] + p
    while True:
        # Get a random vector that only goes in the positive directions
        # for all dimensions
        vector = [0] * len(p)
        for i in range(1, len(p)):
            if posTime:
                vector[i] = abs(np.random.normal(0, 0.1, 1)[0])
            else:
                vector[i] = np.random.normal(0, 0.1, 1)[0]
        # Calculate magnitude of our random vector
        magnitude = 0.0
        for i in range(1, len(p)):
            magnitude += vector[i] ** 2
        magnitude = magnitude ** 0.5

        # Normalize our vector to have a total magnitude of 1
        for i in range(1, len(p)):
            vector[i] /= (magnitude * 1.0)
            p[i] += vector[i]

        if not inSTN(inequalDict, p):
            return count

        count += 1

# This function chooses a random direction and determines how far the
# furthest boundary in that direction is
def rayPerturb(inequalDict, p, posTime):
    p = [coord - p[0] for coord in p[1:]]
    p = [0] + p
    count = 0
    # Get a random vector 
    vector = [0] * len(p)
    for i in range(1, len(p)):
        if posTime:
            vector[i] = abs(np.random.normal(0, 0.1, 1)[0])
        else:
            vector[i] = np.random.normal(0, 0.1, 1)[0]
    # Calculate magnitude of our random vector
    magnitude = 0.0
    for i in range(1, len(p)):
        magnitude += vector[i] ** 2
    magnitude = magnitude ** 0.5
    magnitude *= 1.0
    vector = [c/magnitude for c in vector]
    mult = 1
    while True:
        if count > 100: mult*=2
        for i in range(1, len(p)):
            p[i] += mult*vector[i]
        if not inSTN(inequalDict, p):
            return count
        count += mult*1

def distance(p1, p2):
    sum = 0.0
    p1 = [coord - p1[0] for coord in p1[1:]]
    p2 = [coord - p2[0] for coord in p2[1:]]
    for i in range(len(p1)):
        sum += (p1[i] - p2[i]) ** 2
    dist = sum ** 0.5
    return dist

# Determines whether a schedule is in an STN's solution space


def inSTN(S, p):
    for i in range(len(p)):
        for ineq in S[i]:
            start = ineq[0]
            end = ineq[1]
            constr = ineq[2]
            if p[end] - p[start] > constr:
                return False
    return True


cubestn = STN()
cubestn.addCreatedVertex(cubestn.Z_TIMEPOINT)
for i in range(1, 3):
    cubestn.addVertex(i, i, 0)
for i in range(1, 3):
    cubestn.addEdge(0, i, 0, 10)

trianglestn = STN()
trianglestn.addCreatedVertex(trianglestn.Z_TIMEPOINT)
for i in range(1, 3):
    trianglestn.addVertex(i, i, 0)
for i in range(1, 3):
    trianglestn.addEdge(0, i, 0, 10)
trianglestn.addEdge(1, 2, 0, 10)

