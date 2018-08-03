import random
import randomWalk
import randomSchedule
import plotting
import os
import json
import multiprocessing as mp
import test
import stnConverter
import copy
import matplotlib.pyplot as plt
import numpy as np

# Function to get a list of all the constraints in an STN
def getList(stn):
    lists = stn.values()
    final_list = lists[0]

    for list in lists[1:]:
        
        for ineq in list:
            
            if not ineq in final_list:
                final_list.append(ineq)
    
    return final_list

# Function to check whether a point satisfies all the constraints
# of an STN
def inLiSTN(stn, p):
    
    for ineq in stn:
        start = ineq[0]
        end = ineq[1]
        constr = ineq[2]
        
        if p[end] - p[start] > constr:
            return False
    
    return True

# Implementation of a random "shave". Chooses random constraints,
# one at a time, and tightens that constraint slightly. The goal
# is to see how many times we can do this before the point p is no
# longer a valid schedule 
def randomShave(stn, p):
    stnCopy = copy.deepcopy(stn)
    ineqList = getList(stnCopy)
    count = 0
    shaved = False
    shavedCount = 0

    while True:
        rand = random.randint(0,len(ineqList)-1)

        # We subtract the constraint value by one whether or not the
        # constraint value is negative because if the constraint is
        # positive, we are lowering the upper bound, and if the value 
        # is negative, we are increasing the lower bound. Both serve to
        # tighten the constraints between events. 
        backward = []
        for qeni in ineqList:
            if qeni[1] == ineqList[rand][0] and qeni[0] == ineqList[rand][1]:
                backward = qeni
    
        if (ineqList[rand][2] * -1) < backward[2]:
            ineqList[rand][2] -= 1

        else: 
            shavedCount += 1

        count += 1
        #print count
        if shavedCount >= len(ineqList):
            shaved = True
        
        if not inLiSTN(ineqList, p) or shaved:
            return count

def melt(stn, p):
    stnCopy = copy.deepcopy(stn)

    ineqList = realMinList(stnCopy)
    count = 0
    melted = False
    while True:
        meltedCount = 0
        for ineq in ineqList:
            backward = []

            for qeni in ineqList:
                if qeni[1] == ineq[0] and qeni[0] == ineq[1]:
                    backward = qeni
    
            if (ineq[2] * -1) < backward[2]:
                ineq[2] -= 1
            else: 
                meltedCount += 1
        count += 1

        if meltedCount == len(ineqList):
            melted = True
        
        if (not inLiSTN(ineqList, p)) or melted:
            return count

def advanceTime(stn):
    ineqList = stn
    for ineq in ineqList:
        if ineq[1] == 0:
            ineq[2] -= 1
    return ineqList

def realMinList(stn):
    origList = getList(stn)
    newList = getList(stnConverter.matrixToDict(test.minimal(stnConverter.dictToMatrix(stn))))
    finalList = []
    for ineq in origList:
        if ineq in newList:
            finalList.append(ineq)
    return finalList

# This function is for multiprocessing - it takes in a tuple with an stn, the number of desired shaves
# and the number of random points to consider and then returns a list with six elements:
# the average number of shaves for cheb, cheb's melt, cent's shaves, cent's melt, randos' average shaves
# rando's melt
def mpRandomShaver(tup):
    genstn = tup[0]
    numShaves = tup[1]
    numRandos = tup[2]

    cent = randomSchedule.centroid(genstn)
    cheb = randomWalk.chebyshev(genstn)[0]
    mid = randomSchedule.midPoint(genstn)
    
    centAvg = 0.0
    chebAvg = 0.0
    randoAvg = 0.0
    midAvg = 0.0

    for i in range(numRandos):
        rando = randomSchedule.simpleRandSched(genstn)
        
        for j in range(numShaves):
            randoAvg += randomShave(genstn, rando)
    
    randoAvg /= (numRandos*numShaves)

    for j in range(numShaves):
        centAvg += randomShave(genstn, cent)
        chebAvg += randomShave(genstn, cheb)
        midAvg += randomShave(genstn, mid)
    
    centAvg /= float(numShaves)
    chebAvg /= float(numShaves)
    midAvg /= float(numShaves)
    print "SHAAAVVVVEEEE " + str(tup[3])
    return (chebAvg, centAvg, midAvg, randoAvg)

# Really doesn't belong here but it needs the get list function...
# Also, needs more debugging (dist_to_bound part is off); doesn't currently generate schedules within STN
def uniformRandSched(S, numPoints = 100):
    pointList = []

    # start at arbitrary point, we take cheb
    curr = randomWalk.chebyshev(S)[0]
    # print curr
    curr = [0] + [coord - curr[0] for coord in curr[1:]]
    # print curr
    # retrieve a list of the constraints
    constraints = getList(S)
    # go through the number of points to cycle through
    i = 0
    while i < numPoints + 100:
        # set up a vector in a random direction
        vector = [0 for n in range(len(S.keys()))]

        for j in range(1, len(vector)):
            vector[j] = np.random.normal(0, 0.1, 1)[0]
        # print vector
        # calculate magnitude of this vector]
        magnitude = 0.0
        for j in range(1, len(vector)):
            magnitude += vector[j] ** 2
        magnitude = magnitude ** 0.5
        # print magnitude
        for j in range(1, len(vector)):
            vector[j] /= magnitude
        # print vector

        # initialize thetas for determining the lower and upper bounds for traveling along vector
        theta_pos = 10000
        theta_neg = -10000
        # loop through the bounds and calculate distance to each, updating thetas as we go, to get bounds
        for bound in constraints:
            const = bound[2]
            start = bound[0]
            end = bound[1]
            dist_to_bound = (const - (curr[end] - curr[start]))/(vector[end] - vector[start])
            if dist_to_bound > 0: theta_pos = min(theta_pos, dist_to_bound)
            if dist_to_bound < 0: theta_neg = max(theta_neg, dist_to_bound)
        # print theta_pos
        # print theta_neg
        # print vector
        # initialize next point
        # new_point = [0 for n in range(len(curr))]
        new_point = copy.deepcopy(curr)
        # uniformly choose displacement along vector and calculate the appropriate point
        displacement = np.random.uniform(theta_neg, theta_pos)
        for j in range(1, len(vector)):
            new_point[j] += vector[j] * displacement
        
        curr = new_point
        if i > 99:
            pointList.append(curr)
        i += 1
        # print 'i = ' + str(i)

    return pointList

# cube = randomWalk.getInequal(randomWalk.cubestn)
# xList, yList = uniformRandSched(cube)
# triangle = randomWalk.getInequal(randomWalk.trianglestn)
# xList, yList = uniformRandSched(triangle)
# print rando
# print randomWalk.inSTN(cube, rando)
# plotting.justPlot(xList, 'x', yList, 'y', 'Generating Random Points by Point/Vector')


# with open(os.path.abspath('./Data/' + 'ALL_OF_THE_DATA_30000_WITH_ICE_CREAM' + '.json'), 'r') as fp:
#     datastuff = json.load(fp)

# print len(datastuff)
# print len(datastuff[24])

# plotting.justPlot(datastuff[2], 'sph met', datastuff[25], 'rando shave', 'asdas')

# randomSchedule.getMarg(datastuff, 22, 25)

# stns = datastuff[0][:5]
# tups = []

# for i in range(len(stns)):
#     tups.append((randomSchedule.fixStupidJSON(stns[i]), 100, 100, i))

# procPool = mp.Pool(64)

# results = procPool.map(mpRandomShaver, tups)

# procPool.close()

# datastuff.append([])
# datastuff.append([])
# datastuff.append([])
# datastuff.append([])

# for i in range(len(results)):
#     datastuff[22].append(results[i][0])
#     datastuff[23].append(results[i][1])
#     datastuff[24].append(results[i][2])
#     datastuff[25].append(results[i][3])

# with open(os.path.abspath('./Data/' + 'ALL_OF_THE_DATA_30000_WITH_ICE_CREAM' + '.json'), 'w') as fp:
#     json.dump(datastuff, fp)

# chebShaves = [result[0] for result in results]
# chebMelt = [result[1] for result in results]
# centShaves = [result[2] for result in results]
# centMelt = [result[3] for result in results]
# randoShaves = [result[4] for result in results]
# randoMelt = [result[5] for result in results]

# stuff = [chebShaves, chebMelt, centShaves,centMelt,randoShaves,randoMelt]

# with open(os.path.abspath('./Data/' + 'ICE_CREAM_10000' + '.json'), 'w') as fp:
#     json.dump(stuff, fp)
