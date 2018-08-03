# This file contains functions that make use of random schedules
# for different STNs
# Authors: Viva Ojha and Joon Lee
# Contact Info: vmojha@g.hmc.edu and joolee@g.hmc.edu

import randomWalk
import stnConverter
import test
import random
import json
import metricFunctions
import os
import sys
import multiprocessing as mp
import plotting
import matplotlib.pyplot as plt
import numpy as np
import copy
import math

numWalks = 100

# Function to get a list of all the constraints in an STN
def getList2(stn):
    lists = stn.values()
    final_list = lists[0]

    for list in lists[1:]:
        
        for ineq in list:
            
            if not ineq in final_list:
                final_list.append(ineq)
    
    return final_list

# Function that tells us how far away the centroid and
# the chebyshev center are
def distMetric(S, dim):
    one = randomWalk.chebyshev(S)[0]
    cent = centroid(S)
    return randomWalk.distance(one, cent)


# Convert our STN form to matrix form (required for use
# of code in Polybots repository that we didn't write)
# then call getRandomSchedule function from there
def simpleRandSched(dict):
    matrix_form = stnConverter.dictToMatrix(dict)
    minimal_matrix = test.minimal(matrix_form)
    # print minimal_matrix
    return test.getRandomSchedule(minimal_matrix)



# Same as above function but conversion is from robot brunch
# stn format to matrix format
def simpleRandSchedRBS(S):
    matrix_form = stnConverter.listToMatrix(S)
    minimal_matrix = test.minimal(matrix_form)
    return test.getRandomSchedule(minimal_matrix)

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
    constraints = getList2(S)
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

# Computes an approximation of the centroid by taking an
# avg of numPoints number of random schedules
def centroid(ineqDict, numPoints):
    n = len(ineqDict.keys())
    centroid = [0] * n

    for i in range(numPoints):
        # p = simpleRandSched(ineqDict)
        p = uniformRandSched(numPoints)
        
        for j in range(n):
            centroid[j] += p[j]

    for i in range(n):
        centroid[i] /= float(numPoints)
        # centroid[i] += 0.5

    return centroid


# First convert robot brunch stn to our inequality
# format, then call the above centroid function
def centroidRBS(S):
    ineqDict = randomWalk.getInequal(S)
    return centroid(ineqDict)

def flexPoint(ineqDict, metric):
    minimalSTN = test.minimal(stnConverter.dictToMatrix(ineqDict))
    indices = range(0, len(minimalSTN))
    flexSchedule = [0]*(len(minimalSTN))

    while indices:
        maxFlex = -1
        maxIndex = -1
        time = -1
        
        for index in indices:
            timepoint = ((minimalSTN[0][index] + minimalSTN[index][0])/2) - minimalSTN[index][0]
            newMat = copy.deepcopy(minimalSTN)
            newMat[index][0] = -timepoint
            newMat[0][index] = timepoint
            if len(indices) > 1:
                currFlex = metricFunctions.sphericalMetric(stnConverter.matrixToDict(newMat), len(indices)-1, metric)
            else:
                currFlex = 1
            if currFlex > maxFlex:
                # print 'currFlex: '  + str(currFlex)
                # print 'maxFlex: ' + str(maxFlex)
                maxFlex = currFlex
                maxIndex = index
                time = timepoint
        
        indexToAssign = maxIndex
        # delete the chosen index from the list of indices still to be assigned
        # print "THIS INDEX"
        # print indexToAssign
        indices.remove(indexToAssign)
        # assign the index a random time that can satisfy the STN
        flexSchedule[indexToAssign] = time

        # recalculate the minimal STN with this event nailed down
        minimalSTN[indexToAssign][0] = -time
        minimalSTN[0][indexToAssign] = time
        minimalSTN = test.minimal(minimalSTN)
    return flexSchedule

def rigidPoint(ineqDict, metric):
    minimalSTN = test.minimal(stnConverter.dictToMatrix(ineqDict))
    indices = range(0, len(minimalSTN))
    flexSchedule = [0]*(len(minimalSTN))

    while indices:
        maxFlex = 10000000
        maxIndex = -1
        time = -1
        
        for index in indices:
            timepoint = ((minimalSTN[0][index] + minimalSTN[index][0])/2) - minimalSTN[index][0]
            newMat = copy.deepcopy(minimalSTN)
            newMat[index][0] = -timepoint
            newMat[0][index] = timepoint
            if len(indices) > 1:
                currFlex = metricFunctions.sphericalMetric(stnConverter.matrixToDict(newMat), len(indices)-1, metric)
            else:
                currFlex = 1
            if currFlex < maxFlex:
                # print 'currFlex: '  + str(currFlex)
                # print 'maxFlex: ' + str(maxFlex)
                maxFlex = currFlex
                maxIndex = index
                time = timepoint
        
        indexToAssign = maxIndex
        # delete the chosen index from the list of indices still to be assigned
        # print "THIS INDEX"
        # print indexToAssign
        indices.remove(indexToAssign)
        # assign the index a random time that can satisfy the STN
        flexSchedule[indexToAssign] = time

        # recalculate the minimal STN with this event nailed down
        minimalSTN[indexToAssign][0] = -time
        minimalSTN[0][indexToAssign] = time
        minimalSTN = test.minimal(minimalSTN)
    return flexSchedule


def midPoint(ineqDict):
    minimalSTN = test.minimal(stnConverter.dictToMatrix(ineqDict))
    indices = range(0, len(minimalSTN))
    midSchedule = [0]*(len(minimalSTN))

    while indices:
        
        indexToAssign = random.choice(indices)
        time = ((minimalSTN[0][indexToAssign] + minimalSTN[indexToAssign][0])/2) - minimalSTN[indexToAssign][0]
        # delete the chosen index from the list of indices still to be assigned
        # print "THIS INDEX"
        # print indexToAssign
        indices.remove(indexToAssign)
        # assign the index a random time that can satisfy the STN
        midSchedule[indexToAssign] = time

        # recalculate the minimal STN with this event nailed down
        minimalSTN[indexToAssign][0] = -time
        minimalSTN[0][indexToAssign] = time
        minimalSTN = test.minimal(minimalSTN)

    return midSchedule

def avgMidPoint(ineqDict):
    dim = len(ineqDict.keys())
    p = [0] * dim
    numPoints = 20
    for i in range(numPoints):
        mid = midPoint(ineqDict)
        
        for j in range(dim):
            p[j] += mid[j]
    
    for i in range(dim):
        p[i] /= numPoints
    
    return p

# def jimPoint(ineqDict):
#     minimalSTN = test.minimal(stnConverter.dictToMatrix(ineqDict))
#     indices = range(0, len(minimalSTN))
#     jimSchedule = [0]*(len(minimalSTN))
#     sortedList = sorted(minimalSTN[0])

#     for i in range(len(sortedList)):
#         time = sortedList[i]/2.0
        
#         indexToAssign = minimalSTN[0].index(sortedList[i])
#         print indexToAssign
        
#         jimSchedule[indexToAssign] = time
#         minimalSTN[indexToAssign][0] = -time
#         minimalSTN[0][indexToAssign] = time
#         minimalSTN = test.minimal(minimalSTN)
#     return jimSchedule

def jimPoint(ineqDict):
    minimalSTN = test.minimal(stnConverter.dictToMatrix(ineqDict))
    indices = range(0, len(minimalSTN))
    jimSchedule = [0]*(len(minimalSTN))

    while indices:
        minFlex = 100000000000
        minIndex = -1
        time = -1
        
        for index in indices:
            timepoint = ((minimalSTN[0][index] + minimalSTN[index][0])/2) - minimalSTN[index][0]
            newMat = copy.deepcopy(minimalSTN)
            newMat[index][0] = -timepoint
            newMat[0][index] = timepoint
            if len(indices) > 1:
                currFlex = (minimalSTN[0][index] + minimalSTN[index][0])
            else:
                currFlex = 10000000
            if currFlex < minFlex:
                # print 'currFlex: '  + str(currFlex)
                # print 'maxFlex: ' + str(maxFlex)
                minFlex = currFlex
                minIndex = index
                time = timepoint
        
        indexToAssign = minIndex
        # delete the chosen index from the list of indices still to be assigned
        # print "THIS INDEX"
        # print indexToAssign
        indices.remove(indexToAssign)
        # assign the index a random time that can satisfy the STN
        jimSchedule[indexToAssign] = time

        # recalculate the minimal STN with this event nailed down
        minimalSTN[indexToAssign][0] = -time
        minimalSTN[0][indexToAssign] = time
        minimalSTN = test.minimal(minimalSTN)
    return jimSchedule

# Computes the earliest possible start schedule by going
# through the inequalities and setting each point
# to its earliest possible value
def earlyBird(ineqDict):
    n = len(ineqDict.keys())
    earlyPoint = [0] * n

    matrix = stnConverter.dictToMatrix(ineqDict)
    matrix = test.minimal(matrix)
    ineqDict = stnConverter.matrixToDict(matrix)

    for i in range(n-1, 0, -1):
        earlyPoint[i] = -50000
        for ineq in ineqDict[i]:
            start = ineq[0]
            end = ineq[1]
            constraint = ineq[2]

            if i == start and (start < end or end == 0):
                earlyPoint[i] = max(
                    earlyPoint[end] - constraint, earlyPoint[i])
    return earlyPoint


# Given an STN and a point in that STN, calculate the radius of the
# inscribed sphere centered around that point
def closestBoundary(S, p):
    S = stnConverter.matrixToDict(test.minimal(stnConverter.dictToMatrix(S)))

    # throw out the zeroeth coordinate of the point because that is
    # the zero timepoint and isn't part of a given schedule
    p = [0] + [coord - p[0] for coord in p[1:]]

    r = sys.maxint

    # Loop through all the constraints and adjust r, the radius,
    # effectively finding the distance to the closest boundary of
    # the solution space for the STN
    for l in S.values():
        for inequal in l:
            start = inequal[0]
            end = inequal[1]
            constraint = inequal[2]

            temp = constraint - p[end] + p[start]
            if start != 0 and end != 0:
                temp /= math.sqrt(2)

            r = min(temp, r)
    return r


# Given an STN and a point in that STN, calculate the distance to
# the furthest boundary
def furthestBoundary(S, p):
    # throw out the zeroeth coordinate of the point because that is
    # the zero timepoint and isn't part of a given schedule
    # S = stnConverter.matrixToDict(test.minimal(stnConverter.dictToMatrix(S)))

    p = [0] + [coord - p[0] for coord in p[1:]]

    r = -sys.maxint

    # Loop through all the constraints and adjust r, the radius,
    # effectively finding the distance to the closest boundary of
    # the solution space for the STN
    for l in S.values():
        for inequal in l:
            start = inequal[0]
            end = inequal[1]
            constraint = inequal[2]

            temp = constraint - p[end] + p[start]
            if start != 0 and end != 0:
                temp /= math.sqrt(2)
            r = max(temp, r)
    return r

def furthestBoundaryFix(S, p):
    # throw out the zeroeth coordinate of the point because that is
    # the zero timepoint and isn't part of a given schedule
    # S = stnConverter.matrixToDict(test.minimal(stnConverter.dictToMatrix(S)))

    p = [0] + [coord - p[0] for coord in p[1:]]

    r = -sys.maxint


    # Loop through all the constraints and adjust r, the radius,
    # effectively finding the distance to the closest boundary of
    # the solution space for the STN
    for l in S.values():
        for inequal in l:
            start = inequal[0]
            end = inequal[1]
            constraint = inequal[2]

            temp = constraint - p[end] + p[start]
            if start != 0 and end != 0:
                temp /= math.sqrt(2)
            # fixing frivolous boundaries code!! only designed for axial ones but including others
            # isn't hard
            # obtains point of intersection for checking and then checks if it's in stn
            check_point = copy.deepcopy(p)
            checkers = []

            use = True

            if start == 0:
                check_point[end] -= temp
                checkers = S[end]

            if end == 0:
                check_point[start] += temp
                checkers = S[start]

            for ineq in checkers:
                if check_point[ineq[1]] - check_point[ineq[0]] > ineq[2]:
                    use = False

            if use:
                r = max(temp, r)
    return r

# Given stupid JSON, returns non-stupid STN
# JSON is stupid because it makes all dictionary keys into strings
# instead of ints
def fixStupidJSON(stupidJSON):
    actualData = {}
    for key in stupidJSON.keys():
        actualData[int(key)] = []
        for ineq in stupidJSON[key]:
            actualData[int(key)].append([ineq[0],ineq[1],ineq[2]])
    return actualData

# MP Center Helper Function Random Walker
def mpCenterRandomWalker(tup):
    genstn = tup[0]
    cheb = tup[1]
    cent = tup[2]
    birdie = tup[3]
    posTime = tup[4]
    ray = tup[5]
    return (randomWalk.perturb(genstn, cheb, not ray, posTime, ray), 
            randomWalk.perturb(genstn, cent, not ray, posTime, ray),
            randomWalk.perturb(genstn, birdie, not ray, posTime, ray))


# MP Rando Helper Function Random Walker
def mpRandoRandomWalker(tup):
    genstn = tup[0]
    numWalks = tup[1]
    posTime = tup[2]
    ray = tup[3]
    curr = simpleRandSched(genstn)
    randAvg = 0
    for i in range(numWalks):
        randAvg += randomWalk.perturb(genstn, curr, True, posTime, ray)
    return randAvg


# Function for generating data and computing possible useful quantities that 
# we may want to examine for a given dataset of STNs
def pointBrawl(numTrials, dim, manWalk=False, posTime=False, randDim=False, multi=False, ray=False):

    chebMargin = 0
    centMargin = 0
    birdMargin = 0
    chebVCentMargin = 0
    chebVBirdMargin = 0
    centVBirdMargin = 0

    # Matric to store all the different quantities we may want to compute in a
    # simulation. The numbers on the right refer to indices of the matrix
    lol = [[],  # sphereMetric          0
           [],  # burgerMetric          1
           [],  # wilsonMetric          2
           [],  # radRatio              3
           [],  # naiveMetric           4
           [],  # chebMargin            5
           [],  # centMargin            6
           [],  # birdMargin            7
           [],  # chebVCentMargin       8
           [],  # chebVBirdMargin       9
           [],  # centVBirdMargin       10
           [],  # randomWalkCheb        11
           [],  # randomWalkCent        12
           [],  # randomWalkBird        13
           [],  # randomWalkManyPoints  14
           [],  # dimension             15
           [],  # STNs.                 16
           ]

    # Run the point arena on numTrials number of STNs
    procPool = mp.Pool(2)
    for i in range(numTrials):
        print "\nTrial " + str(i)
        chebAvg = 0
        centAvg = 0
        birdAvg = 0
        randAvg = 0

        numWalks = 100

        # Randomly choose the dimension between 2 and 10 events
        if randDim:
            dim = random.randint(2, 10)

        # generate STN with dim events, each constrained by 0 to 50
        genstn = stnConverter.matrixToDict(test.generateRandomSTN(dim,0,50))

        # minimize the STN because why not
        genstn = stnConverter.matrixToDict(
                test.minimal(stnConverter.dictToMatrix(genstn)))

        lol[15].append(dim)
        lol[16].append(genstn)

        # Metrics for analyzing an STN depending on what we are looking for.
        # This could be how "spherelike" the solution space is,
        sphereMet = metricFunctions.sphericalMetric(genstn, dim, "spherical")
        burgerMet = metricFunctions.sphericalMetric(genstn, dim, "burgers")
        wilsonMet = metricFunctions.sphericalMetric(genstn, dim, "wilson")
        radRatio = metricFunctions.sphericalMetric(genstn, dim, "radRatio")
        naiveMet = metricFunctions.sphericalMetric(genstn, dim, "hell")

        # Compute the possible "centers", or points we think may be optimal
        cheb = randomWalk.chebyshev(genstn)[0]
        cent = centroid(genstn)
        birdie = earlyBird(genstn)

        # Append the computed metric values to their respective rows
        lol[0].append(sphereMet)
        lol[1].append(burgerMet)
        lol[2].append(wilsonMet)
        lol[3].append(radRatio)
        lol[4].append(naiveMet)

        # Take the average of 100 randomWalks using 100 different points
        if not multi:
            for i in range(100):
                curr = simpleRandSched(genstn)
                for i in range(numWalks):
                    randAvg += randomWalk.perturb(genstn, curr, not ray, posTime, ray)
        else:
            tuplist = [(genstn,numWalks,posTime,ray)]*100
            results = procPool.map(mpRandoRandomWalker,tuplist)
            randAvg = sum(results)

        randAvg /= 100*float(numWalks)

        # Take the average of 100 randomWalks starting at each "center"
        if not multi:
            for i in range(numWalks):
                chebAvg += randomWalk.perturb(genstn, cheb, not ray, posTime, ray)
                centAvg += randomWalk.perturb(genstn, cent, not ray, posTime, ray)
                birdAvg += randomWalk.perturb(genstn, birdie, not ray, posTime, ray)
        else:
            tup = [(genstn, cheb, cent, birdie, posTime, ray)] * numWalks
            results = procPool.map(mpCenterRandomWalker,tup)
            chebAvg = sum([res[0] for res in results])
            centAvg = sum([res[1] for res in results])
            birdAvg = sum([res[2] for res in results])
            
        chebAvg /= float(numWalks)
        centAvg /= float(numWalks)
        birdAvg /= float(numWalks)

        print "cheb " + "Average: " + str(chebAvg)
        print "cent " + "Average: " + str(centAvg)
        print "bird " + "Average: " + str(birdAvg)
        print "rand " + "Average: " + str(randAvg)

        # Compute the different margins of the centers vs each other
        if chebAvg == 0 and centAvg == 0:
            chebVCentMargin = 0
        else:
            chebVCentMargin = 2 * (chebAvg - centAvg) / (chebAvg + centAvg)
        if chebAvg == 0 and birdAvg == 0:
            chebVBirdMargin = 0
        else:
            chebVBirdMargin = 2 * (chebAvg - birdAvg) / (chebAvg + birdAvg)
        if centAvg == 0 and centAvg == 0:
            centVBirdMargin = 0
        else:
            centVBirdMargin = 2 * (centAvg - birdAvg) / (centAvg + birdAvg)

        # Compute the different margins of the centers vs random points
        if chebAvg == 0 and randAvg == 0:
            chebMargin = 0
        else:
            chebMargin = 2 * (chebAvg - randAvg) / (chebAvg + randAvg)
        if centAvg == 0 and randAvg == 0:
            centMargin = 0
        else:
            centMargin = 2 * (centAvg - randAvg) / (centAvg + randAvg)
        if birdAvg == 0 and randAvg == 0:
            birdMargin = 0
        else:
            birdMargin = 2 * (birdAvg - randAvg) / (birdAvg + randAvg)

        # Append the computed margins to their respective rows
        lol[5].append(chebMargin)
        lol[6].append(centMargin)
        lol[7].append(birdMargin)
        lol[8].append(chebVCentMargin)
        lol[9].append(chebVBirdMargin)
        lol[10].append(centVBirdMargin)

        # Append the computed randomWalk lengths starting at each center to
        # their respective rows
        lol[11].append(chebAvg)
        lol[12].append(centAvg)
        lol[13].append(birdAvg)
        lol[14].append(randAvg)

        print "SPHEERE " + str(sphereMet)
        print "BURGERS " + str(burgerMet)
        print "WILSOON " + str(wilsonMet)
        print "RAAAADD " + str(radRatio)
        print "NAIIIVE " + str(naiveMet)

    title = "THE_GRAND_BATTLE_ON_"
    if randDim:
        title += "MANY"
    else:
        title += str(dim)
    title += "_DIMENSIONS_"
    title += "AND_" + str(numTrials) + "_MILLION_TRIALS"
    if posTime:
        title += "_AS_TIME_MOVES_FORWARD!"
    if ray:
        title += "...IN A RAY!"

    # Code to store data in a json file
    with open(os.path.abspath('./Data/' + title + '.json'), 'w') as fp:
        json.dump(lol, fp)

    return lol

# This function takes in an inequality dictionary and prints out the
# format we need to input into polymake to visualize our polytope

def printIneq(ineqDict):

    polytope = '$p = new Polytope(INEQUALITIES=>['
    n = len(ineqDict.keys())
    for l in ineqDict.values():
        for ineq in l:
            start = ineq[0]
            end = ineq[1]
            constraint = ineq[2]
            ineqArr = [0] * n

            ineqArr[0] = constraint

            for i in range(n):
                if i != 0:
                    if i == start:
                        ineqArr[i] = 1
                    if i == end:
                        ineqArr[i] = -1

            polytope += str(ineqArr) + ', '

    polytope = polytope[0:-2]
    polytope += ']);'
    print polytope

def mpCentRandomWalker(tup):
    genstn = tup[0]
    numWalks = tup[1]
    posTime = tup[2]
    ray = tup[3]
    cent = centroid(genstn)
    centAvg = 0
    for i in range(numWalks):
        centAvg += randomWalk.perturb(genstn, cent, True, posTime, ray)
    centAvg /= float(numWalks)
    print "Trial " + str(tup[4])
    return centAvg

def mpMidWalker(tup):
    genstn = tup[0]
    numWalks = tup[1]
    posTime = tup[2]
    ray = tup[3]
    mid = flexPoint(genstn, "spherical")
    midAvg = 0
    for i in range(numWalks):
        midAvg += randomWalk.perturb(genstn, mid, True, posTime, ray)
    midAvg /= float(numWalks)
    print "Trial " + str(tup[4])
    return midAvg

def mpRigWalker(tup):
    genstn = tup[0]
    numWalks = tup[1]
    posTime = tup[2]
    ray = tup[3]
    mid = rigidPoint(genstn, "spherical")
    midAvg = 0
    for i in range(numWalks):
        midAvg += randomWalk.perturb(genstn, mid, True, posTime, ray)
    midAvg /= float(numWalks)
    print "Trial " + str(tup[4])
    return midAvg

def calcMargin(avg1, avg2):
    # return 2.0 * (avg1 - avg2)/(avg1 + avg2 +0.000000001)
    # if avg1 < 1 or avg2 < 1:
    #     # print avg1
    #     # print avg2
    #     return 0.0
        # print avg1 / avg2
    if avg1 < 0.1 or avg2 < 0.1:
        return 0.0
    if avg2 == 0:
        return 0.0
    # if (avg1 / avg2) > 100.0:
    #     return 0.0 
    return avg1 / avg2

def listMargin(list1, list2):
    retList = []
    for i in range(len(list1)):
        retList.append(calcMargin(list1[i],list2[i]))
    return retList

def mpPointBrawl(tup):
    genstn = tup[0]
    dim = len(genstn.keys()) - 1

    sphereMet = metricFunctions.sphericalMetric(genstn, dim, "spherical")
    burgerMet = metricFunctions.sphericalMetric(genstn, dim, "burgers")/(dim*(dim-1)/2)
    wilsonMet = metricFunctions.sphericalMetric(genstn, dim, "wilson")/dim
    radRatio = metricFunctions.sphericalMetric(genstn, dim, "radRatio")
    naiveMet = metricFunctions.sphericalMetric(genstn, dim, "hell")
    
    cheb = randomWalk.chebyshev(genstn)[0]
    cent = centroid(genstn)
    birdie = earlyBird(genstn)
    mid = midPoint(genstn)
    
    randAvg = 0.0
    randAvgRay = 0.0
    randAvgRayPos = 0.0
    #randAvgShaves = 0.0

    for i in range(100):
        curr = simpleRandSched(genstn)

        for j in range(numWalks):
            randAvg += randomWalk.perturb(genstn, curr, 1, 0, 0)
            randAvgRay += randomWalk.perturb(genstn, curr, 0, 0, 1)
            randAvgRayPos += randomWalk.perturb(genstn, curr, 0, 1, 1)
            #andAvgShaves += randomMelt.randomShave(genstn, curr)

    randAvg /= i*j
    randAvgRay /= i*j
    randAvgRayPos /= i*j
    #randAvgShaves /= i*j

    chebAvg = 0.0
    chebAvgRay = 0.0
    chebAvgRayPos = 0.0
    #chebAvgShaves = 0.0

    centAvg = 0.0
    centAvgRay = 0.0
    centAvgRayPos = 0.0
    #centAvgShaves = 0.0

    birdAvg = 0.0
    birdAvgRay = 0.0
    birdAvgRayPos = 0.0
    #birdAvgShaves = 0.0

    midAvg = 0.0
    midAvgRay = 0.0
    midAvgRayPos = 0.0
    #midAvgShaves = 0.0

    for i in range(numWalks):
        chebAvg += randomWalk.perturb(genstn, cheb, 1, 0, 0)
        chebAvgRay += randomWalk.perturb(genstn, cheb, 0, 0, 1)
        chebAvgRayPos += randomWalk.perturb(genstn, cheb, 0, 1, 1)
        #chebAvgShaves += randomMelt.randomShave(genstn, cheb)

        centAvg += randomWalk.perturb(genstn, cent, 1, 0, 0)
        centAvgRay += randomWalk.perturb(genstn, cent, 0, 0, 1)
        centAvgRayPos += randomWalk.perturb(genstn, cent, 0, 1, 1)
        #centAvgShaves += randomMelt.randomShave(genstn, cent)

        birdAvg += randomWalk.perturb(genstn, birdie, 1, 0, 0)
        birdAvgRay += randomWalk.perturb(genstn, birdie, 0, 0, 1)
        birdAvgRayPos += randomWalk.perturb(genstn, birdie, 0, 1, 1)
        #birdAvgShaves += randomMelt.randomShave(genstn, bird)

        midAvg += randomWalk.perturb(genstn, mid, 1, 0, 0)
        midAvgRay += randomWalk.perturb(genstn, mid, 0, 0, 1)
        midAvgRayPos += randomWalk.perturb(genstn, mid, 0, 1, 1)
        #midAvgShaves += randomMelt.randomShave(genstn, mid)

    chebAvg /= numWalks
    chebAvgRay /= numWalks
    chebAvgRayPos /= numWalks
    #chebAvgShaves /= numWalks

    centAvg /= numWalks
    centAvgRay /= numWalks
    centAvgRayPos /= numWalks
    #centAvgShaves /= numWalks

    birdAvg /= numWalks
    birdAvgRay /= numWalks
    birdAvgRayPos /= numWalks
    #birdAvgShaves /= numWalks

    midAvg /= numWalks
    midAvgRay /= numWalks
    midAvgRayPos /= numWalks
    #midAvgShaves /= numWalks

    print "TRIAL!! " + str(tup[1])

    return [genstn, dim, sphereMet, burgerMet, wilsonMet, radRatio, naiveMet, 
            chebAvg, centAvg, birdAvg, midAvg, randAvg,
            chebAvgRay, centAvgRay, birdAvgRay, midAvgRay, randAvgRay,
            chebAvgRayPos, centAvgRayPos, birdAvgRayPos, midAvgRayPos, randAvgRayPos]
            #chebAvgShaves, centAvgShaves, midAvgShaves, randAvgShaves]

def getMarg(datastuff, num1, num2):
    largMarg = 0.00
    bigWalk = 0.0
    sphereMet = 5
    burgerMet = 10
    wilson = 10
    naive = 500
    dim = 0

    n = 0
    negn = 0
    neutraln = 0
    numstn = 0

    margins = listMargin(datastuff[num1], datastuff[num2])

    for i in range(len(datastuff[0])):
        if datastuff[1][i] > dim and (datastuff[num1][i] > bigWalk and datastuff[num2][i] > bigWalk): 
            numstn += 1.0
            if margins[i] > 1.0:
                n += 1
            if margins[i] < 1.0:
                negn += 1
            # if chebVecMargins[i] == largMarg:
            #     neutraln += 1
    # print numstn
    # print n
    print n / (n + negn + 0.0), negn / (n + negn + 0.0), numstn, (n + negn)

# cube = randomWalk.getInequal(randomWalk.cubestn)
# triangle = randomWalk.getInequal(randomWalk.trianglestn)
# xList = []
# yList = []

# for i in range(100):
#     sched = simpleRandSched(triangle)
#     xList.append(sched[1])
#     yList.append(sched[2])

# plotting.justPlot(xList, 'x', yList, 'y', 'Generating Random Points by Dimension')
# cube = randomWalk.getInequal(randomWalk.cubestn)
# cheb = randomWalk.chebyshev(cube)[0]
# rad = randomWalk.chebyshev(cube)[1]
# print rad
# cheb = [p - cheb[0] for p in cheb]
# print randomWalk.inSTN(cube, cheb)
# cent = centroid(cube)
# chebAv = 0
# centAv = 0
# for i in range(100):
#     chebAv += randomWalk.perturb(cube, cheb, 1, 0, 0)
#     centAv += randomWalk.perturb(cube, cent, 1, 0, 0)
# chebAv /= 100
# centAv /= 100
# print cheb
# print chebAv

# print cent
# print centAv

# with open(os.path.abspath('./Data/' + 'ALL_OF_THE_DATA_30000' + '.json'), 'r') as fp:
#     datastuff = json.load(fp)

# metVals = []

# for i in range(30):
#     stn = datastuff[0][i]
#     stn = fixStupidJSON(stn)
#     cheb = randomWalk.chebyshev(stn)[0]
#     metVal = schedMetrics.wilsonCubeMet(stn, cheb)
#     metVals.append(metVal)

# plotting.justPlot(metVals, "metricVals", datastuff[7][:30], "random walks from cheb", "title")

# print len(datastuff[21])
# first = randomWalk.getInequal(randomWalk.cubestn)
# first = datastuff[0][34]
# first = fixStupidJSON(first)

# cent = centroid(first)
# avgMid = avgMidPoint(first)

# print cent
# print avgMid

# dist = randomWalk.distance(cent, avgMid)
# rad = closestBoundary(first, cent)

# print dist
# print rad 
# print furthestBoundary(first, cent)
# print  dist / rad 
# print 'cheb'
# getMarg(7, 11)

# getMarg(12, 16)

# getMarg(17, 21)

# print 'cent'
# getMarg(8, 11)

# getMarg(13, 16)

# getMarg(18, 21)

# print 'mid'
# getMarg(10, 11)

# getMarg(15, 16)

# getMarg(20, 21)

# margins = listMargin(datastuff[7], datastuff[8])
# plotting.justPlot(datastuff[2], 'Spherical Metric Value', margins, 'Chebyshev vs. Centroid Margin', 'Chebyshev/Centroid Margin v. Spherical Metric Value')

# print 'chebVcent'
# getMarg(7, 8)

# getMarg(12, 13)

# getMarg(17, 18)

# print 'chebVmid'
# getMarg(7, 10)

# getMarg(12, 15)

# getMarg(17, 20)

# print 'centVmid'
# getMarg(8, 10)

# getMarg(13, 15)

# getMarg(18, 20)

# marg = listMargin(datastuff[21], datastuff[19])

# plotting.justPlot(datastuff[2], 'Spherical Metric Value', marg, 'Random Walk Length', 'Random Walk from Random Points v. Spherical Metric Value')

# # n = 100
# getMarg(19,21)
# centMin = 1.0
# centMax = 0.0
# centAv = 0.0

# midMin = 1.0
# midMax = 0.0
# midAv = 0.0

# boundLengths = []
# midDists = []
# centDists = []
# sphereMets = datastuff[2][:n]
# for stn in datastuff[0][:n]:
#     cheby = randomWalk.chebyshev(stn)
#     stn = fixStupidJSON(stn)
#     cent = centroid(stn)
#     cheb = cheby[0]
#     mid = midPoint(stn)
#     far = furthestBoundary(stn, cheb)
#     boundLengths.append(far)

#     chebCent = randomWalk.distance(cheb, cent)
#     chebCent = chebCent/far
#     centDists.append(chebCent)

#     chebMid = randomWalk.distance(cheb, mid)
#     chebMid = chebMid/far
#     midDists.append(chebMid)

#     centMax = max(chebCent, centMax)
#     centMin = min(chebCent, centMin)
#     midMax = max(chebMid, midMax)
#     midMin = min(chebMid, midMin)
#     midAv += chebMid
#     centAv += chebCent
#     print "FURTHEST CHEB BOUNDARY " + str(far)
#     print "CHEB CENT DIST " + str(chebCent)
#     print "CHEB MID DIST " + str(chebMid)

# print "centMin " + str(centMin)
# print "centMax " + str(centMax)
# print "centAv " + str(centAv/n)
# plotting.justPlot(centDists, 'centDists', sphereMets, 'sphereMets', 'title')

# print "midMin " + str(midMin)
# print "midMax " + str(midMax)
# print "midAv " + str(midAv/n)
# plotting.justPlot(midDists, 'midDists', sphereMets, 'sphereMets', 'title')


# plotting.justPlot(listMargin(datastuff[10], datastuff[11]), "cent", listMargin(datastuff[7], datastuff[11]), "mid", "FIGHT")

# chebBeatsCent = sorted([(fixStupidJSON(datastuff[0][i]),
#                          calcMargin(datastuff[7][i], datastuff[8][i]))
#                         for i in range(30000)],
#                         key=lambda x: x[1]
# )

# onestn = chebBeatsCent[5][0]

# chebAvg = 0
# randAvg = 0

# for i in range(100):
#     print "HI " + str(i)
#     rando = simpleRandSched(loser)
#     for j in range(100):
#         chebAvg += randomWalk.perturb(loser, cheb, 1, 0, 0)
#         randAvg += randomWalk.perturb(loser, rando, 1, 0, 0)
# chebAvg /= 10000.0
# randAvg /= 10000.0
# print chebAvg
# print randAvg

# onestn = fixStupidJSON(datastuff[0][5723])
# cheb = randomWalk.chebyshev(onestn)[0]
# cent = centroid(onestn)
# mid = midPoint(onestn)
# chebAvg = 0
# centAvg = 0
# midAvg = 0

# for i in range(500):
#     print "HI " + str(i)
#     chebAvg += randomWalk.perturb(onestn, cheb, 1, 0, 0)
#     centAvg += randomWalk.perturb(onestn, cent, 1, 0, 0)   
#     midAvg += randomWalk.perturb(onestn, mid, 1, 0, 0)   

# print chebAvg / 500.0
# print centAvg / 500.0
# print midAvg / 500.0
# print randomWalk.distance(cheb, cent)
# print randomWalk.distance(cent, mid)
# print randomWalk.distance(cheb, mid)
# print cheb
# print cent
# print mid

# print "chebvrandwins: "
# wins = 0
# for marg in listMargin(datastuff[12], datastuff[16]):
#     if marg > 0: wins += 1
# print wins/ 30000.0

# print "centvrandwins: "
# wins = 0
# for marg in listMargin(datastuff[13], datastuff[16]):
#     if marg > 0: wins += 1
# print wins/ 30000.0

# print "birdvrandwins: "
# wins = 0
# for marg in listMargin(datastuff[14], datastuff[16]):
#     if marg > 0: wins += 1
# print wins/ 30000.0

# print "midvrandwins: "
# wins = 0
# for marg in listMargin(datastuff[15], datastuff[16]):
#     if marg > 0: wins += 1
# print wins/ 30000.0

# with open(os.path.abspath('./Data/' + 'THE_GRAND_BATTLE_ON_MANY_DIMENSIONS_AND_10000_MILLION_TRIALS...IN A RAY!_FOR_REAL' + '.json'), 'r') as fp:
#     datastuffRay = json.load(fp)

# with open(os.path.abspath('./Data/' + 'THE_GRAND_BATTLE_ON_MANY_DIMENSIONS_AND_10000_MILLION_TRIALS_AS_TIME_MOVES_FORWARD!...IN A RAY!_FOR_REAL' + '.json'), 'r') as fp:
#     datastuffRayPos = json.load(fp)

# procPool = mp.Pool(64)
# stns = datastuff[16]
# stnsRay = datastuffRay[16]
# stnsRayPos = datastuffRayPos[16]
# superSTNs = stns + stnsRay + stnsRayPos

# superDuperSTNs = [stnConverter.matrixToDict(
#                 test.minimal(stnConverter.dictToMatrix(fixStupidJSON(stn)))) for stn in superSTNs]

# tups = []

# for i in range(len(superDuperSTNs)):
#     tups.append((superDuperSTNs[i], i))

# results = procPool.map(mpPointBrawl, tups)

# lol = [[0 for i in range(len(results))] for j in range(len(results[0]))]
# for i in range(len(lol)):
#     for j in range(len(lol[0])):
#         lol[i][j] = results[j][i]

# with open(os.path.abspath('./Data/ALL_OF_THE_DATA_10.json'), 'w') as fp:
#     json.dump(lol, fp)

# with open(os.path.abspath('./Data/ALL_OF_THE_DATA_10.json'), 'r') as fp:
#     datastuff = json.load(fp)
# print datastuff[2]
# print datastuff[7]
# plotting.justPlot(datastuff[7], 'spherical', datastuff[8], 'randomWalk Cheb', 'title')

# numTrials = 100
# vecTups = []
# rayTups = []
# rayPosTups = []

# for i in range(numTrials):
#     stn = stns[i]
#     stn = fixStupidJSON(stn)
#     stnRay = stnsRay[i]
#     stnRay = fixStupidJSON(stnRay)
#     stnRayPos = stnsRayPos[i]
#     stnRayPos = fixStupidJSON(stnRayPos)
#     vecTups.append((stn, 100, 0, 0, i))
#     rayTups.append((stnRay, 100, 0, 1, i))
#     rayPosTups.append((stnRayPos, 100, 1, 1, i))

# vecResult = procPool.map(mpMidWalker, vecTups)

# for i in range(numTrials):
#     stn = stns[i]
#     stn = fixStupidJSON(stn)
#     stnRay = stnsRay[i]
#     stnRay = fixStupidJSON(stnRay)
#     stnRayPos = stnsRayPos[i]
#     stnRayPos = fixStupidJSON(stnRayPos)
#     vecTups.append((stn, 100, 0, 0, i))
#     rayTups.append((stnRay, 100, 0, 1, i))
#     rayPosTups.append((stnRayPos, 100, 1, 1, i))

# rigResult = procPool.map(mpRigWalker, vecTups)

# rayResult = procPool.map(mpMidWalker, rayTups)
# rayPosResult = procPool.map(mpMidWalker, rayPosTups)
    
# with open(os.path.abspath('./Data/midPointVectorWalks.json'), 'w') as fp:
#     json.dump(vecResult, fp)

# with open(os.path.abspath('./Data/midPointRayWalks.json'), 'w') as fp:
#     json.dump(rayResult, fp)

# # with open(os.path.abspath('./Data/midPointRayWalksPosTime.json'), 'w') as fp:
# #     json.dump(rayPosResult, fp)


# with open(os.path.abspath('./Data/' + 'midPointVectorWalks' + '.json'), 'r') as fp:
#     vecWalks = json.load(fp)
    
# # with open(os.path.abspath('./Data/' + 'midPointRayWalks' + '.json'), 'r') as fp:
# #     rayWalks = json.load(fp)

# # with open(os.path.abspath('./Data/' + 'midPointRayWalksPosTime' + '.json'), 'r') as fp:
# #     rayPosWalks = json.load(fp)

# margins = listMargin(vecResult[:50], vecWalks)

# plotting.justPlot(datastuff[0][:50], 'sphereMet', margins, 'margins', 'TITLEEE')

# margins = listMargin(vecResult, datastuff[11][:numTrials])

# plotting.justPlot(datastuff[0][:numTrials], 'sphereMet', margins, 'margins', 'TITLEEE')

# margins = listMargin(vecResult, datastuff[12][:numTrials])

# plotting.justPlot(datastuff[0][:numTrials], 'sphereMet', margins, 'margins', 'TITLEEE')

# margins = listMargin(vecResult, datastuff[14][:numTrials])

# plotting.justPlot(datastuff[0][:numTrials], 'sphereMet', margins, 'margins', 'TITLEEE')

# margins = listMargin(vecResult, rigResult)

# plotting.justPlot(datastuff[0][:numTrials], 'sphereMet', margins, 'margins', 'TITLEEE')
