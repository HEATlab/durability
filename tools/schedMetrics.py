import os
import sys
import json
import math
import copy
import itertools
import numpy as np

import randomSchedule
import randomWalk
import stnConverter
import plotting 
import randomMelt
import test
import metricFunctions

# cube = randomWalk.getInequal(randomWalk.cubestn)
# cheb = randomWalk.chebyshev(cube)[0]

# with open(os.path.abspath('./Data/' + 'ALL_OF_THE_DATA_30000' + '.json'), 'r') as fp:
#     datastuff = json.load(fp)



def chebMet(stn, p):
    cheb = randomWalk.chebyshev(stn)[0]
    p = [coord - p[0] for coord in p[1:]]
    p = [0] + p
    return -1*randomWalk.distance(cheb, p)

def closeFarMet(stn, p):
    return randomSchedule.closestBoundary(stn, p)/randomSchedule.furthestBoundary(stn, p)

def secondClosestBoundary(stn, p):
    S = stnConverter.matrixToDict(test.minimal(stnConverter.dictToMatrix(stn)))

    # throw out the zeroeth coordinate of the point because that is
    # the zero timepoint and isn't part of a given schedule
    p = [0] + [coord - p[0] for coord in p[1:]]

    r = sys.maxint - 1
    nex = sys.maxint

    # Loop through all the constraints and adjust r, the radius,
    # effectively finding the distance to the closest boundary of
    # the solution space for the STN
    for l in S.values():
        for inequal in l:
            start = inequal[0]
            end = inequal[1]
            constraint = inequal[2]

            temp = constraint - p[end] + p[start]
            if temp < r:
                nex = r
                r = temp
            elif temp == r:
                continue
            elif temp < nex:
                nex = temp
    return nex


# This metric is the schedule flexibility extension of the sphere metric
# In practice it just returns the distance to the closest boundary
def sphereMet(stn, p):
    return randomSchedule.closestBoundary(stn, p)

# This metric is like sphereMet but finds the distance to the average 
# boundary instead of the minimal one
# Now geometric! And a whole lot better!
def geoMet(stn, p, useDim=True):
    #S = stnConverter.matrixToDict(test.minimal(stnConverter.dictToMatrix(stn)))
    S = stn
    power = 10.0

    if useDim: 
        power = float(len(stn.keys()) - 1)
    # throw out the zeroeth coordinate of the point because that is
    # the zero timepoint and isn't part of a given schedule
    p = [0] + [coord - p[0] for coord in p[1:]]

    r = 1

    constraints = randomMelt.getList(S)

    # print len(constraints)
    # Loop through all the constraints and multiply to r each time
    # the distance to the each boundary
    for inequal in constraints:
        start = inequal[0]
        end = inequal[1]
        constraint = inequal[2]
        temp = constraint - p[end] + p[start]
        if start != 0 and end != 0:
            temp /= math.sqrt(2)
        r *= temp

    if r < 0:
        return 0

    # Take the geometric mean and then raise to the power of 10
    r = math.pow(r, 1.0/len(constraints))  # / (10000000000.0)
    return r


# THEORETICALLY COMPLETE BUT HAVEN'T HAD THE CHANCE TO TEST THIS YET (MAY NOT EVEN COMPILE)
def geoMetFix(stn, p, useDim=True):
    #S = stnConverter.matrixToDict(test.minimal(stnConverter.dictToMatrix(stn)))
    S = stn
    power = 10.0

    if useDim: 
        power = float(len(stn.keys()) - 1)
    # throw out the zeroeth coordinate of the point because that is
    # the zero timepoint and isn't part of a given schedule
    p = [0] + [coord - p[0] for coord in p[1:]]

    r = 1

    constraints = randomMelt.getList(S)

    # keep track of how many constraints are actually used
    unused = 0

    # print len(constraints)
    # Loop through all the constraints and multiply to r each time
    # the distance to the each boundary
    for inequal in constraints:
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
        unused_bool = False

        if start == 0:
            check_point[end] += temp
            checkers = S[end]

        if end == 0:
            check_point[start] -= temp
            checkers = S[start]

        for ineq in checkers:
            if check_point[ineq[1]] - check_point[ineq[0]] > ineq[2]:
                use = False
                unused_bool = True

        if unused_bool: unused += 1



        if use:
            r *= temp


    if r < 0:
        return 0

    # Take the geometric mean and then raise to the power of 10
    r = math.pow(r, 1.0/(len(constraints)-unused)) # / (10000000000.0)
    return r

# This metric is a schedule flexibility metric inspired from the Wilson
# metric 
# It finds the largest inscribed cube centered around the given point
# Ratio to sphereMet should be between 1 and sqrt(n)
def wilsonCubeMet(stn, p):
    # adjust first coordinate if necessary
    p = [coord - p[0] for coord in p[1:]]
    p = [0] + p

    dim = len(stn.keys()) - 1
    
    # create dictionary to keep track of the 2^n cube vectors
    tracker = {}
    keys = itertools.product([1,-1], repeat=dim)

    for vec in keys:
        # normalize each vector to have magnitude 1 and fill
        # in the dictionary to make each vector map to the
        # original point initially
        norm = [0]*len(vec)
        for i in range(len(vec)):
            norm[i] = vec[i]/math.sqrt(len(p) - 1)
        tracker[tuple(norm)] = copy.deepcopy(p)
    
    dist = 0

    # loop through and move the original point by 1 in each
    # cube vector direction
    while True:
        for key in tracker.keys():
            for i in range(len(key)):
                tracker[key][i+1] += key[i]
            if not randomWalk.inSTN(stn, tracker[key]):
                for i in range(len(key)):
                    tracker[key][i+1] -= key[i]
                dist += randomSchedule.closestBoundary(stn, tracker[key])
                return dist
        dist += 1

# This schedule flexibility metric is inspired from the Naive metric
# for STNs
# For each event dimension, this metric looks at the "slack" in either
# direction and adds the smaller one
# Takes an optional argument to take the larger slack each time instead
def naiveMet(stn, p, takeMax = True):
    S = stnConverter.matrixToDict(test.minimal(stnConverter.dictToMatrix(stn)))
    # adjust first coordinate if necessary
    p = [coord - p[0] for coord in p[1:]]
    p = [0] + p

    # get STN's matrix form 
    mat = test.minimal(stnConverter.dictToMatrix(stn))
    
    # initialize metric value to zero
    met = 1
    n = 0
    # loop through each event and add the higher or lower slack as appropriate
    for i in range(1,len(p)):
        if takeMax: 
            r =  max(p[i], mat[0][i] - p[i])
            met *= abs(r)
            n += 1
        else:
            # print "MAT0I"
            # print mat[0][i]
            # print "PI"
            # print p[i]
            r = min(p[i], mat[0][i] - p[i])
            met *= abs(r)
            n += 1
    met = math.pow(met, 1.0 / n)
    return met

# This schedule flexibility metric is inspired from the Naive metric
# for STNs
# For each event dimension, this metric looks at the "slack" in either
# direction and adds the smaller one
# Takes an optional argument to take the larger slack each time instead
def hunsMet(stn, p, takeMax = True):
    S = stnConverter.matrixToDict(test.minimal(stnConverter.dictToMatrix(stn)))
    # adjust first coordinate if necessary
    p = [coord - p[0] for coord in p[1:]]
    p = [0] + p

    # get STN's matrix form 
    mat = test.minimal(stnConverter.dictToMatrix(stn))
    
    # initialize metric value to zero
    met = 0

    # loop through each event and add the higher or lower slack as appropriate
    for i in range(1,len(p)):
        if takeMax: met += max(p[i], mat[0][i] - p[i])
        else:
            # print "MAT0I"
            # print mat[0][i]
            # print "PI"
            # print p[i]
            met += min(p[i], mat[0][i] - p[i])
    for l in S.values():
        for inequal in l:
            start = inequal[0]
            end = inequal[1]
            constraint = inequal[2]

            temp = constraint - p[end] + p[start]
            met += temp
    return met / len(S.keys())


# This function takes in an STN and a number of desired points
# It outputs several lists that provide the durability metric values
# for all of these random points as well as the random walk 
# and random shave values for these points
def randomPointsMetricsWalksShaves(stn, numPoints):
    print 'CALL FUNCTION'
    stn = randomSchedule.fixStupidJSON(stn)
    # initializes lists to store walk/shave/metric values
    walk_vals = []
    const_tight_vals = []
    sph_met_vals = []
    wil_met_vals = []
    nai_met_vals = []
    geo_met_vals = []

    print 'BEFORE FIRST LOOP'
    # generates data for numPoints points
    for i in range(numPoints):

        # generates random point, initializes walk and shave values
        point = randomSchedule.simpleRandSched(stn)
        average_walk = 0.0
        average_shave = 0.0

        print 'BEFORE SECOND LOOP'
        # performs 100 random walks and shaves
        for j in range(100):
            average_shave += randomMelt.randomShave(stn, point)
            average_walk += randomWalk.perturb(stn, point, 1, 0, 0)

        # adds in the point's metric values
        sph_met_vals.append(sphereMet(stn, point))
        wil_met_vals.append(wilsonCubeMet(stn, point))
        nai_met_vals.append(naiveMet(stn, point, 0))
        geo_met_vals.append(geoMet(stn, point))
        
        # adds in the point's shave and walk values
        average_shave /= 100
        const_tight_vals.append(average_shave)
        average_walk /= 100
        walk_vals.append(average_walk)
    
    return sph_met_vals, wil_met_vals, nai_met_vals, geo_met_vals, walk_vals, const_tight_vals

# stn = randomSchedule.fixStupidJSON(datastuff[0][57])
# print len(stn.keys())
# lists = randomPointsMetricsWalksShaves(stn, 200)
# geomVals = lists[3]
# walkVals = lists[4]
# plotting.justPlot(geomVals, 'geomVals', walkVals, 'walkVals', 'SOMETHING TITLEISH')

# print 'GEOMETRIC MEAN VALUES'
# for i in range(200):
#   print geomVals[i]

# print 'WALK VALUES'
# for i in range(200):
#   print walkVals[i]

# This function takes in a center of interest, a number of STNs to test with,
# a durability metric and a number of points
# It returns a list containing the durability values for the given center on
# the numSTN STNs tested, as well as a list containing the average durability
# values of the numPoints points tested for this STN
def centerDurabilityComparison(center, numSTN, metric, numPoints):

    centerMets = []
    averagePointsMets = []

    for i in range(numSTN):
        print i 
        stn = randomSchedule.fixStupidJSON(datastuff[0][i])
        CoI = None
        marg = 0.0
        index = 0

        if metric == 'sphMet': index = 0
        elif metric == 'wilMet': index = 1
        elif metric == 'naiMet': index = 2
        else: index = 3         
        
        average_point_durability = sum(randomPointsMetricsWalksShaves(stn, numPoints)[index])/numPoints

        if center == 'cheb':
            CoI = randomWalk.chebyshev(stn)[0]
        elif center == 'cent':
            CoI = randomSchedule.centroid(stn)
        elif center == 'greed':
            CoI = randomSchedule.avgMidPoint(stn)
        elif center == 'jim':
            CoI = randomSchedule.jimPoint(stn)
        
        if metric == 'sphMet':
            met = sphereMet(stn, CoI)
        elif metric == 'wilMet':
            met = wilsonCubeMet(stn, CoI)
        elif metric == 'naiMet':
            met = naiveMet(stn, CoI)
        else:
            met = geoMet(stn, CoI)
        
        centerMets.append(met)
        averagePointsMets.append(average_point_durability)
    return centerMets, averagePointsMets

# Code to call and test centerDurabilityComparison()
#
# lists = centerDurabilityComparison('jim', 30, 'geoMet', 100)
# margs = randomSchedule.listMargin(lists[0], lists[1])
# print float(sum(margs))/len(margs)
# plotting.justPlot(lists[0], 'centerMets', lists[1], 'averageMets', 'SOMETHING TITLEISH')

def centerPerformanceComparison(center, numSTN, eict, numPoints):

    centerWalksShaves = []
    averagePointsWalksShaves = []

    for i in range(numSTN):
        print i 
        stn = randomSchedule.fixStupidJSON(datastuff[0][i])
        CoI = None
        marg = 0.0
        index = 0
        walk = 0

        if eict == 'walk': index = 4
        elif eict == 'shave': index = 5         
        
        average_point_walk_shave = sum(randomPointsMetricsWalksShaves(stn, numPoints)[index])/numPoints

        if center == 'cheb':
            CoI = randomWalk.chebyshev(stn)[0]
        elif center == 'cent':
            CoI = randomSchedule.centroid(stn)
        elif center == 'greed':
            CoI = randomSchedule.avgMidPoint(stn)
        elif center == 'jim':
            CoI = randomSchedule.jimPoint(stn)
        
        if eict == 'walk':
            walk = randomWalk.perturb(stn, CoI, 1, 0, 0)
        elif eict == 'shave':
            walk = randomMelt.randomShave(stn, CoI)
        
        centerWalksShaves.append(walk)
        averagePointsWalksShaves.append(average_point_walk_shave)
    return centerWalksShaves, averagePointsWalksShaves

# Code to call and test centerPerformanceComparison()
#
# lists = centerPerformanceComparison('jim', 10, 'walk', 100)
# margs = randomSchedule.listMargin(lists[0], lists[1])
# print float(sum(margs))/len(margs)

def durabilityRandsSingleSTN(stn, metric, numPoints):
    durabilityMets = []
    walksList = []
    shavesList = []

    for i in range(numPoints):
        # print i
        sched = randomSchedule.simpleRandSched(stn)
        walk_val = 0.0
        shave_val = 0.0

        for i in range(100):
            walk_val += randomWalk.perturb(stn, sched, 1, 0, 0)
            shave_val += randomMelt.randomShave(stn, sched)
        
        walksList.append(walk_val/100)
        shavesList.append(shave_val/100)

        if metric == 'sphMet': 
            durabilityMets.append(sphereMet(stn, sched))
        elif metric == 'wilMet': 
            durabilityMets.append(wilsonCubeMet(stn, sched))
        elif metric == 'naiMet': 
            durabilityMets.append(naiveMet(stn, sched))
        else: 
            durabilityMets.append(geoMet(stn, sched))

    return durabilityMets, walksList, shavesList

# Code to call and test durabilityRandsSingleSTN()

# stn = randomSchedule.fixStupidJSON(datastuff[0][16065])
# stn = randomSchedule.fixStupidJSON(datastuff[0][0])
# lists = durabilityRandsSingleSTN(stn, 'geoMet', 20)
# plotting.justPlot(lists[0], 'durabilityMet', lists[1], 'randomWalkLength', 'FIGHT')

def durabilityRandsManySTN(metric, numSTN, numPoints):
    corr_coeff_list = []

    for i in range(numSTN):
        print "CLONING VIVA NUMBER " + str(i) 
        stn = randomSchedule.fixStupidJSON(datastuff[0][i])
        lists = durabilityRandsSingleSTN(stn, metric, numPoints)
        durabilityList = lists[0]
        walk_vals = lists[1]
        corr_coeff = np.corrcoef(durabilityList, walk_vals)[1][0]

        print corr_coeff
        corr_coeff_list.append(corr_coeff)
    
    return corr_coeff_list

# Code to call and test durabilityRandsManySTN()

# for i in range(50):
#   stn = randomSchedule.fixStupidJSON(datastuff[0][i])
#   jimPoint = randomSchedule.jimPoint(stn)
#   midPoint = randomSchedule.midPoint(stn)
#   print randomWalk.distance(jimPoint, midPoint) / randomSchedule.furthestBoundary(stn, midPoint)
# corrCoeffList = durabilityRandsManySTN('geoMet', 20, 100)
# print sum(corrCoeffList) / float(len(corrCoeffList))


def mpRandomMetricsHelper(tup):
    genstn = tup[0]
    numWalks = tup[1]
    dim = len(genstn.keys()) - 1

    # Compute the STN flexibility metric values for the current STN
    sphereMet = metricFunctions.sphericalMetric(genstn, dim, "spherical")
    burgerMet = metricFunctions.sphericalMetric(genstn, dim, "burgers")/(dim*(dim-1)/2)
    wilsonMet = metricFunctions.sphericalMetric(genstn, dim, "wilson")/dim
    naiveMet = metricFunctions.sphericalMetric(genstn, dim, "hell")
    
    # Compute the centers of interest for the current STN
    cheb = randomWalk.chebyshev(genstn)[0]
    cent = randomSchedule.centroid(genstn)
    avgMid = randomSchedule.avgMidPoint(genstn)
    
    # Lists to hold the disturbance lengths
    randWalks = []
    randShaves = []

    # Lists to hold the durability metric values
    randclosevals = []
    randcubevals = []
    randgeovals = []

    # Compute the disturbance lengths and durability metric values for 100
    # random schedules in the STN
    for i in range(100):
        walk = 0.0
        shave = 0.0
        curr = simpleRandSched(genstn)

        randclosevals.append(sphereMet(genstn, curr))
        randcubevals.append(wilsonCubeMet(genstn, curr))
        randgeovals.append(geoMet(genstn, curr))

        for j in range(numWalks):
            walk += randomWalk.perturb(genstn, curr, 1, 0, 0)
            shave += randomMelt.randomShave(genstn, curr)

        walk /= numWalks
        shave /= numWalks

        randWalks.append(walk)
        randShaves.append(shave)

    # Compute averages for the disturbance lengths and durability metric
    # values for the 100 STNs
    randAvg = sum(randWalks)/len(randWalks)
    randAvgShaves = sum(randShaves)/len(randShaves)
    randClose = sum(randclosevals)/len(randclosevals)
    randCube = sum(randcubevals)/len(randcubevals)
    randGeo = sum(randgeovals)/len(randgeovals)

    # Compute correlation coefficients for the disturbance lengths vs each of
    # our durability metrics
    walk_close_corr = np.corrcoef(randclosevals, randWalks)[1][0]
    walk_cube_corr = np.corrcoef(randcubevals, randWalks)[1][0]
    walk_geo_corr = np.corrcoef(randgeovals, randWalks)[1][0]
    shave_close_corr = np.corrcoef(randclosevals, randShaves)[1][0]
    shave_cube_corr = np.corrcoef(randcubevals, randShaves)[1][0]
    shave_geo_corr = np.corrcoef(randgeovals, randShaves)[1][0]

    # Placeholders for disturbance lengths for centers of interest
    chebAvg = 0.0
    chebAvgShaves = 0.0

    centAvg = 0.0
    centAvgShaves = 0.0

    midAvg = 0.0
    midAvgShaves = 0.0

    # Compute disturbance lengths for centers of interest multiple times
    for j in range(numWalks):
        chebAvg += randomWalk.perturb(genstn, cheb, 1, 0, 0)
        chebAvgShaves += randomMelt.randomShave(genstn, cheb)

        centAvg += randomWalk.perturb(genstn, cent, 1, 0, 0)
        centAvgShaves += randomMelt.randomShave(genstn, cent)

        midAvg += randomWalk.perturb(genstn, avgMid, 1, 0, 0)
        midAvgShaves += randomMelt.randomShave(genstn, avgMid)

    # Compute average disturbance lengths for centers of interest
    chebAvg /= numWalks
    chebAvgShaves /= numWalks

    centAvg /= numWalks
    centAvgShaves /= numWalks

    midAvg /= numWalks
    midAvgShaves /= numWalks

    # Compute durability metric values for centers of interest
    cheb_close = sphereMet(genstn, cheb)
    cent_close = sphereMet(genstn, cent)
    mid_close = sphereMet(genstn, avgMid)
    cheb_cube = wilsonCubeMet(genstn, cheb)
    cent_cube = wilsonCubeMet(genstn, cent)
    mid_cube = wilsonCubeMet(genstn, avgMid)
    cheb_geo = geoMet(genstn, cheb)
    cent_geo = geoMet(genstn, cent)
    mid_geo = geoMet(genstn, avgMid)

    print "TRIAL!! " + str(tup[2])

    return [genstn, dim, sphereMet, burgerMet, wilsonMet, naiveMet, cheb,cent,
            avgMid, randAvg, randAvgShaves, randClose, randCube, randGeo,
            walk_close_corr, walk_cube_corr, walk_geo_corr, shave_close_corr,
            shave_cube_corr, shave_geo_corr, chebAv, chebAvgShaves, centAvg,
            centAvgShaves, midAvg, midAvgShaves, cheb_close, cent_close,
            mid_close, cheb_cube, cent_cube, mid_cube, cheb_geo, cent_geo, 
            mid_geo]


##############
# walk50s = []
# walk100s = []
# walk500s = []
# walks = []

# for j in range(100):
#     stn = randomSchedule.fixStupidJSON(datastuff[0][j])
#     cheb = randomWalk.chebyshev(stn)[0]
#     walk = 0.0
#     walk50 = 0.0
#     walk100 = 0.0
#     walk500 = 0.0
#     for i in range(1000):
#         if i == 50 or i == 100 or i == 500:
#             if i == 50:
#                 walk50s.append(walk/i)
#             elif i == 100:
#                 walk100s.append(walk/i)
#             else:
#                 walk500s.append(walk/i)
#             print str(i) + " WALKS: " + str(walk/(i))
#         walk += randomWalk.perturb(stn, cheb, 1, 0, 0)
#     walks.append(walk/1000)
#     print "STN #" + str(j) + str(walk/1000)
# rat50 = 0.0
# rat100 = 0.0
# rat500 = 0.0
# print walk50s
# print walk100s
# print walk500s
# print walks
# walk50s = [49.78, 64.8, 22.74, 112.22, 41.54, 27.36, 29.16, 13.54, 30.12, 31.56, 148.74, 38.42, 101.74, 31.94, 52.88, 30.18, 21.38, 11.46, 10.62, 42.86, 22.74, 15.94, 107.26, 7.5, 59.12, 85.32, 321.74, 146.1, 0.16, 23.86, 145.32, 346.26, 452.62, 153.44, 7.98, 135.8, 24.64, 39.4, 57.9, 196.78, 12.36, 702.88, 3.96, 2.48, 285.6, 298.16, 368.94, 86.34, 62.76, 705.78, 27.94, 519.42, 45.3, 364.38, 13.06, 9.36, 26.28, 271.54, 98.54, 107.34, 8.26, 26.84, 416.52, 94.38, 275.1, 61.28, 13.04, 38.82, 12.22, 15.66, 29.78, 22.92, 24.5, 108.04, 13.58, 26.18, 21.14, 2.18, 542.42, 28.56, 133.18, 23.26, 168.36, 221.2, 4.78, 22.92, 34.46, 0.0, 1131.8, 14.08, 338.34, 33.4, 0.0, 16.14, 17.88, 0.58, 147.36, 51.46, 73.98, 132.8]
# walk100s = [47.66, 68.1, 22.81, 124.06, 38.79, 26.81, 28.87, 12.27, 31.01, 31.08, 133.77, 37.77, 99.16, 31.08, 56.05, 32.68, 21.47, 11.3, 10.31, 45.01, 23.46, 14.76, 107.14, 7.65, 58.41, 91.13, 277.22, 154.14, 0.11, 26.72, 147.53, 379.28, 454.19, 151.25, 8.35, 154.19, 25.04, 38.21, 59.5, 190.69, 12.6, 716.63, 3.54, 2.43, 273.59, 323.72, 362.87, 85.74, 74.6, 643.25, 27.6, 507.89, 43.0, 330.56, 12.53, 9.36, 28.94, 256.24, 91.71, 117.95, 8.37, 28.92, 409.16, 95.35, 292.14, 58.47, 14.04, 35.11, 12.12, 14.81, 29.49, 24.29, 24.69, 112.99, 13.98, 28.22, 20.87, 2.05, 551.66, 27.69, 144.61, 22.52, 161.3, 206.41, 5.22, 24.05, 34.13, 0.0, 1012.28, 12.84, 389.9, 34.66, 0.0, 16.39, 18.82, 0.57, 155.05, 50.92, 77.2, 128.69]
# walk500s = [47.8, 61.002, 23.248, 129.604, 39.144, 29.184, 26.824, 11.87, 30.876, 33.444, 129.034, 35.594, 96.744, 29.992, 58.062, 31.042, 20.4, 11.124, 9.006, 39.49, 23.584, 13.528, 98.464, 7.642, 57.074, 86.282, 286.478, 184.788, 0.072, 24.644, 142.736, 348.208, 484.566, 146.052, 9.138, 160.226, 29.848, 36.352, 58.344, 183.254, 12.616, 685.372, 3.534, 2.672, 278.748, 328.778, 401.252, 82.404, 75.512, 596.502, 29.114, 494.216, 44.452, 343.278, 11.68, 8.802, 29.488, 244.104, 82.592, 131.208, 7.658, 29.512, 393.716, 84.538, 315.694, 64.866, 14.698, 34.77, 11.628, 13.452, 29.72, 22.694, 23.922, 122.904, 13.436, 27.326, 23.384, 2.078, 549.776, 27.026, 146.906, 23.274, 156.856, 195.286, 5.438, 22.638, 31.732, 0.0, 877.216, 12.708, 382.672, 34.628, 0.0, 16.49, 17.544, 0.596, 149.348, 47.84, 79.866, 124.846]
# walks = [45.4, 62.96, 23.233, 126.231, 38.253, 29.088, 25.934, 11.442, 31.812, 34.361, 129.586, 35.361, 96.219, 30.103, 58.813, 29.922, 19.698, 10.769, 8.856, 38.418, 24.675, 13.676, 96.864, 7.65, 56.832, 88.214, 289.509, 191.401, 0.071, 24.115, 145.649, 354.867, 504.324, 146.243, 9.119, 157.105, 29.722, 36.019, 58.789, 179.268, 12.508, 677.017, 3.528, 2.682, 275.612, 338.11, 400.252, 83.053, 74.966, 590.097, 29.133, 488.585, 45.561, 335.404, 11.931, 8.941, 29.07, 257.839, 80.564, 131.785, 7.556, 30.618, 378.202, 85.741, 310.525, 62.982, 14.766, 35.114, 11.273, 13.292, 29.952, 23.081, 23.913, 121.988, 13.607, 26.62, 23.027, 2.116, 556.238, 28.058, 148.171, 23.329, 155.858, 195.98, 5.321, 22.765, 30.887, 0.0, 859.324, 12.681, 368.657, 35.064, 0.0, 16.299, 17.417, 0.585, 154.989, 47.914, 81.838, 126.087]
# j = 0
# for i in range(100):
#     if walks[i] != 0.0:
#         rat50 += walk50s[i]/walks[i]
#         rat100 += walk100s[i]/walks[i]
#         rat500 += walk500s[i]/walks[i]
#     else:
#         j += 1
# print rat50/(100 - j)
# print rat100/(100 - j)
# print rat500/(100 - j)

# w_rat10 = 0
# s_rat10 = 0
# w_rat50 = 0
# s_rat50 = 0
# w_rat100 = 0
# s_rat100 = 0
# for i in range(100):
#     stn = randomSchedule.fixStupidJSON(datastuff[0][i])
#     c10 = randomSchedule.centroid(stn, 10)
#     c50 = randomSchedule.centroid(stn, 50)
#     c100 = randomSchedule.centroid(stn, 100)
#     c500 = randomSchedule.centroid(stn, 500)
#     c10_walk = 0.0
#     c50_walk = 0.0
#     c100_walk = 0.0
#     c500_walk = 0.0
#     c10_walk = randomWalk.distance(c10, c500)
#     c50_walk = randomWalk.distance(c50, c500)
#     c100_walk = randomWalk.distance(c100, c500)
#     c500_walk = randomWalk.distance(c100, c500)
#     c10_shave = 0.0
#     c50_shave = 0.0
#     c100_shave = 0.0
#     c500_shave = 0.0
#     # for j in range(100):
#     #     c10_walk += randomWalk.perturb(stn, c10, 1, 0, 0)
#     #     c50_walk += randomWalk.perturb(stn, c50, 1, 0, 0)
#     #     c100_walk += randomWalk.perturb(stn, c100, 1, 0, 0)
#     #     c500_walk += randomWalk.perturb(stn, c500, 1, 0, 0)
#     #     c10_shave += randomMelt.randomShave(stn, c10)
#     #     c50_shave += randomMelt.randomShave(stn, c50)
#     #     c100_shave += randomMelt.randomShave(stn, c100)
#     #     c500_shave += randomMelt.randomShave(stn, c500)
#     print "WALK10 " + str(i) + " " + str(c500_walk/(c10_walk+0.0001))
#     print "SHAVE10 " + str(i) + " " + str(c500_shave/(c10_shave+0.0001))
#     w_rat10 += c500_walk/(c10_walk+0.0001)
#     s_rat10 += c500_shave/(c10_shave+0.0001)
#     print "WALK50 " + str(i) + " " + str(c500_walk/(c50_walk+0.0001))
#     print "SHAVE50 " + str(i) + " " + str(c500_shave/(c50_shave+0.0001))
#     w_rat50 += c500_walk/(c50_walk+0.0001)
#     s_rat50 += c500_shave/(c50_shave+0.0001)
#     print "WALK100 " + str(i) + " " + str(c500_walk/(c100_walk+0.0001))
#     print "SHAVE100 " + str(i) + " " + str(c500_shave/(c100_shave+0.0001))
#     w_rat100 += c500_walk/(c100_walk+0.0001)
#     s_rat100 += c500_shave/(c100_shave+0.0001)
# print str(w_rat10/100) + " WRAT10!!"
# print str(s_rat10/100) + " SRAT10!!"
# print str(w_rat50/100) + " WRAT10!!"
# print str(s_rat50/100) + " SRAT10!!"
# print str(w_rat100/100) + " WRAT10!!"
# print str(s_rat100/100) + " SRAT10!!"
# with open(os.path.abspath('./Data/FINAL_DATA_500.json'), 'r') as fp:
#     data = json.load(fp)

# stns = data[0][:500]
# walkCorrList = []
# shaveCorrList = []

# for k in range(500):
#     genstn = randomSchedule.fixStupidJSON(stns[k])
#     numWalks = 100
#     dim = len(genstn.keys()) - 1

#     # Compute the centers of interest for the current STN
#     # cheb = randomWalk.chebyshev(genstn)[0]
#     # cent = randomSchedule.centroid(genstn, 100)
#     # avgMid = randomSchedule.jimPoint(genstn)

#     # Lists to hold the disturbance lengths
#     randWalks = []
#     randShaves = []

#     # Lists to hold the durability metric values
    

#     # Compute the disturbance lengths and durability metric values for 100
#     # random schedules in the STN
#     for i in range(100):
#         walk = 0.0
#         shave = 0.0
#         curr = randomSchedule.simpleRandSched(genstn)

#         randfurthvals.append(randomSchedule.furthestBoundary(genstn, curr))

#         for j in range(numWalks):
#             walk += randomWalk.perturb(genstn, curr, 1, 0, 0)
#             shave += randomMelt.randomShave(genstn, curr)

#         walk /= numWalks
#         shave /= numWalks

#         randWalks.append(walk)
#         randShaves.append(shave)

#     # Compute averages for the disturbance lengths and durability metric
#     # values for the 100 STNs
#     randAvg = sum(randWalks)/len(randWalks)
#     randAvgShaves = sum(randShaves)/len(randShaves)
#     randFurth = sum(randfurthvals)/len(randfurthvals)

#     # Compute correlation coefficients for the disturbance lengths vs each of
#     # our durability metrics
#     try:
#         walk_furth_corr = np.corrcoef(randfurthvals, randWalks)[1][0]
#     except:
#         walk_furth_corr = 0.0

#     try:
#         shave_furth_corr = np.corrcoef(randfurthvals, randShaves)[1][0]
#     except:
#         shave_furth_corr = 0.0

#     walkCorrList.append(walk_furth_corr)
#     shaveCorrList.append(shave_furth_corr)
#     print k
    # Placeholders for disturbance lengths for centers of interest
    # chebAvg = 0.0
    # chebAvgShaves = 0.0

    # centAvg = 0.0
    # centAvgShaves = 0.0

    # midAvg = 0.0
    # midAvgShaves = 0.0

    # # Compute disturbance lengths for centers of interest multiple times
    # for i in range(numWalks):
    #     chebAvg += randomWalk.perturb(genstn, cheb, 1, 0, 0)
    #     chebAvgShaves += randomMelt.randomShave(genstn, cheb)

    #     centAvg += randomWalk.perturb(genstn, cent, 1, 0, 0)
    #     centAvgShaves += randomMelt.randomShave(genstn, cent)

    #     midAvg += randomWalk.perturb(genstn, avgMid, 1, 0, 0)
    #     midAvgShaves += randomMelt.randomShave(genstn, avgMid)

    # # Compute average disturbance lengths for centers of interest
    # chebAvg /= numWalks
    # chebAvgShaves /= numWalks

    # centAvg /= numWalks
    # centAvgShaves /= numWalks

    # midAvg /= numWalks
    # midAvgShaves /= numWalks

    # # Compute durability metric values for centers of interest
    # cheb_furth = randomSchedule.furthestBoundary(genstn, cheb)
    # cent_furth = randomSchedule.furthestBoundary(genstn, cent)
    # mid_furth = randomSchedule.furthestBoundary(genstn, avgMid)

# print sum(walkCorrList) / len(walkCorrList)
# print sum(shaveCorrList) / len(shaveCorrList)
