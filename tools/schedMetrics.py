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
