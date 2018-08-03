import os
import sys
import json
import math
import copy
import itertools
import numpy as np
import multiprocessing as mp
import random

import randomSchedule
import randomWalk
import stnConverter
import plotting 
import randomMelt
import test
import metricFunctions
import schedMetrics

def mpHelper(tup):
    genstn = tup[0]
    numWalks = tup[1]
    dim = len(genstn.keys()) - 1

    # Compute the STN flexibility metric values for the current STN
    # Hunsberger, Wilson, and Naive are all adjusted to account for the
    # different dimensions
    sphereMet = metricFunctions.sphericalMetric(genstn, dim, "spherical")
    burgerMet = metricFunctions.sphericalMetric(genstn, dim, "burgers")/(dim*
                                                (dim-1)/2)
    wilsonMet = metricFunctions.sphericalMetric(genstn, dim, "wilson")/dim
    naiveMet = metricFunctions.sphericalMetric(genstn, dim, "hell")/dim
    
    # Compute the centers of interest for the current STN
    cheb = randomWalk.chebyshev(genstn)[0]
    cent = randomSchedule.centroid(genstn, 100)
    avgMid = randomSchedule.jimPoint(genstn)
    
    # Lists to hold the disturbance lengths
    randWalks = []
    randShaves = []

    # Lists to hold the durability metric values
    randclosevals = []
    randfurthvals = []
    # randcubevals = []
    randgeovals = []

    # Compute the disturbance lengths and durability metric values for 100
    # random schedules in the STN
    pointList = randomMelt.uniformRandSched(genstn)

    for i in range(100):
        walk = 0.0
        shave = 0.0
        # curr = randomSchedule.simpleRandSched(genstn)
        # print len(pointList)
        curr = pointList[i]

        randclosevals.append(schedMetrics.sphereMet(genstn, curr))
        randfurthvals.append(randomSchedule.furthestBoundary(genstn, curr))
        # randcubevals.append(schedMetrics.wilsonCubeMet(genstn, curr))
        randgeovals.append(schedMetrics.geoMet(genstn, curr))

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
    randFurth = sum(randfurthvals)/len(randfurthvals)
    # randCube = sum(randcubevals)/len(randcubevals)
    randGeo = sum(randgeovals)/len(randgeovals)

    # Compute correlation coefficients for the disturbance lengths vs each of
    # our durability metrics
    try:
        walk_close_corr = np.corrcoef(randclosevals, randWalks)[1][0] 
    except: 
        walk_close_corr = 0.0
    try:
        walk_furth_corr = np.corrcoef(randfurthvals, randWalks)[1][0]
    except:
        walk_furth_corr = 0.0
    # try:
    #     walk_cube_corr = np.corrcoef(randcubevals, randWalks)[1][0]
    # except:
    #     walk_cube_corr = 0.0
    try:
        walk_geo_corr = np.corrcoef(randgeovals, randWalks)[1][0]
    except:
        walk_geo_corr = 0.0

    try:    
        shave_close_corr = np.corrcoef(randclosevals, randShaves)[1][0]
    except: 
        shave_close_corr = 0.0
    try:
        shave_furth_corr = np.corrcoef(randfurthvals, randShaves)[1][0]
    except:
        shave_furth_corr = 0.0
    # try:
    #     shave_cube_corr = np.corrcoef(randcubevals, randShaves)[1][0]
    # except:
    #     shave_cube_corr = 0.0
    try:
        shave_geo_corr = np.corrcoef(randgeovals, randShaves)[1][0]
    except:
        shave_geo_corr = 0.0

    # walk_close_corr = 0.0
    # walk_cube_corr = 0.0 
    # walk_geo_corr = 0.0 
    # shave_close_corr = 0.0
    # shave_cube_corr = 0.0 
    # shave_geo_corr = 0.0
    if math.isnan(walk_close_corr):
        walk_close_corr = 0.0
    if math.isnan(walk_furth_corr):
        walk_furth_corr = 0.0
    if math.isnan(walk_geo_corr):
        walk_geo_corr = 0.0
    
    if math.isnan(shave_close_corr):
        shave_close_corr = 0.0
    if math.isnan(shave_furth_corr):
        shave_furth_corr = 0.0
    if math.isnan(shave_geo_corr):
        shave_geo_corr = 0.0

    # Placeholders for disturbance lengths for centers of interest
    chebAvg = 0.0
    chebAvgShaves = 0.0

    centAvg = 0.0
    centAvgShaves = 0.0

    midAvg = 0.0
    midAvgShaves = 0.0

    # Compute disturbance lengths for centers of interest multiple times
    for i in range(numWalks):
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
    cheb_close = schedMetrics.sphereMet(genstn, cheb)
    cent_close = schedMetrics.sphereMet(genstn, cent)
    mid_close = schedMetrics.sphereMet(genstn, avgMid)
    cheb_furth = randomSchedule.furthestBoundary(genstn, cheb)
    cent_furth = randomSchedule.furthestBoundary(genstn, cent)
    mid_furth = randomSchedule.furthestBoundary(genstn, avgMid)
    # cheb_cube = schedMetrics.wilsonCubeMet(genstn, cheb)
    # cent_cube = schedMetrics.wilsonCubeMet(genstn, cent)
    # mid_cube = schedMetrics.wilsonCubeMet(genstn, avgMid)
    cheb_geo = schedMetrics.geoMet(genstn, cheb)
    cent_geo = schedMetrics.geoMet(genstn, cent)
    mid_geo = schedMetrics.geoMet(genstn, avgMid)

    print "TRIAL!! " + str(tup[2])

    return [genstn, dim,  cheb, cent, avgMid, sphereMet, burgerMet, wilsonMet,
            naiveMet, randAvg, randAvgShaves, randClose, randFurth, randGeo,
            walk_close_corr, walk_furth_corr, walk_geo_corr, shave_close_corr,
            shave_furth_corr, shave_geo_corr, chebAvg, chebAvgShaves, centAvg,
            centAvgShaves, midAvg, midAvgShaves, cheb_close, cent_close,
            mid_close, cheb_furth, cent_furth, mid_furth,  cheb_geo, cent_geo,
            mid_geo]

# STNs = []
# for i in range(30000):
#     dim = random.randint(2, 10)
#     genstn = stnConverter.matrixToDict(test.generateSTN(dim, 'seq'))
#     genstn = stnConverter.matrixToDict(test.minimal(stnConverter.dictToMatrix(genstn)))
#     STNs.append(genstn)

# with open(os.path.abspath('./Data/' + 'ALL_OF_THE_STICK_STNS_30000' + '.json'), 'w') as fp:
#     json.dump(STNs, fp)

# ============================================================================
# Code to run multiprocessing
# ============================================================================
# with open(os.path.abspath('./Data/' + 'ALL_OF_THE_STNS_30000' + '.json'), 'r') as fp:
#     datastuff = json.load(fp)

# with open(os.path.abspath('./Data/' + 'ALL_OF_THE_STICK_STNS_20' + '.json'), 'r') as fp:
#     datastuff = json.load(fp)

# STNS = datastuff
# procPool = mp.Pool(2)
# STNS = [stnConverter.matrixToDict(test.minimal(stnConverter.dictToMatrix(
#         randomSchedule.fixStupidJSON(stn)))) for stn in STNS]

# tups = []

# for i in range(len(STNS)):
#     tups.append((STNS[i], 100, i))

# results = procPool.map(mpHelper, tups)

# lol = [[0 for i in range(len(results))] for j in range(len(results[0]))]
# for i in range(len(lol)):
#     for j in range(len(lol[0])):
#         lol[i][j] = results[j][i]

# with open(os.path.abspath('./Data/FINAL_DATA_30000_3.json'), 'w') as fp:
#     json.dump(lol, fp)

# with open(os.path.abspath('./Data/STICK_DATA_20.json'), 'w') as fp:
#     json.dump(lol, fp)
# ============================================================================
# Code to Plot Single STN durability metrics
# ============================================================================

# genstn = randomSchedule.fixStupidJSON(datastuff[13214])
# numWalks = 100
# randclosevals = []
# randWalks = []

# for i in range(100):
#     walk = 0.0
#    # curr = randomSchedule.simpleRandSched(genstn)
#     curr = randomSchedule.uniformRandSched(genstn)
#     randclosevals.append(schedMetrics.sphereMet(genstn, curr))

#     for j in range(numWalks):
#         walk += randomWalk.perturb(genstn, curr, 1, 0, 0)

#     walk /= numWalks

#     randWalks.append(walk)

# plotting.justPlot(randclosevals, 'Closest Boundary Metric Value', randWalks, 
#                   'EI Durability', '')

# ============================================================================
# Code to Manipulate and Plot Data
# ============================================================================

with open(os.path.abspath('./Data/FINAL_DATA_30000_3.json'), 'r') as fp:
    data = json.load(fp)
# # nadj = [data[8][i]/float(data[1][i]) for i in range(len(data[0]))]
# plotting.justPlot(data[8], 'Naive Metric Value', data[10], 
#                   'Random Shave Length', '')

# ============================================================================
# Code for Correlation Coefficients
# ============================================================================

# print 'Correlation for Walk VS Closest Boundary'
# print sum(data[14]) / len(data[14])

# print 'Correlation for Walk VS Furthest Boundary'
# print sum(data[15]) / len(data[15])

# print 'Correlation for Walk VS Geometric Mean'
# print sum(data[16]) / len(data[16])

# print 'Correlation for Shave VS Closest Boundary'
# print sum(data[17]) / len(data[17])

# print 'Correlation for Shave VS Furthest Boundary'
# print sum(data[18]) / len(data[18])

# print 'Correlation for Shave VS Geometric Mean'
# print sum(data[19]) / len(data[19])

# ============================================================================
# Code for Correlation Coefficients by Dimension
# ============================================================================

# dimList = [2, 3, 4, 5, 6, 7, 8, 9, 10]
# dimCorrList = []

# list2 = []
# list3 = []
# list4 = []
# list5 = []
# list6 = []
# list7 = []
# list8 = []
# list9 = []
# list10 = []

# corr = data[19]

# for i in range(30000):
#     if data[1][i] == 2:
#         list2.append(corr[i])
#     if data[1][i] == 3:
#         list3.append(corr[i])
#     if data[1][i] == 4:
#         list4.append(corr[i])
#     if data[1][i] == 5:
#         list5.append(corr[i])
#     if data[1][i] == 6:
#         list6.append(corr[i])
#     if data[1][i] == 7:
#         list7.append(corr[i])
#     if data[1][i] == 8:
#         list8.append(corr[i])
#     if data[1][i] == 9:
#         list9.append(corr[i])
#     if data[1][i] == 10:
#         list10.append(corr[i])

# dimCorrList.append(sum(list2) / len(list2))
# dimCorrList.append(sum(list3) / len(list3))
# dimCorrList.append(sum(list4) / len(list4))
# dimCorrList.append(sum(list5) / len(list5))
# dimCorrList.append(sum(list6) / len(list6))
# dimCorrList.append(sum(list7) / len(list7))
# dimCorrList.append(sum(list8) / len(list8))
# dimCorrList.append(sum(list9) / len(list9))
# dimCorrList.append(sum(list10) / len(list10))

# plotting.justPlot(dimList, 'dimList', dimCorrList, 'dimCorrList', 'title')

# ============================================================================
# Code for Correlation Coefficients across STNs
# ============================================================================

STNs = data[0]

sphereMets = []

randclosevals = []
randfurthvals = []
randgeovals = []

randWalks = []
randShaves = []

for i in range(len(STNs)):
    print 'STN ' + str(i)

    genstn = randomSchedule.fixStupidJSON(STNs[i])
    dim = len(genstn.keys()) - 1
    sphereMets.append(metricFunctions.sphericalMetric(genstn, dim, "spherical"))
    # curr = randomSchedule.simpleRandSched(genstn)
    curr = randomSchedule.uniformRandSched(genstn)[99]

    randclosevals.append(schedMetrics.sphereMet(genstn, curr))
    randfurthvals.append(randomSchedule.furthestBoundaryFix(genstn, curr))
    randgeovals.append(schedMetrics.geoMetFix(genstn, curr))

    walk = 0.0
    shave = 0.0
    
    for j in range(100):
        walk += randomWalk.perturb(genstn, curr, 1, 0, 0)
        shave += randomMelt.randomShave(genstn, curr)

    walk /= 100.0
    shave /= 100.0

    randWalks.append(walk)
    randShaves.append(shave)

print 'Walk Correlations'
plotting.justPlot(sphereMets, 'sphereMets', randWalks, 'randWalks', '')
plotting.justPlot(randclosevals, 'randClose', randWalks, 'randWalks', '')
plotting.justPlot(randfurthvals, 'randFurth', randWalks, 'randWalks', '')
plotting.justPlot(randgeovals, 'randGeo', randWalks, 'randWalks', '')

print 'Shave Correlations'
plotting.justPlot(sphereMets, 'sphereMets', randShaves, 'randShaves', '')
plotting.justPlot(randclosevals, 'randClose', randShaves, 'randShaves', '')
plotting.justPlot(randfurthvals, 'randFurth', randShaves, 'randShaves', '')
plotting.justPlot(randgeovals, 'randGeo', randShaves, 'randShaves', '')

# ============================================================================
# Code for Correlation Coefficients for Single STN
# ============================================================================

# STNs = data[0]

# walk_furth_corrs = []
# shave_furth_corrs = []


# for i in range(len(STNs)):
#     print 'STN ' + str(i)

#     genstn = randomSchedule.fixStupidJSON(STNs[i])
#     dim = len(genstn.keys()) - 1
    
#     randfurthvals = []
#     randWalks = []
#     randShaves = []

#     pointList = randomMelt.uniformRandSched(genstn)

#     for i in range(100):
#         walk = 0.0
#         shave = 0.0
#         # curr = randomSchedule.simpleRandSched(genstn)
#         curr = pointList[i]
#         randfurthvals.append(randomSchedule.closestBoundary(genstn, curr))

#         for j in range(100):
#             walk += randomWalk.perturb(genstn, curr, 1, 0, 0)
#             shave += randomMelt.randomShave(genstn, curr)

#         walk /= 100.0
#         shave /= 100.0

#         randWalks.append(walk)
#         randShaves.append(shave)

#     # Compute averages for the disturbance lengths and durability metric
#     # values for the 100 STNs

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

#     if math.isnan(walk_furth_corr):
#         walk_furth_corr = 0.0
    
#     if math.isnan(shave_furth_corr):
#         shave_furth_corr = 0.0

#     walk_furth_corrs.append(walk_furth_corr)
#     shave_furth_corrs.append(shave_furth_corr)

# print sum(walk_furth_corrs) / len(walk_furth_corrs)
# print sum(shave_furth_corrs) / len(shave_furth_corrs)

# ============================================================================
# Code for Margins Between Lists
# ============================================================================

# print len(data[0])
# randomSchedule.getMarg(data, 20, 22)
# randomSchedule.getMarg(data, 21, 23)
# randomSchedule.getMarg(data, 20, 24)
# randomSchedule.getMarg(data, 21, 25)
# randomSchedule.getMarg(data, 22, 24)
# randomSchedule.getMarg(data, 23, 25)
# margins = randomSchedule.listMargin(data[23], data[25])
# print sum(margins) / float(len(margins))


# ============================================================================
# Code for Margins of Disturbance Lengths
# ============================================================================

# print 'RANDOM WALK MARGINS'
# # Cheb random walk vs randos
# print 'Cheb'
# margins = randomSchedule.listMargin(data[20], data[9])
# print sum(margins) / float(len(margins))

# # Centroid random walk vs randos
# print 'Cent'
# margins = randomSchedule.listMargin(data[22], data[9])
# print sum(margins) / float(len(margins))

# # Greedy random walk vs randos
# print 'Greedy'
# margins = randomSchedule.listMargin(data[24], data[9])
# print sum(margins) / float(len(margins))

# print 'RANDOM SHAVE MARGINS'
# # Cheb random walk vs randos
# print 'Cheb'
# margins = randomSchedule.listMargin(data[21], data[10])
# print sum(margins) / float(len(margins))

# # Centroid random walk vs randos
# print 'Cent'
# margins = randomSchedule.listMargin(data[23], data[10])
# print sum(margins) / float(len(margins))

# # Greedy random walk vs randos
# print 'Greedy'
# margins = randomSchedule.listMargin(data[25], data[10])
# print sum(margins) / float(len(margins))

# ============================================================================
# Code for Margins of Durability Metrics
# ============================================================================

# print 'CLOSEST BOUND MARGINS'
# # Cheb random walk vs randos
# print 'Cheb'
# margins = randomSchedule.listMargin(data[26], data[11])
# print sum(margins) / float(len(margins))

# # Centroid random walk vs randos
# print 'Cent'
# margins = randomSchedule.listMargin(data[27], data[11])
# print sum(margins) / float(len(margins))

# # Greedy random walk vs randos
# print 'Greedy'
# margins = randomSchedule.listMargin(data[28], data[11])
# print sum(margins) / float(len(margins))


# print 'FURTHEST BOUND MARGINS'
# # Cheb random walk vs randos
# print 'Cheb'
# margins = randomSchedule.listMargin(data[29], data[12])
# # margins[174] = 0.0
# print sum(margins) / float(len(margins))

# # Centroid random walk vs randos
# print 'Cent'
# margins = randomSchedule.listMargin(data[30], data[12])
# print sum(margins) / float(len(margins))

# # Greedy random walk vs randos
# print 'Greedy'
# margins = randomSchedule.listMargin(data[31], data[12])
# print sum(margins) / float(len(margins))


# print 'GEOMETRIC MEAN MARGINS'
# # Cheb random walk vs randos
# print 'Cheb'
# margins = randomSchedule.listMargin(data[32], data[13])
# # print margins
# print sum(margins) / float(len(margins))

# # Centroid random walk vs randos
# print 'Cent'
# margins = randomSchedule.listMargin(data[33], data[13])
# print sum(margins) / float(len(margins))

# # Greedy random walk vs randos
# print 'Greedy'
# margins = randomSchedule.listMargin(data[34], data[13])
# print sum(margins) / float(len(margins))

