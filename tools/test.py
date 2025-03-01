# This file contains code from Polybots that helps us run tests
# Authors: Amy Huang and Liam Lloyd
# 
# File created by: Viva Ojha and Joon Lee
# Contact Info: vmojha@g.hmc.edu and joolee@g.hmc.edu

import copy
import math
import random
import randomWalk
import stnConverter
from numpy.linalg import norm
from pulp import *

TARGET_RATIO = .5
EPSILON = .03
NUM_SAMPLES = 100

# ============================================================================
# MATRIX FUNCTIONS
# ============================================================================
def minimal(STN):
    ''' Takes distance matrix of an STN,
        Returns the minimal form of STN
        Returns [] if STN is inconsistent
    '''

    n = len(STN)
    stn = copy.deepcopy(STN)
    # print stn
    for k in range(n):
        for i in range(n):
            for j in range(n):
                if stn[i][j] > stn[i][k] + stn[k][j]:
                    stn[i][j] = stn[i][k] + stn[k][j]
                if stn[i][j] + stn[j][i] < 0:
                    print stn
                    print "i: " + str(i)
                    print "j: " + str(j)
                    print "K: " + str(k)
                    print stn[i][j]
                    print stn[j][i]
                    print "Inconsistent STN: cycle sums to: ",
                    print stn[i][j]+stn[j][i]
                    return []

    return stn


def printMatrix(matrix, maxDigitLen=2):
    ''' Displays matrix in pretty form
    '''
    if matrix == []:
        print "matrix is empty\n"
        return
    elif not isinstance(matrix[0], list):
        print matrix
        return

    n = len(matrix)
    m = len(matrix[0])

    for i in range(n):
        for j in range(m):
            val = matrix[i][j]
            space = ""
            if val == float("inf") or val == -float("inf"):
                space = " " * (maxDigitLen - 1)
                print "x",
                print space,
            else:
                s = str(matrix[i][j])
                space = " " * (maxDigitLen - len(s))
                print str(matrix[i][j]),
                print space,
        print
    print
    print

# ============================================================================
# GENERATORS FOR STNS
# ============================================================================
def generateSTN(dim, typ, interval=5):
    ''' Generates concurrent, sequential, and half sequential STNs
        dim: the dimension of the STN
        interval: the size of every event's domain
    '''
    numRows = dim + 1
    stn = []
    if typ == "seq" or typ == "seq2":
        for i in range(numRows):
            row = [0] + [0] * i + [interval] * (dim - i)
            stn.append(row)

    elif typ == "con":
        for i in range(numRows):
            if i == 0:
                row = [0] + [interval] * dim
            else:
                row = [0] + [float("inf")]*(i-1) + [0] + \
                    [float("inf")]*(dim - i)
            stn.append(row)

    if typ == "seq2":
        rowsOfFives = dim / 2
        if dim % 2 == 1:
            rowsOfFives += 1

        for i in range(rowsOfFives):
            for j in range(dim/2):
                stn[i+1][numRows-j-1] = interval
        for i in range(dim/2):
            for j in range(rowsOfFives):
                stn[numRows-j-1][i+1] = interval

    return minimal(stn)


def generateRectSTN(dim, sideLens):
    ''' Generates an STN whose poltyope is a rectangular prism
        sideLens is an array of domain sizes for each axis
    '''
    numRows = dim + 1
    stn = []
    for i in range(numRows):
        if i == 0:
            row = [0] + sideLens
        else:
            row = [0] + [float("inf")]*(i-1) + [0] + [float("inf")]*(dim-i)
        stn.append(row)
    return minimal(stn)


def generateRandomSTN(dim, lowerBound, upperBound):
    ''' Takes in number of events, lower and upper bound
        Returns distance matrix with random numbers
    '''
    stn = []
    while True:
        attemptSTN = []
        for i in range(dim+1):
            row = []
            for j in range(dim+1):
                if i == j:
                    row.append(0)
                else:
                    value = random.randint(lowerBound, upperBound)
                    row.append(value)
            attemptSTN.append(row)

        mstn = minimal(attemptSTN)
        if mstn:
            stn = mstn
            break

    return stn


def generateNestedSTN(prevSTN):
    ''' Takes an STN, generates a new STN that is
        contained in it
        Assumes integer values; generates integer STN

        Limits decrease in constraint to a constant=10
        Any constraint has a .3 chance of decreasing
    '''
    limit = 10
    changeProb = 0.3

    stn = []
    n = len(prevSTN)
    while True:
        newSTN = copy.deepcopy(prevSTN)
        for i in range(n):
            for j in range(i+1, n):
                prob = random.random()
                if prob < changeProb:
                    # integer values
                    newSTN[i][j] = random.randint(
                        max(newSTN[i][j]-10, -newSTN[j][i]), newSTN[i][j])
                    newSTN[j][i] = random.randint(
                        max(newSTN[j][i]-10, -newSTN[i][j]), newSTN[j][i])
        mstn = minimal(newSTN)
        if mstn:
            stn = mstn
            break
    return stn

# ============================================================================
# MISCELLANEOUS UTILITIES
# ============================================================================
def reducedForm(STN):
    ''' given an STN in matrix form,
        finds all redundant edges in the
        graph and eliminates them, leaving
        an equivalent STN with the fewest
        explicit constraints'''

    reducedSTN = copy.deepcopy(STN)
    n = len(reducedSTN)  # get matrix dimension

    # create matrix that flags whether edges have been changed
    flags = []
    for i in range(n):
        flags.append([False] * n)

    # no edge from point to itself
    for i in range(n):
        reducedSTN[i][i] = float("inf")

    # delete redundant edges
    for k in range(n):
        for i in range(n):
            for j in range(n):
                if k != i and k != j:
                    if reducedSTN[i][j] >= reducedSTN[i][k] + reducedSTN[k][j]:
                        flags[i][j] = True

    # delete flagged edges
    for i in range(n):
        for j in range(n):
            if flags[i][j]:
                reducedSTN[i][j] = float("inf")

    return reducedSTN


def STNtoPolytopeMatrix(STN):
    '''Takes an STN in the matrix representation
       of its graph form and converts it to the
       matrix representation of the polytope it
       defines. Returns a matrix representing the
       constraints and a vector representing their
       upper bounds'''

    constraintMatrix = []
    boundsVector = []

    for i in range(len(STN)):
        for j in range(len(STN[i])):
            if not i == j and not (STN[i][j] == float("inf")):
                newRow = [0]*len(STN)

                # in the STN matrix, t_j - t_i <= STN[i][j]
                newRow[i] = -1
                newRow[j] = 1

                # delete the first column of the matrix, since it represents the
                # zero timepoint, which is only necessary in the graph form.
                del newRow[0]
                constraintMatrix.append(newRow)
                boundsVector.append(STN[i][j])

    return constraintMatrix, boundsVector


def fixValue(stn, fixedVar, fixedTime):
    ''' Takes an STN with n dimensions, and fixes
        one of its variables.
        Returns an n-1 dimensional matrix
    '''
    newSTN = copy.deepcopy(stn)

    # set event to the timepoint
    newSTN[0][fixedVar] = fixedTime
    newSTN[fixedVar][0] = -fixedTime
    crossSectionSTN = minimal(newSTN)

    # remove the extra dimension
    crossSectionSTN.pop(fixedVar)
    for row in crossSectionSTN:
        row.pop(fixedVar)

    return crossSectionSTN


def sphereVolume(radius, dim):
    ''' Returns the volume of a sphere, given radius and dimension 
    '''
    return math.pi ** (dim/2.0) * radius**dim / math.gamma(dim/2.0 + 1)

# ============================================================================
# CODE RELATED TO TESTING OF SCHEDULES INSIDE AN STN
# ============================================================================
def getRandomSchedule(minimalSTN):
    '''Takes a minimal STN and returns
       a random schedule satisfying it'''
    '''
    # scale the STN matrix up by 10^10 so we can deal with long ints instead 
    # of floats
    for i in range(len(minimalSTN)):
        for j in range(len(minimalSTN[i])):

    '''
    # keep track of indices still to be assigned
    indices = range(0, len(minimalSTN))
    randSchedule = [0]*(len(minimalSTN))

    for i in range(len(minimalSTN)):
        for j in range(len(minimalSTN[i])):
            minimalSTN[i][j] *= 10000000000

    while indices:  # while indices is non-empty
        # print minimalSTN
        indexToAssign = random.choice(indices)  # pick an index to assign
        # print "INDEX!!!" + str(indexToAssign)
        # delete the chosen index from the list of indices still to be assigned
        indices.remove(indexToAssign)
        # print(indexToAssign)
        x = random.random()*(minimalSTN[0][indexToAssign] + minimalSTN[indexToAssign][0]) - minimalSTN[indexToAssign][0]
        # print "0i"
        # print minimalSTN[0][indexToAssign]
        # print "i0"
        # print minimalSTN[indexToAssign][0]
        # print "weird rhing"
        # print x
        
        # timepoint = long(x)
        
        timepoint = long(x)

        # assign the index a random time that can satisfy the STN
        randSchedule[indexToAssign] = timepoint

        # recalculate the minimal STN with this event nailed down
        minimalSTN[indexToAssign][0] = -timepoint
        minimalSTN[0][indexToAssign] = timepoint
        minimalSTN = minimal(minimalSTN)
 
    for i in range(len(randSchedule)):
        randSchedule[i] = float(randSchedule[i])/10000000000

    return randSchedule


def randomStep(schedule, constraintMatrix, boundsVector):
    '''Takes in a point inside the polytope defined by an STN,
       chooses a direction from a uniform distribution and a
       point on the line segment in that direction from the
       starting point to the boundary of the polytope from a
       uniform distribution. Returns this point'''

    direction = random.randint(0, len(constraintMatrix[0])-1)
    maximum = float('inf')
    minimum = -float('inf')

    # find bounds to ensure the point selected is in our polytope
    for i in range(len(constraintMatrix)):
        if not constraintMatrix[i][direction] == 0:
            bound = (boundsVector[i] - sum([a*b for a, b in zip(
                constraintMatrix[i], schedule)]))/constraintMatrix[i][direction]
        else:
            bound = float('inf')
        if constraintMatrix[i][direction] > 0 and bound < maximum:
            maximum = bound
        elif constraintMatrix[i][direction] < 0 and bound > minimum:
            minimum = bound

    # select a point
    change = minimum + random.random()*(maximum-minimum)
    newSample = copy.deepcopy(schedule)
    newSample[direction] += change
    return newSample


def wiggle(schedule, perturbAll, distribution, parameter, constraintMatrix):
    '''Selects an event in the given schedule at random and
       adds a random number between -perturbationSize and
       perturbationSize to it, then returns the schedule'''

    if perturbAll:
        for i in range(0, len(constraintMatrix[0])):
            direction.append(random.normalvariate(
                random.random(), random.random()))
        directionNorm = norm(direction)
        direction = [i/directionNorm for i in direction]
    else:
        direction = [0]*len(constraintMatrix[0])
        # start at 1 so we don't move the 0 timepoint
        direction[random.randint(0, len(direction)-1)] = 1

    if distribution == "uniform":
        magnitude = (random.random()*2*parameter - parameter)
    if distribution == "normal":
        magnitude = random.normalvariate(0, parameter)
    if distribution == "exact":
        magnitude = parameter
    direction = [i*magnitude for i in direction]

    schedule = [i+j for i, j in zip(schedule, direction)]
    return schedule, direction


def inPolytope(schedule, constraintMatrix, boundsVector):
    '''Takes a schedule and a constraint matrix and bounding vector
       and determines whether the schedule is in the polytope defined
       by the matrix and the vector'''

    for i in range(len(constraintMatrix)):
        if sum([a*b for a, b in zip(constraintMatrix[i], schedule)]) > boundsVector[i]:
            return False
    return True


def firstNonZeroIndex(list):
    '''Takes a list and returns the index of the first nonzero value'''

    index = 0
    while list[index] == 0:
        index += 1
    return index


def resilience(STN, numSamples, perturbAll, distribution, parameter, printFailureAxes):
    '''Takes an STN in the matrix representation
       of its graph form and a number of samples
       to take. Takes numSamples uniformly distributed
       points in the polytope defined by the STN and
       perturbs them slightly according to a normal or
       uniform distribution along either all dimensions
       or only one, then checks whether they
       are still in the polytope. Returns the fraction
       of the sampled points that remain in the polytope
       after perturbation. If distribution is normal,
       parameter should be standard deviation, if distribution
       is uniform, parameter should be the radius of the interval.'''

    reducedSTN = reducedForm(STN)
    constraintMatrix, boundsVector = STNtoPolytopeMatrix(reducedSTN)
    minimalSTN = minimal(STN)
    startPoint = getRandomSchedule(minimalSTN)
    samples = [startPoint]
    directionDict = {}

    # take samples
    for i in range(numSamples-1):
        samples.append(randomStep(samples[-1], constraintMatrix, boundsVector))
    # perturb samples, check whether they are still in the polytope
    numPerturbedInPolytope = 0
    for sample in samples:
        perturbedSample, perturbationDirection = wiggle(
            sample, perturbAll, distribution, parameter, constraintMatrix)

        # track the axes along which things fail
        if inPolytope(perturbedSample, constraintMatrix, boundsVector):
            numPerturbedInPolytope += 1
        if printFailureAxes and not perturbAll and not inPolytope(perturbedSample, constraintMatrix, boundsVector):
            directionKey = firstNonZeroIndex(perturbationDirection)
            if directionKey in directionDict:
                directionDict[directionKey] += 1
            else:
                directionDict[directionKey] = 1

    if printFailureAxes and not perturbAll:
        for directionKey in directionDict:
            print directionDict[directionKey], " samples fell out of the bounding polytope along the ", directionKey, "vector"

    return float(numPerturbedInPolytope)/numSamples


def randomWalkk(STN, radius, samples):
    '''Takes an STN, generates a random satisfying
       schedule, counts the number of steps it can take
       in a random walk of with step size of radius 
       before falling out of the STN's
       solution space. Repeats this samples times and
       and averages the results.
    '''

    reducedSTN = reducedForm(STN)
    constraintMatrix, boundsVector = STNtoPolytopeMatrix(reducedSTN)
    minimalSTN = minimal(STN)
    stepSum = 0

    for i in range(samples):
        startPoint = getRandomSchedule(copy.deepcopy(minimalSTN))[1:]
        while inPolytope(startPoint, constraintMatrix, boundsVector):
            startPoint = wiggle(startPoint, False, "uniform",
                                radius, constraintMatrix)[0]
            stepSum += 1

    averageSteps = float(stepSum)/samples
    # averageIndividualFlex = averageSteps/(len(minimalSTN)-1) #divide out by dimension to get how much freedom each variable had on average

    return averageSteps

# ============================================================================
# STN FLEXIBILITY METRICS
# ============================================================================
def naiveFlex(d):
    ''' input:  stn, a distance matrix
        output: the naive flexibility metric  
    '''
    stn = stnConverter.dictToMatrix(d)
    n = len(stn)

    flex = 0
    for i in range(1, n):
        flex += stn[i][0] + stn[0][i]

    return flex


def naiveVolFlex(d):
    stn = stnConverter.dictToMatrix(d)
    n = len(stn)
    flex = 1
    for i in range(1, n):
        flex *= stn[i][0] + stn[0][i]
    return flex


def hunsbergerFlex(stn):
    ''' Takes distance matrix of STN,
        returns Hunsberger flexibility
    '''
    n = len(stn.keys())
    flex = naiveFlex(stn)
    stnn = stnConverter.dictToMatrix(stn)
    for i in range(1, n):
        for j in range(i+1, n):
            flex += stnn[i][j] + stnn[j][i]
    return flex


def wilsonFlex(stn):
    ''' input:  stn, a distance matrix
        output: the Wilson metric for flexibility
    '''
    stn = stnConverter.dictToMatrix(stn)
    n = len(stn)
    tplus = [0] * n
    tminus = [0] * n
    prob = LpProblem("LP", LpMaximize)

    # add variables and t- <= t+ constraints
    for i in range(0, n):
        tplus[i] = LpVariable("tp" + str(i))
        tminus[i] = LpVariable("tm" + str(i))
        prob += tminus[i] - tplus[i] <= 0

    # maximize all (t+ - t-)
    prob += sum([tplus[i]-tminus[i] for i in range(1, n)])

    # add constraints
    prob += tplus[0] == 0
    prob += tminus[0] == 0
    for i in range(n):
        for j in range(n):
            if i != j and stn[i][j] != float("inf"):
                prob += tplus[i] - tminus[j] <= stn[i][j]

    status = prob.solve()

    flex = sum([value(x) for x in tplus]) - sum([value(x) for x in tminus])
    return flex

# WARNING: THIS TAKES FOREVER


def pointXRadius(stn):
    '''Using sampling within the polytope defined
       by the stn and a region around it defined by
       its surface area to find its surface area to
       volume ratio'''

    stn = stnConverter.dictToMatrix(stn)
    minimall = minimal(stn)
    largestInterval = 0
    for i in range(len(minimall[0])):
        if minimall[i][0] - minimall[0][i] > largestInterval:
            largestInterval = minimall[i][0] - minimall[0][i]

    radius = largestInterval/2.0
    ratio = 0
    leftRatio = float('inf')
    leftRadius = 0
    rightRatio = 0
    rightRadius = largestInterval
    while ratio < TARGET_RATIO-EPSILON or ratio > TARGET_RATIO+EPSILON:
        ratio = resilience(stn, NUM_SAMPLES, False, "exact", radius, False)
        if ratio > TARGET_RATIO:
            leftRadius = radius
        else:
            rightRadius = radius
        radius = (leftRadius+rightRadius)/2.0
    return radius


def spherical(STN, dim):
    rad = randomWalk.chebyshev(STN)[1]
    return sphereVolume(rad, dim) ** (1.0/dim)
# matrixCube = stnConverter.listToMatrix(randomWalk.cubestn)
# sched = getRandomSchedule(minimal(matrixCube))
# print sched
