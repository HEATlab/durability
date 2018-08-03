# This file contains code that converts between the matrix STN format used
# in the Polybots code from Summer 2017
# and the inequality dictionary format that is based on the Robot Brunch code
# Authors: Viva Ojha and Joon Lee
# Contact Info: vmojha@g.hmc.edu and joolee@g.hmc.edu

# We import this file because it contains code to convert the robot
# brunch STN format to inequalities, which are easy to deal with
import randomWalk

# This function converts robot brunch STNs to Polybots STNs


def listToMatrix(S):
    inequalDict = randomWalk.getInequal(S)

    # initialize a default matrix of distances between nodes
    n = len(inequalDict.keys())
    matrixStn = [[float("inf")] * n for i in range(n)]

    for i in range(n):
        currRow = matrixStn[i]
        currEdges = inequalDict[i]

        # Distance from node to itself is 0
        currRow[i] = 0.0

        # Loop through all edges
        for edge in currEdges:
            start = edge[0]
            end = edge[1]
            constraint = float(edge[2])

            # Fill in appropriate entry of distance matrix
            if i == start:
                currRow[end] = constraint

    return matrixStn

# This function converts our STN format to Polybots STN format


def dictToMatrix(D):
    # initialize a default matrix of distances between nodes
    n = len(D.keys())
    matrixStn = [[float("inf")] * n for i in range(n)]

    for i in range(n):
        currRow = matrixStn[i]
        currEdges = D[i]

        # Distance from node to itself is 0
        currRow[i] = 0.0

        # Loop through all edges
        for edge in currEdges:
            start = edge[0]
            end = edge[1]
            constraint = float(edge[2])

            # Fill in appropriate entry of distance matrix
            if i == start:
                currRow[end] = constraint

    return matrixStn

# This function converts the Polybots STN format to our STN format


def matrixToDict(M):

    # initialize dictionary of inequalities and empty lists for the edges
    # relating to each node
    inequalDict = {}

    for i in range(len(M)):
        inequalDict[i] = []

    # loop through each node and fill in edges relating to it appropriately
    for i in range(len(M)):
        currRow = M[i]

        # loop through all the possible end nodes
        for j in range(len(currRow)):

            # only include nodes in dictionary if the connections are
            # nontrivial
            if currRow[j] != float("inf") and i != j:
                # append the edge to the lists of the start node and the end
                # node
                inequalDict[i].append([i, j, currRow[j]])
                inequalDict[j].append([i, j, currRow[j]])

    return inequalDict
