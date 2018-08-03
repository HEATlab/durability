## This file contains code that converts the Robot Brunch STN format
## to the Polybots STN format

# We import this file because it contains code to convert the robot
# brunch STN format to inequalities, which are easy to deal with
import randomWalk 

def listToMatrix(S):
	inequalDict = getInequal(S)

	# initialize a default matrix of distances between nodes
	n = len(inequalDict.keys())
	matrixStn = [[float("inf")] * n] * n

	for i in range(n):
		currRow = matrixStn[i]
		currEdges = inequalDict[i]

		# Distance from node to itself is 0
		currRow[i] = 0.0

		# Loop through all edges
		for edge in currEdges:
			start = edge[0]
			end = edge[1]
			constraint = edge[2]

			# Fill in appropriate entry of distance matrix
			if i == start:
				currRow[end] = constraint
