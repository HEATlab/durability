## This file contains code to convert a JSON file to a DC_STN.
import json
import sys

from dc_stn import *

## \fn loadSTN(fn)
#  \brief loads an STN from a JSON file
def loadSTNfromJSON(filename):
  # obtain data from the json file 
  with open(filename) as data_file:
    data = json.load(data_file)

  # Create a new STN
  outputSTN = DC_STN()

  # Add Node
  for node in data["nodes"]:
    outputSTN.addVertex(node["node_id"])

  # Add Edges
  for edge in data["constraints"]:
    i = int(edge["first_node"])
    j = int(edge["second_node"])
    Cij = float(edge["max_duration"])
    Cji = float(edge["min_duration"])

    outputSTN.addEdge(i, j, Cij)
    outputSTN.addEdge(j, i, -Cji)

    # If there is a distribution, then edge is contingent
    # Thus there is Lower and Upper edge
    if "distribution" in edge:
      outputSTN.addEdge(i, j, Cji, edgeType.LOWER, j)
      outputSTN.addEdge(j, i, -Cij, edgeType.UPPER, j)

  return outputSTN

if __name__ == '__main__':

  # handle command line input
  filename = sys.argv[-1]
  print loadSTNfromJSON(filename)