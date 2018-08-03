from sympy import *
from sympy.functions import Min as min
from simpleLinearDecoupling import floydWarshallReduction
import json, sys, os, pprint

## \file
#
#  \brief currently runs Floyd-Warshall symbolically in terms of makespan
#  \details doesn't work -- failed idea

## \fn floydWarshallReduction(verts,edges)
#
#  \brief Takes vertices and edges of an STN loaded from a JSON file and runs
#  the Floyd Warshall algorithm to reduce the constraints to only include
#  the valid solution space.
#
#  \returns the new verts and edges of the minimal STN
#
#  \warning This function in its current state modifies the input verts and
#  edges. Don't assume that the input variables are unchanged. They won't be.
def floydWarshallReduction(verts,edges):
  numVerts = len(verts) + 1
  stnMatrix = [[float('inf') for i in range(numVerts)] for j in range(numVerts)]

  for vert in verts:
    i = vert['node_id']
    stnMatrix[0][i] = vert['max_domain']
    stnMatrix[i][0] = -vert['min_domain']

  for edge in edges:
    i = edge['first_node']
    j = edge['second_node']
    minTime = edge['min_duration']
    maxTime = edge['max_duration']
    if maxTime == 'inf':
      maxTime = float('inf')
    stnMatrix[i][j] = maxTime
    stnMatrix[j][i] = -minTime

  for k in range(numVerts):
    for i in range(numVerts):
      for j in range(numVerts):
        stnMatrix[i][j] = min(stnMatrix[i][j],stnMatrix[i][k]+stnMatrix[k][j])
        print "\n\n", stnMatrix[i][j], "\n\n"
  for vert in verts:
    i = vert['node_id']
    vert['max_domain'] = stnMatrix[0][i]
    vert['min_domain'] = -stnMatrix[i][0]

  for edge in edges:
    i = edge['first_node']
    j = edge['second_node']
    edge['max_duration'] = stnMatrix[i][j]
    edge['min_duration'] = -stnMatrix[j][i]

  return verts,edges

## \fn parseCommandOptions
#  \brief parses the command line options and returns the parsed input
#  \details currently returns the json filename
def parseCommandOptions():
  usage = "USAGE: python symbolicSTNdecoupler.py JSONFILE"

  if len(sys.argv) != 2:
    if len(sys.argv) > 2:
      print "Too many arguments"
    else:
      print "Not enough arguments"
    print usage
    sys.exit(1)

  if not os.path.isfile(sys.argv[1]):
    print "Invalid filename: {}".format(sys.argv[1])
    print usage
    sys.exit(1)

  try:
    with open(sys.argv[1],"r") as f:
      STN = json.load(f)
      return STN
  except ValueError:
    print "File not in json format: {}".format(sys.argv[1])
    print usage
    sys.exit(1)
  except:
    raise

## \fn main
#  \brief loads a json file, replaces the makespan with a sympy variable, and
#  runs Floyd-Warshall
def main():
  STN = parseCommandOptions()
  makespan = Symbol("makespan_val", positive = True)
  decouplingPt = Symbol("x", positive = True)
  for node in STN["nodes"]:
    if node["node_id"] == 5:
      node["max_domain"] = decouplingPt
    elif node["node_id"] == 8:
      node["min_domain"] = decouplingPt
  #  node["max_domain"] = makespan
  floydWarshallReduction(STN["nodes"],STN["constraints"])
  p = pprint.PrettyPrinter(indent = 2)
  p.pprint(STN)


if __name__ == "__main__":
  main()
