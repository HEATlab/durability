import json, sys
from pulp import *

## \file flexDecomposition.py
#
#  \brief Solves a modified version of the Wilson et al decomposition that
#  maximizes the flexibility in contingent timepoints
#
#  \note You'll need pulp to be able to run this. On a linux machine run
#  `sudo pip install pulp` or `sudo easy_install -U pulp`.
#  THEN RUN `sudo pulptest`, otherwise it won't work.

## \fn main
#  \brief Sets up and solves the LP
def main():
  data = []

  with open(sys.argv[1]) as f:
    data = json.load(f)



  originalVerts = data["nodes"]
  originalEdges = data["constraints"]

  doubleSTN = [LpVariable("Z'"), LpVariable("t0Up"), LpVariable("t0Lo")]
  numVerts = 0
  for vert in originalVerts:
    numVerts += 1
    doubleSTN.append(LpVariable("t" + str(vert["node_id"]) + "Up"))
    doubleSTN.append(LpVariable("t" + str(vert["node_id"]) + "Lo"))

  print doubleSTN;


  prob = LpProblem("STN Flex Finder",LpMaximize)

  intervalSum = doubleSTN[1] - doubleSTN[2]
  intervalIndex = 3
  for x in range(0,len(originalVerts)):

    # This if statement weights the flexibility to cluster around uncertain verticies,
    # comment it out if you just want the straignt up flexibility
    if(originalVerts[x]["start"] == True and originalVerts[x]["local_id"] != 0):
      intervalSum += (doubleSTN[intervalIndex] - doubleSTN[intervalIndex + 1])
      intervalIndex += 2
      continue;

    intervalIndex+=2;


  prob += intervalSum, "The sum of the intervals, weighted to prioratize uncertain intervals"

  print prob
  vertStart = 1
  for x in range(0,1+len(originalVerts)):
    prob += doubleSTN[vertStart] >= doubleSTN[vertStart+1]
    #print "adding constraint that " + doubleSTN[vertStart].name + " >= " +  doubleSTN[vertStart + 1].name
    vertStart += 2


  prob += doubleSTN[1] == 0
  prob += doubleSTN[2] == 0
  prob += doubleSTN[0] == 0
  #add the constraints imposed by edges from the origin
  #We know that each event must start at least some time after the origin,
  #so we make sure that the lower bound of the event is at least that amount
  #of time after the latest end of the origin (0)
  variableIndex = 4
  vertIndex = 0
  for vert in originalVerts:
    prob += doubleSTN[1] - doubleSTN[variableIndex] <= originalVerts[vertIndex]["min_domain"]
    print "adding constraint that " + doubleSTN[1].name + " - " +  doubleSTN[variableIndex].name + " < " + str(originalVerts[vertIndex]["min_domain"])

    variableIndex += 2
    vertIndex += 1


  #add the constraints imposed by edges from the origin
  #We also know that the longest time from the origin to any vertex,
  #so we make sure the earliest time the origin can occur (0) is within
  #that amount of time of the latest any event can occur
  variableIndex = 3
  vertIndex = 0
  for vert in originalVerts:
    prob += doubleSTN[variableIndex] - doubleSTN[2] <= originalVerts[vertIndex]["max_domain"]
    print "adding constraint that " + doubleSTN[variableIndex].name + " - " +  doubleSTN[2].name + " < " + str(originalVerts[vertIndex]["max_domain"])

    variableIndex += 2
    vertIndex += 1



  #next we add the constraints between the verticies themselves, start by
  #imposing the forward direction. If there is a constraint between two
  #edges A-->B of magnitude [0,5], then B cannot have in its interval any
  #times more than five seconds after the earliest start of A. Simillarly,
  #the upper bound of A cannot go later 0 before the lower bound of b
  for edge in originalEdges:
    integerFrom = int(edge["first_node"])
    integerTo = int(edge["second_node"])


    if(edge["min_duration"] == "inf"):
      #print "found inf"
      edgeMin = sys.float_info.max
    else:
      edgeMin = float(edge["min_duration"])

    if(edge["max_duration"] == "inf"):
      #print "found inf"
      edgeMax = sys.float_info.max
    else:
      edgeMax = float(edge["max_duration"])
    print("Edge from " + str(integerFrom) + " to " + str(integerTo) + " of ["+ str(edgeMin)+ ","+ str(edgeMax)+"]")
    print(doubleSTN[2+integerFrom*2].name + "-" + doubleSTN[2+integerTo*2 -1].name + "< 0")

    prob += doubleSTN[2+integerTo*2 - 1] - doubleSTN[2+integerFrom*2] <= edgeMax
    print "adding constraint that " + doubleSTN[2+integerTo*2 -1].name + " - " +  doubleSTN[2+integerFrom*2].name + " < " + str(edgeMax)

    prob += doubleSTN[2+integerFrom*2 -1] - doubleSTN[2+integerTo*2  ] <= -edgeMin
    print "adding constraint that " + doubleSTN[2+integerFrom*2 -1].name + " - " +  doubleSTN[2+integerTo*2].name + " < " + "-" +str(edgeMin)


  prob.writeLP("STN.lp")
  print "before solve"
  prob.solve()
  print "Status:", LpStatus[prob.status]


  # Each of the variables is printed with it's resolved optimum value
  for v in prob.variables():
    print v.name, "=", v.varValue

  total = 0
  for x in range(1, len(prob.variables())):
    if(x%2 == 1):
      total += float(prob.variables()[x+ 1].varValue) - float(prob.variables()[x].varValue)
      #print "adding " + prob.variables()[x].name + " - " + prob.variables()[x+1].name + " to the sum"

  print "Total flexibility is: " + str(total)


if __name__ == '__main__':
  main()
