import json, sys
from pulp import *

## \file simpleLinearDecoupling.py
#
#  \brief Solves the naive linear program for decoupling a schedule
#
#  \details The program first runs Floyd-Warshall on the given STN, then sets up
#  and solves the naive LP
#
#  \note You'll need pulp to be able to run this. On a linux machine run
#  `sudo pip install pulp` or `sudo easy_install -U pulp`.
#  THEN RUN `sudo pulptest`, otherwise it won't work.

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

## \fn addConstraint(constraint,problem)
#  \brief Adds an LP constraint to the given LP
def addConstraint(constraint,problem):
    problem += constraint
    print 'adding constraint', constraint

## \fn main
#  \brief Sets up and solves the LP
def main():
  data = []

  with open(sys.argv[1]) as f:
    data = json.load(f)

  originalVerts = data['nodes']
  originalEdges = data['constraints']

  originalVerts,originalEdges = floydWarshallReduction(originalVerts,originalEdges)
  print json.dumps(originalVerts, sort_keys=True,
                    indent=4, separators=(',', ': '))
  print json.dumps(originalEdges, sort_keys=True,
                    indent=4, separators=(',', ': '))


  doubleSTN = {0:{'+':LpVariable('t0Up'),'-':LpVariable('t0Lo')}}

  initialTimes = {}
  for vert in originalVerts:
    initialTimes[vert['node_id']] = {'+':vert['max_domain'],'-':vert['min_domain']}
    doubleSTN[vert['node_id']] = {'+':LpVariable('t'+str(vert['node_id'])+'Up'),
                                 '-':LpVariable('t'+str(vert['node_id'])+'Lo'),
                                 'agent': vert['owner_id']}
  numVerts = len(doubleSTN)

  prob = LpProblem('STN Flex Finder',LpMinimize)

  for i in range(numVerts):
    addConstraint(doubleSTN[i]['+'] >= doubleSTN[i]['-'], prob)

  prob += doubleSTN[0]['+'] == 0
  prob += doubleSTN[0]['-'] == 0

  # add the constraints imposed by edges from the origin
  # We know that each event must start within a closed range of times, so we add
  # the upper and lower constraints for each edge
  for vert in originalVerts:
    i = vert['node_id']
    addConstraint( doubleSTN[i]['-'] >= doubleSTN[0]['+'] + vert['min_domain'], prob)
    addConstraint( doubleSTN[i]['+'] <= doubleSTN[0]['-'] + vert['max_domain'], prob)

  #next we add the constraints between the verticies themselves, start by
  #imposing the forward direction. If there is a constraint between two
  #edges A-->B of magnitude [0,5], then B cannot have in its interval any
  #times more than five seconds after the earliest start of A. Simillarly,
  #the upper bound of A cannot go later 0 before the lower bound of b

  epsilons = []

  for edge in originalEdges:
    i = int(edge['first_node'])
    j = int(edge['second_node'])

    if(edge['min_duration'] == 'inf'):
      edgeMin = sys.float_info.max
    else:
      edgeMin = float(edge['min_duration'])

    if(edge['max_duration'] == 'inf'):
      edgeMax = sys.float_info.max
    else:
      edgeMax = float(edge['max_duration'])

    #Interagent constraint
    if doubleSTN[i]['agent'] != doubleSTN[j]['agent']:
      addConstraint( doubleSTN[j]['+'] - doubleSTN[i]['-'] <= edgeMax, prob)
      addConstraint( doubleSTN[i]['+'] - doubleSTN[j]['-'] <= -edgeMin, prob)
      continue


    #Contingent constraint
    elif 'distribution' in edge:
      e_i_up = LpVariable('e_'+str(i)+'_Up')
      e_i_lo = LpVariable('e_'+str(i)+'_Lo')
      e_j_up = LpVariable('e_'+str(j)+'_Up')
      e_j_lo = LpVariable('e_'+str(j)+'_Lo')
      addConstraint( doubleSTN[j]['+'] == initialTimes[j]['+'] - e_j_up, prob)
      addConstraint( doubleSTN[j]['-'] == initialTimes[j]['-'] + e_j_lo, prob)
      addConstraint( doubleSTN[i]['+'] == initialTimes[i]['+'] - e_i_up, prob)
      addConstraint( doubleSTN[i]['-'] == initialTimes[i]['-'] + e_i_lo, prob)
      epsilons += [e_i_up, e_i_lo, e_j_up, e_j_lo]

    # Arc consistency
    addConstraint( doubleSTN[j]['+'] - doubleSTN[i]['+'] <= edgeMax, prob)
    addConstraint( doubleSTN[i]['+'] - doubleSTN[j]['+'] <= -edgeMin, prob)
    addConstraint( doubleSTN[j]['-'] - doubleSTN[i]['-'] <= edgeMax, prob)
    addConstraint( doubleSTN[i]['-'] - doubleSTN[j]['-'] <= -edgeMin, prob)


  max_e = LpVariable('max_e')
  for e1 in epsilons:
    addConstraint(e1 <= max_e,prob)

  # ##
  # Generate the objective function.
  #   We are currently using (SUM e_ij)
  # ##
  intervalSum = len(epsilons)*max_e
  for epsilon in epsilons:
    addConstraint(epsilon >= 0, prob)
    intervalSum += epsilon

  prob += intervalSum, 'Minimize change from shortest path STN while decoupling'

  prob.writeLP('STN.lp')
  print 'before solve'
  prob.solve()
  print 'Status:', LpStatus[prob.status]


  # Each of the variables is printed with it's resolved optimum value
  for v in prob.variables():
    print v.name, '=', v.varValue


if __name__ == '__main__':
  main()
