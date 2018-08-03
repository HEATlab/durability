#!/usr/bin/python
import json, copy, os, re, sys, time
from json_stn_converter import convertJSON_old2new
from subprocess import Popen,PIPE
from multiprocessing import Pool
from itertools import product

## \file optimalDecouplingSearch_single.py
#
#  \brief Performs a brute force search to find the optimal decoupling for a
#  schedule over one makespan.
#
#  \details The algorithm first runs a search on an n-dimensional grid (where n
#  is the number of finite inter-agent constraints), which is constructed with
#  INITIAL_STEPS points in each dimension. This means it checks INITIAL_STEPS^n
#  points, so keep this reasonably small (like 20). The algorithm then will check around
#  the maximum point in a binary search, with the given number of iterations.
#  This can take a **very** long time on a large problem, so we suggest using a
#  computer with many cores, like Knuth. (see optimalDecouplingSearch_Knuth.py).
#  To run this program, use the command:
#  \code{.unparsed}
#  python optimalDecouplingSearch.py JSONFILE INITIAL_STEPS ITERATIONS [--intervalPlot PLOTFILE]
#  \endcode
#  This code utilizes a thread pool to decrease runtime. The size of the pool
#  should be changed depending on the machine. Remember that the simulator will
#  also run multiple threads, so be careful to not open a crazy number of
#  processes and kill your computer.


## \fn floydWarshallReduction(verts,edges)
#
#  \brief Takes vertices and edges of an STN loaded from a JSON file and runs
#  the Floyd Warshall algorithm to reduce the constraints to only include
#  the valid solution space.
#
#  \details Expects the vertices and edges to be in dictionaries indexed by
#  timepoints. See the main function for an example of how to convert input
#  JSON to the dictionary format
#
#  \returns the new verts and edges of the minimal STN
#
#  \warning This function in its current state modifies the input verts and
#  edges. Don't assume that the input variables are unchanged. They won't be.
def floydWarshallReduction(initialVerts,initialEdges):
  verts = initialVerts #copy.deepcopy(initialVerts)
  edges = initialEdges #copy.deepcopy(initialEdges)

  numVerts = len(verts) + 1
  stnMatrix = [[float('inf') for i in range(numVerts)] for j in range(numVerts)]

  for i,vert in verts.items():
    stnMatrix[0][i] = vert['max_domain']
    stnMatrix[i][0] = -vert['min_domain']

  for (i,j),edge in edges.items():
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

  for i,vert in verts.items():
    vert['max_domain'] = stnMatrix[0][i]
    vert['min_domain'] = -stnMatrix[i][0]

  for (i,j),edge in edges.items():
    edge['max_duration'] = stnMatrix[i][j]
    edge['min_duration'] = -stnMatrix[j][i]

  return verts,edges

## \fn getRobustness(verts,edges)
#  \brief Calls the Rust simulator to get the robustness of a schedule
#
#  \details Converts the STN to the new JSON format and then writes it to a file
#  with the process ID appended on the end to differentiate from other threads.
#  Calls the robustness simulator on this STN and returns the robustness.
#
#  \param verts
#  the vertices of the STN, in the dictionary format (not the original JSON)
#  \param edges
#  the edges of the STN, in the dictionary format (not the original JSON)
#
#  \note Assumes that there is a json folder in the directory that it is called
#  from
def getRobustness(verts, edges, numAgents):
  decouple_temp = os.path.dirname(os.path.abspath(__file__))+'/temp/decouple_temp_{0}.json'.format(os.getpid())
  jsonSTN = {'num_agents':numAgents,'nodes':verts.values(),'constraints':edges.values()}
  with open(decouple_temp, 'w+') as f:
    json.dump(jsonSTN, f, indent=2, separators=(',', ':'))

  # Find the robustness of the decoupling
  simulation = Popen(['../stpsimulator/target/release/simulator_stp',
                                '--samples', str(10000),
                                '--threads', '2',
                                '--sample_directory', '../stpsimulator/samples/',
                                decouple_temp],
                                stdout=PIPE, stderr=PIPE)
  simulation.wait()

  dataRegex = re.compile(ur'result:\n([0-9\.]+)')

  # extract the robustness from simulator output
  simData = simulation.stdout.read()
  match = dataRegex.search(simData)

  # simulator failed -- 0 robustness
  if match is None:
    return 0

  return float(match.group(1))

## \fn maxDecoupledRobustness(verts,edges,initialSteps,iterations,makespan)
#
#  \brief Runs the brute force search on a schedule for a given makespan
#         and returns (robustness, (pt, verts, edges))
def getDecoupledRobustness(tup):
  verts, edges, C_x, numAgents, pt = tup
  vertsCopy = copy.deepcopy(verts)
  edgesCopy = copy.deepcopy(edges)
  for index,(i,j) in enumerate(C_x):
    if (i,j) in edgesCopy:
      vertsCopy[i]['max_domain'] = pt[index]
      vertsCopy[j]['min_domain'] = pt[index] + edgesCopy[(i,j)]['min_duration']
    else:        
      vertsCopy[i]['max_domain'] = pt[index]
      vertsCopy[j]['min_domain'] = pt[index] - edgesCopy[(j,i)]['max_duration']
  robustness = getRobustness(vertsCopy,edgesCopy, numAgents)
  return robustness, (pt, vertsCopy, edgesCopy)


## \fn maximizeRobustness(makespan)
#
#  \brief Sets up the STN, runs maxDecoupledRobustness and writes the output to
#  a file
#
#  \details Also runs Floyd-Warshall before running maxDecoupledRobustness
#
def maximizeRobustness(verts, edges, initialSteps, iterations, C_x, numAgents, makespan = None):
  if makespan is not None:
    print 'makespan is {0}'.format(makespan)
    for node in verts:
      verts[node]['max_domain'] = makespan

  verts,edges = floydWarshallReduction(verts,edges)

  ptRanges = []
  stepsizes = []
  for i,j in C_x:
    minVal = verts[i]['min_domain']
    maxVal = verts[i]['max_domain']
    stepsize = (maxVal-minVal)/initialSteps
    stepsizes += [stepsize]
    ptRanges += [range(minVal,maxVal+stepsize,stepsize)]

  points = product(*ptRanges)
  tup = ((verts, edges, C_x, numAgents, pt) for pt in points)

  p = Pool(5)
  result = dict(p.map(getDecoupledRobustness,tup))

  best_robust = max(result.keys())
  best_result = result[best_robust]
  best_pt = best_result[0]
  if best_robust == 100.0:
    v, e = floydWarshallReduction(best_result[1], best_result[2]) 
    best_result = (best_pt, v, e)
    return 100, best_result

  #perform a binary search to improve the initial search
  for iter in range(iterations):
    for index,(i,j) in enumerate(C_x):
      stepsizes[index] = stepsizes[index]/2
      if stepsizes[index] <= 0:
        break
      ptRanges[index] = range( best_pt[index]-stepsizes[index],
                               best_pt[index]+2*stepsizes[index],
                               2*stepsizes[index])
    points = product(*ptRanges)
    tup = ((verts, edges, C_x, numAgents, pt) for pt in points)
    result = dict(p.map(getDecoupledRobustness,tup))
    robust = max(result.keys())
    best_result_pt = result[robust][0]

    if robust > best_robust:
        best_robust = robust
        best_pt = best_result_pt
        best_result = result[robust]
        if best_robust == 100:
          break

  v, e = floydWarshallReduction(best_result[1], best_result[2])  
  best_result = (best_pt, v, e)
  return best_robust, best_result

## \fn main
#  \brief Takes in command line arguements, initializes the STN, and uses a
#  thread pool to perform the robustness sweep.
#
#  \note The pool size should be changed depending on the number of cores on the
#  computer running this code.
def main(jsonFile, initialsteps, it, ms = None):
  usage = "python optimalDecouplingSearch.py JSONFILE "+\
          "INITIAL_STEPS ITERATIONS [MAKESPAN]"
  try:
    with open(jsonFile) as f:
      jsonSTN = json.load(f)

    initialSteps = int(initialsteps)
    iterations = int(it)

  except:
    print usage
    raise

  makespan = None
  if ms is not None:
    makespan = int(ms)

  numAgents = jsonSTN['num_agents']
  verts = dict((vert['node_id'],vert) for vert in jsonSTN["nodes"])
  edges = dict(((edge['first_node'],edge['second_node']),edge) for edge in jsonSTN["constraints"])

  agent = {}
  for i,vert in verts.items():
    agent[i] = vert['owner_id']

  C_x = []
  for (i,j), edge in edges.items():
    if agent[i] != agent[j]:
      if edge['min_duration'] != '-inf':
        C_x += [(i,j)]
      if edge['max_duration'] != 'inf':
        C_x += [(j,i)]

  robust, schedule = maximizeRobustness(verts, edges, initialSteps, iterations, C_x, numAgents, makespan)
  return robust, schedule




if __name__ == '__main__':
  main(*sys.argv[1:])
