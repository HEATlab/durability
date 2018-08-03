#!/usr/bin/python
import json, copy, os, re, sys, time
from json_stn_converter import convertJSON_old2new
from subprocess import Popen,PIPE
from multiprocessing import Pool
from itertools import product

## \file optimalDecouplingSearch.py
#
#  \brief Performs a brute force search to find the optimal decoupling for a
#  schedule over multiple makespans.
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
#  python optimalDecouplingSearch.py JSONFILE MAX_TIME MIN_TIME TIMESTEP INITIAL_STEPS ITERATIONS [--intervalPlot PLOTFILE]
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
def getRobustness(verts,edges):
  global numAgents
  decouple_temp = 'json/decouple_temp_{0}.json'.format(os.getpid())

  jsonSTN = {'num_agents':numAgents,'nodes':verts.values(),'constraints':edges.values()}
  with open(decouple_temp, 'w+') as f:
    json.dump(jsonSTN, f, indent=2, separators=(',', ':'))

  # Find the robustness of the decoupling
  simulation = Popen(['../stpsimulator/target/release/simulator_stp',
                                '--samples', str(10000),
                                '--threads', '4',
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
def maxDecoupledRobustness(verts,edges,initialSteps,iterations,makespan):
  global C_x

  ptRanges = []
  stepsizes = []
  for i,j in C_x:
    minVal = verts[i]['min_domain']
    maxVal = verts[i]['max_domain']
    stepsize = (maxVal-minVal)/initialSteps
    stepsizes += [stepsize]
    ptRanges += [range(minVal,maxVal+stepsize,stepsize)]

  maxRobustness = -1
  best_pt = None
  for pt in product(*ptRanges):
    for index,(i,j) in enumerate(C_x):
      if (i,j) in edges:
        verts[i]['max_domain'] = pt[index]
        verts[j]['min_domain'] = pt[index] + edges[(i,j)]['min_duration']
      else:        
        verts[i]['max_domain'] = pt[index]
        verts[j]['min_domain'] = pt[index] - edges[(j,i)]['max_duration']
    robustness = getRobustness(verts,edges)
    if robustness > maxRobustness:
      maxRobustness = robustness
      best_pt = pt
    if maxRobustness == 100:
          return 100

  #perform a binary search to improve the initial search
  for iter in range(iterations):
    for index,(i,j) in enumerate(C_x):
      stepsizes[index] = stepsizes[index]/2
      if stepsizes[index] <= 0:
        break
      ptRanges[index] = range( best_pt[index]-stepsizes[index],
                               best_pt[index]+2*stepsizes[index],
                               2*stepsizes[index])
    for pt in product(*ptRanges):
      for index,(i,j) in enumerate(C_x):
        if (i,j) in edges:
          verts[i]['max_domain'] = pt[index]
          verts[j]['min_domain'] = pt[index] + edges[(i,j)]['min_duration']
        else:        
          verts[i]['max_domain'] = pt[index]
          verts[j]['min_domain'] = pt[index] - edges[(j,i)]['max_duration']
      robustness = getRobustness(verts,edges)
      if robustness > maxRobustness:
        maxRobustness = robustness
        best_pt = pt
      if maxRobustness == 100:
            return 100

  return maxRobustness

## \fn maximizeRobustness(makespan)
#
#  \brief Sets up the STN, runs maxDecoupledRobustness and writes the output to
#  a file
#
#  \details Writes to the file `./data/optimalSweep_#.dat` where # is the
#  process ID. This allows multiple threads to run and not attempt to write to
#  the same file. Also runs Floyd-Warshall before running maxDecoupledRobustness
#
#  \note Assumes that there is a data folder in the directory that it is called
#  from
def maximizeRobustness(makespan):
  global verts,edges
  global initialSteps,iterations
  global C_x, C_c
  global intervalPlot, intervalPlotFile

  start_time = time.time()
  outfile = './data/optimalSweep_{0}.dat'.format(os.getpid())
  local_verts = copy.deepcopy(verts)
  local_edges = copy.deepcopy(edges)

  for node in local_verts:
    local_verts[node]['max_domain'] = makespan

  newVerts,newEdges = floydWarshallReduction(local_verts,local_edges)
  robustness = maxDecoupledRobustness(newVerts,newEdges,initialSteps,iterations,makespan)

  with open(outfile,"a+") as f:
    f.write("{0} {1}\n".format(makespan,robustness))
  if intervalPlot:
    newVerts,newEdges = floydWarshallReduction(newVerts,newEdges)
    with open(intervalPlotFile,"a+") as f:
      f.write("{0} ".format(makespan))
      for i,j in C_c:
        f.write("{0},{1}, {2} ".format(j,i,newEdges[(i,j)]["min_duration"]))
        f.write("{0},{1}, {2} ".format(i,j,newEdges[(i,j)]["max_duration"]))
      for i,j in C_x:
        f.write("{0},{1}, {2} ".format(i,0,newVerts[i]["min_domain"]))
        f.write("{0},{1}, {2} ".format(0,i,newVerts[i]["max_domain"]))
      f.write("\n")

  end_time = time.time()
  print makespan,robustness
  print 'That took %0.3f seconds' % (end_time - start_time)

## \fn main
#  \brief Takes in command line arguements, initializes the STN, and uses a
#  thread pool to perform the robustness sweep.
#
#  \note The pool size should be changed depending on the number of cores on the
#  computer running this code.
def main():
  global verts,edges
  global initialSteps,iterations,outfile
  global C_x, C_c
  global intervalPlot, intervalPlotFile
  global numAgents

  usage = "python optimalDecouplingSearch.py JSONFILE MAX_TIME MIN_TIME "+\
          "TIMESTEP INITIAL_STEPS ITERATIONS [--intervalPlot PLOTFILE]"
  try:
    with open(sys.argv[1]) as f:
      jsonSTN = json.load(f)

    min_time = int(sys.argv[2])
    max_time = int(sys.argv[3])
    timestep = int(sys.argv[4])
    initialSteps = int(sys.argv[5])
    iterations = int(sys.argv[6])

    intervalPlot = False
    intervalPlotFile = None
    if len(sys.argv) > 7:
      if sys.argv[7] == "--intervalPlot":
        intervalPlotFile = sys.argv[8]
        try:
          with open(intervalPlotFile,"w+") as f:
            intervalPlot = True
        except:
          print "invalid file name: {0}".format(intervalPlotFile)
          raise
  except:
    print usage
    raise

  numAgents = jsonSTN['num_agents']
  verts = dict((vert['node_id'],vert) for vert in jsonSTN["nodes"])
  edges = dict(((edge['first_node'],edge['second_node']),edge) for edge in jsonSTN["constraints"])

  agent = {}
  for i,vert in verts.items():
    agent[i] = vert['owner_id']

  C_x = []
  C_c = []
  for (i,j), edge in edges.items():
    if agent[i] != agent[j]:
      if edge['min_duration'] != '-inf':
        C_x += [(i,j)]
      if edge['max_duration'] != 'inf':
        C_x += [(j,i)]
    if "distribution" in edge:
      C_c += [(i,j)]

  p = Pool(20)
  makespan_range = range(min_time, max_time + timestep, timestep)
  p.map(maximizeRobustness,makespan_range)


if __name__ == '__main__':
  main()
