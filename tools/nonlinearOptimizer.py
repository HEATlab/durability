from guaranteedRobustness import guaranteedRobustness
import sys, os, re
from stntools import *
from itertools import permutations
import scipy.optimize
from subprocess import Popen, PIPE
from schedulase import schedulase
from Queue import Queue

## \file
#  \brief Usages a non-linear optimizer to maximize the guaranteed robustness of
#  a schedule
#  \details Currently does not optimize because it would require a better
#  starting point (and possibly algorithm)

def toOptimize(constraint_values, mySTN):
  mySTN = mySTN.copy()
  constraint_values = iter(constraint_values)

  for edge in mySTN.getAllEdges():
    edge.Cij = constraint_values.next()
    edge.Cji = constraint_values.next()

  robustness = guaranteedRobustness(mySTN)
  if robustness == None:
    return 1
  else:
    return -robustness

def updateSTNwithResult(constraint_values,mySTN):
  constraint_values = iter(constraint_values)

  for edge in mySTN.getAllEdges():
    edge.Cij = constraint_values.next()
    edge.Cji = constraint_values.next()

  return mySTN

def getConstraints(mySTN):
  constraints = []
  bounds = {}
  indexes = {}
  for (index,edge) in enumerate(mySTN.getAllEdges()):
    i = edge.i
    j = edge.j

    bounds[(i,j)] = edge.Cij
    bounds[(j,i)] = edge.Cji

    indexes[(i,j)] = 2*index
    indexes[(j,i)] = 2*index+1


  for i,j in bounds:
    #valid bounds
    constraints += [ {'type':'ineq',
                      'fun':lambda x,*args:x[indexes[(i,j)]]+x[indexes[(j,i)] ] } ]
    #Initial constraints
    constraints += [ {'type':'ineq',
                      'fun':lambda x,*args: bounds[(i,j)]-x[indexes[(i,j)]] }, \
                     {'type':'ineq',
                      'fun':lambda x,*args: bounds[(j,i)]-x[indexes[(j,i)]] } ]


  for t in mySTN.tris:
    for i,j,k in permutations( (t.i,t.j,t.k) ):
      b_ij = indexes[(i,j)]
      b_ik = indexes[(i,k)]
      b_kj = indexes[(k,j)]

      constraints += [{'type':'ineq',
                        'fun':lambda x,*args: x[b_ik]+x[b_kj]-x[b_ij] } ]

  return constraints

def getBounds(mySTN):
  bounds = []
  for edge in mySTN.getAllEdges():
    bounds += [(None,edge.Cij),(None,edge.Cji)]
  return bounds

def getInitialGuess(mySTN):
  mySTN = mySTN.copy()
  x = []
  updateQueue = Queue()
  triangleQueue = Queue()

  for edge in mySTN.getAllEdges():
    x += [-edge.Cji,edge.Cji]
    updateQueue.put({"i":edge.i,"j":edge.j,"Wij":-edge.Cji,"Wji":-edge.Cji})
    DItriSTP(mySTN,updateQueue,triangleQueue)
  return x

def main():
  usage = 'USAGE: {} JSON_FILE'.format(sys.argv[0])
  try:
    stn_dict = loadSTNfromJSONfile(sys.argv[1])
    mySTN = stn_dict['stn']
    numAgents = stn_dict['agent_count']
  except:
    print 'Could not load mySTN from JSON file: {}'.format(sys.argv[1])
    print usage
    raise
  mySTN.minimize()
  #output = schedulase(mySTN.forJSON(),0.45)
  #mySTN,numAgents = loadSTNfromJSON(output['STN'])
  print mySTN

  constraints = getConstraints(mySTN)
  bounds = getBounds(mySTN)
  init_guess = getInitialGuess(mySTN)

  result = scipy.optimize.minimize( toOptimize,
                                      init_guess,
                                      method = 'SLSQP',
                                      bounds = bounds,
                                      constraints = constraints,
                                      args = (mySTN,),
                                      options = {'disp':True})


  if not result.success:
    print "Simulation failed"
    print result.message
    print "\nInitial\tFinal"
    for i in range(len(init_guess)):
      print "{}\t{:6.1f}".format(init_guess[i],result.x[i])
    sys.exit(1)

  updateSTNwithResult(result.x,mySTN)
  print mySTN

  print 'guaranteedRobustness:{}'.format(guaranteedRobustness(mySTN))
  # The simulator requires file input
  outFile = os.path.dirname(os.path.abspath(__file__))+'/json/decouple_temp.json'
  with open(outFile,'w') as f:
    json.dump(mySTN.forJSON(), f, indent = 2, separators=(',',':'))

  # Find the robustness of the decoupling
  simulation = Popen(['../stpsimulator/target/release/simulator_stp',
                                '--samples', str(10000),
                                '--sample_directory', '../stpsimulator/samples/',
                                outFile],
                                env=os.environ, stdout=PIPE)
  simulation.wait()

  # extract the robustness from simulator output
  simData = simulation.stdout.read()
  dataRegex = re.compile(ur'result:\n([0-9\.]+)')
  match = dataRegex.search(simData)

  if match is None:
    print 'NO REGEX MATCH. This was the output'
    print simData
    sys.exit(1)

  print 'Robustness: {:6.2f}'.format(float(match.group(1)))


if __name__ == '__main__':
  main()
