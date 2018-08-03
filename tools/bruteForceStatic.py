from loadInverseCDFs import invCDF
from stntools import *
from itertools import product
from schedulase_stntools import *
from multiprocessing import Pool

ALPHA_RESOLUTION = 0.001

def formatAlpha(alpha):
  return str(round(float(alpha),3))

def updateSTN(STN,alphas):
  if alphas == None:
    return

  i = 0
  for edge in STN.contingentEdges.values():
    low_alpha = formatAlpha(alphas[i])
    high_alpha = formatAlpha(alphas[i+1])
    edge.Cji = -1000*invCDF[edge.distribution][low_alpha]
    edge.Cij = 1000*invCDF[edge.distribution][high_alpha]
    i += 2


def binarySearch(STN, initialSteps, iterations):
  if not STN.possiblyStronglyControllable:
    return 0.0,None

  STN = STN.copy()
  #start with uniform alpha alpha_lists
  stepsize = 1.0/(initialSteps-1)
  alpha_lists = [[round(float(i)*stepsize,3) for i in range(initialSteps)] for _ in range(2*len(STN.contingentEdges))]

  maxRobustness = -1
  best_pt = None
  for alphas in product(*alpha_lists):

    good_alphas = True
    for i in range(0,len(alphas),2):
      if alphas[i]>alphas[i+1]:
        good_alphas = False
        break
    if not good_alphas:
      continue

    updateSTN(STN,alphas)
    #print STN
    if STN.isStronglyControllable():
      robustness = 100
      for i in range(0,len(alphas),2):
        robustness *= alphas[i+1]-alphas[i]
      if robustness > maxRobustness:
        maxRobustness = robustness
        best_alphas = alphas
      if maxRobustness == 100:
            return 100,alphas
  #print "Finished first search"

  #perform a binary search to improve the initial search
  for i in range(iterations):
    #print "Starting secondary search #{}".format(i+1)
    stepsize = int(stepsize/(2*ALPHA_RESOLUTION))*ALPHA_RESOLUTION
    for i in range(len(alpha_lists)):
      alpha_lists[i] = [max(best_alphas[i]-stepsize,0.0),best_alphas[i],min(best_alphas[i]+stepsize,1.0)]

    for alphas in product(*alpha_lists):
      updateSTN(STN,alphas)
      if STN.isStronglyControllable():
        robustness = 100
        for i in range(0,len(alphas),2):
          robustness *= alphas[i+1]-alphas[i]
        if robustness > maxRobustness:
          maxRobustness = robustness
          best_alphas = alphas
        if maxRobustness == 100:
              return 100,alphas

  return maxRobustness,best_alphas

def main():
  usage = "python optimalDecouplingSearch.py INITIAL_STEPS ITERATIONS THREADS JSONFILE(S)"
  try:
    initialSteps = int(sys.argv[1])
    iterations = int(sys.argv[2])
    threads = int(sys.argv[3])
    STNs = map(lambda jsonfile: loadSTNfromJSONfile(jsonfile)['stn'],sys.argv[4:])
  except:
    print usage
    raise

  p = Pool(threads)
  p.map(compareToLP, [(STN,initialSteps,iterations) for STN in STNs] )
    

def compareToLP(args):
  STN,initialSteps,iterations = args
  print STN
  robustness,alphas = binarySearch(STN,initialSteps,iterations)

  LPoutput = schedulaseBinarySearch(STN.copy())
  if LPoutput != None:
    alpha,LP_STN = LPoutput
    LProbustness = getRobustness(LP_STN)
  else:
    LProbustness = 0.0

  updateSTN(STN,alphas)
  print STN
  print STN.isStronglyControllable(debug=True)
  print robustness, LProbustness

if __name__ == '__main__':
  main()
