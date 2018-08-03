from bruteForceStatic import *
import sys

def main():
  jsonfile = sys.argv[1]
  initialSteps = int(sys.argv[2])
  iterations = int(sys.argv[3])
  lower = int(sys.argv[4])
  upper = int(sys.argv[5])
  step = int(sys.argv[6])
  makespans = range(lower,upper+step,step)

  STN = loadSTNfromJSONfile(jsonfile)['stn']

  for makespan in makespans:
    print "Running makespan {}".format(makespan)
    STN.setMakespan(makespan)
    robustness,alphas = binarySearch(STN,initialSteps,iterations)
    print >> sys.stderr, makespan, robustness

if __name__ == "__main__":
  main()
