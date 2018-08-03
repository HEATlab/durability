#!/usr/bin/python
import sys
from loadInverseCDFs import invCDF
from stntools import *

## \file guaranteedRobustness.py
#  \brief Calculates a strict lower bound for the robustness of a schedule,
#         regardless of dispatch strategy
#  \details To run, use the command `python guaranteedRobustness.py JSON_FILE`

## \fn getCDF(distribution,val)
#  \brief Gets the CDF of the given distribuation evaluated at val
#  \details Currently runs a binary search and uses the inverse CDF. This should
#  probably be switched out for a direct calculation of some sort...
def getCDF(distribution,val):
    # dictionary of alphas for binary search
    alphas = { i : str(round(i/1000.0,3)) for i in range(1001)}
    # bounds for binary search
    lower  = 0
    upper = 1001
    # Run binary search on alpha
    while upper - lower > 1:
      alpha = alphas[(upper + lower) // 2]
      # try lower alpha
      if invCDF[distribution][alpha] > val/1000.0:
        upper = (upper + lower) // 2
      # try higher alpha
      else:
        lower = (upper + lower) // 2

      # finished our search
      if upper - lower <= 1:
        return float(alphas[(upper + lower) // 2])

    #should never happen
    return None

## \fn guaranteedRobustness(STN)
#  \brief Calculates a lower bound on the robustness of an STN
#  \returns Robustness in the range 0-100 if we can guarantee any robustness,
#  and None otherwise
def guaranteedRobustness(STN, debug = False):
    if not STN.floydWarshall():
      if debug:
        print STN
        print "Given STN was inconsistent"
      return None

    if debug:
      print STN

    robustness = 100
    for edge in STN.contingentEdges.values():
      i = edge.i
      j = edge.j
      distribution = edge.distribution
      worst_case_max_dur = STN.getEdgeWeight(0,j) - STN.getEdgeWeight(0,i)
      worst_case_min_dur = -( STN.getEdgeWeight(j,0) - STN.getEdgeWeight(i,0) )

      if debug:
        print "Worst case max for {} => {}: {}".format(i,j,worst_case_max_dur)
        print "Worst case min for {} => {}: {}".format(i,j,worst_case_min_dur)

      #This should never happen but better safe than sorry
      if worst_case_max_dur < 0 or worst_case_min_dur < 0:
        if debug:
          print "Edge {} => {} is inconsistent".format(i,j)
        return None

      if worst_case_min_dur > worst_case_max_dur:
        if debug:
          print "Cannot guarantee any robustness. Failed on edge {} => {}".format(i,j)
        return None

      CDF_max = getCDF(distribution,worst_case_max_dur)
      CDF_min = getCDF(distribution,worst_case_min_dur)

      if debug:
        print "\t{} -- {}".format(CDF_min,CDF_max)

      robustness *= (CDF_max-CDF_min)
    return robustness

## \fn main
#  \brief Calculates a lower bound for robustness and prints the result
def main():
  usage = "USAGE: {} JSON_FILE".format(sys.argv[0])
  if len(sys.argv) != 2:
    print usage
    sys.exit(1)

  try:
    STN = loadSTNfromJSONfile(sys.argv[1])['stn']
  except:
    print "Could not load STN from JSON file: {}".format(sys.argv[1])
    print usage
    raise

  print STN

  robustness = guaranteedRobustness(STN, debug = True)
  if robustness is not None:
    print "\nWe can guarantee {}% robustness".format(robustness)
  else:
    print "\nWe cannot guarantee any robustness"

if __name__ == "__main__":
  main()
