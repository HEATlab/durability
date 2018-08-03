import json, sys, os, re
import pulp
import argparse
import stntools

## \file determine_STN_qualitites.py
#
#  \brief This file is intended to make life easier when analyzing an STN.
#   It takes in an STN and then lists properties of the STN that may be useful
#   to know such as if the STN is minimal, an interval schedule, etc...
#
#  \details To run this file on an STN, use the command:
#  \code{.unparsed}
#  python SREA.py JSONFILE ALPHA
#  \endcode
#
#
#  \note You'll need pulp to be able to run this. On a linux machine run
#  `sudo pip install pulp` or `sudo easy_install -U pulp`.
#  THEN RUN `sudo pulptest`, otherwise it won't work.

## \fn isIntervalSchedule(STN)
#  \brief determines if the given STN is an interval schedule
#  \returns A boolean value indicating if the STN is an interval schedule
def isIntervalSchedule(STN):
  bounds = {}

  # stores the upper and lower bounds of each timepoint in a dictionary
  for i in STN.verts:
    bounds[(i,'+')] = STN.getEdgeWeight(0,i)
    bounds[(i,'-')] = STN.getEdgeWeight(i,0)

  # loops through all the edges while checking to make sure they
  for (i,j) in STN.edges:
    if (i,j) in STN.contingentEdges:
      raise ValueError("this STN has contingent edges. Please convert these to intervals")

    else:
      if bounds[(j,'+')] - bounds[(i,'-')] > STN.getEdgeWeight(i,j):
        return False

      if -bounds[(j,'-')] - bounds[(i,'+')] < -STN.getEdgeWeight(j,i):
        return False

  return True


## \fn isControllableSchedule(STN, risk)
#  \brief determines if the contingent edges for a given STN are controllable.  It
#   uses the risk budget to create intervals from the distributions.  These intervals are
#   used to test for controllability.
#  \returns A boolean value indicating if the STN is a controllable schedule
def isControllableSchedule(STN, risk):
  bounds = {}

  for i in STN.verts:
    bounds[(i,'+')] = STN.getEdgeWeight(0,i)
    bounds[(i,'-')] = STN.getEdgeWeight(i,0)

  for (i,j) in STN.edges:
    maxEdge = STN.getEdgeWeight(i,j)
    minEdge = -STN.getEdgeWeight(j,i)
    i_hi = bounds[(i,'+')]
    i_lo = bounds[(i,'-')]
    j_hi = bounds[(j,'+')]
    j_lo = bounds[(j,'-')]

    

    if (i,j) in STN.contingentEdges:
      if (j_lo > i_lo + minEdge or j_hi < i_lo + minEdge) and (i_lo > j_lo - minEdge or i_hi < j_lo - minEdge):
        return False

      if (i_lo > i_hi + maxEdge or j_hi < i_hi + maxEdge) and (i_lo > j_hi - maxEdge or i_hi < j_hi - maxEdge):
        return False

  return True



def main():
# handle command line input

  parser = argparse.ArgumentParser()

  parser.add_argument('jsonFile', metavar='JSON', type=str,
          help='STN json file')

  options = parser.parse_args()

  if options.jsonFile == None:
    parser.error("Incorrect number of arguments")

  try:
    STN = stntools.loadSTNfromJSONfile(options.jsonFile)['stn']
  except IOError as e:
    sys.stderr.write("I/O error({0}): {1}\n".format(e.errno, e.strerror))
    print(sys.exc_info()[0])
    sys.exit(1)
  except TypeError as e:
    sys.stderr.write("Type error: {0}\n".format(e))
    print(traceback.format_exc())
    sys.exit(1)

  interval = isIntervalSchedule(STN)
  controllable = isControllableSchedule(STN)

  print "interval schedule: ", interval
  print "controllable schedule: ", controllable


if __name__ == '__main__':
  main()
