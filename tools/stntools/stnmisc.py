#!/usr/bin/env python

## \fn replaceReceivedTimepoints(STN_1,STN_2)
#  \brief modifies the zero timepoint bounds of the received timepoints
#   of STN_2 so that they are the same as STN_1
def replaceReceivedTimepoints(STN_1, STN_2):
  for timepoint in STN_1.receivedTimepoints:
    STN_2.verts[timepoint] = STN_1.verts[timepoint]
    STN_2.edges[(0,timepoint)] = STN_1.edges[(0,timepoint)]
  return STN_2