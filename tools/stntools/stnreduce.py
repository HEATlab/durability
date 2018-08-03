import sys
import stn
import math

##
# \fn stnreduce(curSTN, agent, origSTN, debug = False)
# \brief Runs one iteration of triangulation on an STN
# \note This function must be run in a loop to triangulate an STN
# \details To run this function, you need an STN and a copy of that
#     STN. Both will be changed (The original is triangulated, and
#     all vertices are eventually removed from the copy). Loop while
#     the copy is not empty, and for every time through the loop run
#     reduce once on each agent. Here is example code for the loop,
#     with a simple timer to ensure all agents have approximately
#     the same number of triangles:
#
#     Also note, really, this is rubbish code. Legacy... :(
#     FIXME: This file is interdependent on stn, and stn is
#     interdependent on this file.
#     That's all rubbish. Fix that.
#
# \code{.py}
# #######################
# # Triangulate the STN #
# #######################
#
# stnCopy = mySTN.copy()
# agentWait = []
# # Set up the wait timer for agent load balancing
# for a in mySTN.agents:
#   agentWait.append(0)
#
# # Perform the reduction
# while len(stnCopy.verts) > 1:
#   for a in range(len(mySTN.agents)):
#     if agentWait[a] == 0:
#       created = stnreduce(stnCopy, mySTN.agents[a], mySTN)
#       agentWait[a] = created
#     else:
#       agentWait[a] -= 1
# \endcode
#
#  \param curSTN The current STN during triangulation. This STN has
#      nodes removed from previous iterations. Nodes are removed after
#      the function checks for triangles with that node.
#      Note: Keep in mind that curSTN and origSTN is passed by reference.      
#  \param agent The agent whose vertices we are using to triangulate.
#      All triangles created while running will belong to this agent
#  \param origSTN The STN that we are triangluating. This STN is
#      changed during triangulation.
def stnreduce(curSTN, agent, origSTN, debug = False):
  verts = curSTN.getAgentVerts(agent)
  # Find the vertex with minimum degree

  # ID of the vertex with the lowest degree.
  minVert = None
  # Degree of the minimum vertex
  minDegree = sys.maxint

  for vert in verts:
    vertDegree = curSTN.getDegree(vert.nodeID)
    if  vertDegree < minDegree:
      minVert = vert.nodeID
      minDegree = vertDegree

  # Now that we have the ID of the node to remove we must find all the triangles
  # it is going to create.
  adjacentVerts = curSTN.getAdjacent(minVert)
  triCount = 0
  # for each pair of vertices i,j adjacent to the min degree vertex, 
  #          check if an edge exists between i and j, because that implies a 
  #          triangle in the original graph.
  for indx, i in enumerate(adjacentVerts):
    for j in adjacentVerts[indx+1:]:
      triCount += 1
      if curSTN.edgeExists(i, j):
        # If an edge (i.e. constraint) exists between i and j, 
        #    add the triangle that is formed to the original STN
        origSTN.addTriangle(agent, i, j, minVert)
      else:
        Eik = curSTN.getEdge(i, minVert)
        Ejk = curSTN.getEdge(j, minVert)

        # If an edge (i.e. constraint) does not exist between i and j, 
        #    infer an artificial constraint using 
        #    the edges between i and minvert, minvert and j
        maxVal = Eik.getWeight(i, minVert) + Ejk.getWeight(minVert, j)
        minVal = Ejk.getWeight(j, minVert) + Eik.getWeight(minVert, i)

        curSTN.addEdge(i, j,\
                        -minVal,\
                        maxVal)

        if math.isnan(Eik.getWeight(i, minVert) + Ejk.getWeight(minVert, j)) or math.isnan(Ejk.getWeight(j, minVert) + Eik.getWeight(minVert, i)):
          print Eik.forJSON()
          print Ejk.forJSON()
          print "%d->%d:%f" %(i, minVert, Eik.getWeight(i, minVert))
          print "%d->%d:%f" %(minVert, j, Ejk.getWeight(minVert, j))
          print "%d->%d:%f" %(j, minVert, Ejk.getWeight(j, minVert))
          print "%d->%d:%f" %(minVert, i, Eik.getWeight(minVert, i))
          print "%d->%d->%d:%f"%(i, minVert, j, Eik.getWeight(i, minVert) + Ejk.getWeight(minVert, j))
          print "%d->%d->%d:%f"%(j, minVert, i, Ejk.getWeight(j, minVert) + Eik.getWeight(minVert, i))
          sys.exit(1)

        if debug:
          print "ADDING: %d->%d: [%f, %f]" % (i, j, minVal, maxVal)

        # add this inferred constraint to the original STN and the triangle 
        # that is formed.
        origSTN.addEdge(i, j,\
                        -minVal,\
                        maxVal, fake=True)
        e = origSTN.getEdge(i, j)
        if debug:
          print "[%f,%f]" % (e.getWeightMin(), e.getWeightMax())
        origSTN.addTriangle(agent, i, j, minVert)
  
  # All of the triangles minvert could be part of are now added to origSTN.
  # Remove minvert from curSTN, which keeps track of the progress of triangulation
  curSTN.removeVertex(minVert)
  return triCount
