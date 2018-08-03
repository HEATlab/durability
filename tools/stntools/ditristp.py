from Queue import Queue
import itertools

import stn

##
# \fn DItriSTP(STN, updateQueue, triangleQueue)
# \brief runs DItriSTP on an STN
#
# @param STN The Simple Temporal Network to check for
# @param updateQueue
# @param triangleQueue
# @return Returns a boolean indicating whether or not the STN is consistent.
def DItriSTP(STN, updateQueue, triangleQueue):
  while not updateQueue.empty() or not triangleQueue.empty():
    while not updateQueue.empty():
      update = updateQueue.get()
      e = STN.getEdge(update['i'], update['j'])
      #print e.getWeightMax()
      forward = STN.updateEdge(update['i'], update['j'], update['Wij'])
      reverse = STN.updateEdge(update['j'], update['i'], update['Wji'])

      if e.getWeightMin() > e.getWeightMax():
        # If we somehow set the minimum weight to be larger than the
        # maximum, we have an inconsistent STN.
        return False

      if forward or reverse:
        adjTri = STN.getAdjacentTriangles(update['i'], update['j'])
        for tri in adjTri:
          triangleQueue.put(tri)

    # Early termination if we didn't make any updates
    if triangleQueue.empty():
      break

    curTri = triangleQueue.get()
    updatedEdges = []
    # running through this for-loop minimizes the weight on edges (i,j), (j,k), (k,i)
    #   (both ways) by inferrring the total weight of the edge by using a path 
    #   over the other two edges in the triangle.
    for order in itertools.permutations(curTri.nodes()):
      i = order[0]
      j = order[1]
      k = order[2]

      weightSum = STN.getEdgeWeight(i, k) + \
                  STN.getEdgeWeight(k, j)
      if STN.updateEdge(i, j, weightSum):
        updatedEdges.append([i, j])

    for edge in updatedEdges:
      # Get the triangles adjacent to vertices i, j that were updated previously 
      adjTri = STN.getAdjacentTriangles(edge[0], edge[1])
      for tri in adjTri:
        triangleQueue.put(tri)
      # get edge between i, j
      sendEdge = STN.getEdge(edge[0], edge[1])
      updateQueue.put({'i':sendEdge.i, 'j':sendEdge.j, \
                       'Wij':sendEdge.Cij, 'Wji':sendEdge.Cji})

  return True


## \fn DItriSTP_alpha(STN, updateQueue, triangleQueue, alpha)
#  \brief runs a modified version of DItriSTP on an STN that guarantees 1-alpha
#  chance of success on each contingent edge
def DItriSTP_alpha(STN, updateQueue, triangleQueue, alpha):
  alpha = str(round(float(alpha),3))
  one_minus_alpha = str(round(1-float(alpha),3))
  p_alpha = {}
  for (i,j),edge in STN.contingentEdges.items():
    # TODO: Where does invCDF come from? It isn't referenced in this file or
    #       its dependencies.
    p_alpha[(i,j)] = 1000*invCDF[edge.distribution][one_minus_alpha]
    p_alpha[(j,i)] = -1000*invCDF[edge.distribution][alpha]

  numTris = 0
  while not updateQueue.empty() or not triangleQueue.empty():
    while not updateQueue.empty():
      update = updateQueue.get()
      e = STN.getEdge(update['i'], update['j'])
      #print e.getWeightMax()
      forward = STN.updateEdge(update['i'], update['j'], update['Wij'])
      reverse = STN.updateEdge(update['j'], update['i'], update['Wji'])

      if e.getWeightMin() > e.getWeightMax():
        #print "Inconsistant: %d->%d: [%f, %f]" %(e.i, e.j, e.getWeightMin(), e.getWeightMax())
        return False

      if forward or reverse:
        adjTri = STN.getAdjacentTriangles(update['i'], update['j'])
        for tri in adjTri:
          triangleQueue.put(tri)

    # Early termination if we didn't make any updates
    if triangleQueue.empty():
      break

    curTri = triangleQueue.get()
    numTris += 1

    # This is really sketchy but it works...
    if numTris > len(STN.tris)**2:
      return False

    updatedEdges = []
    for order in itertools.permutations(curTri.nodes()):
      i = order[0]
      j = order[1]
      k = order[2]

      weightSum = STN.getEdgeWeight(i, k) + \
                  STN.getEdgeWeight(k, j)

      if i==0 and STN.getEdge(k,j).isContingent() and STN.getEdgeWeight(j,k)>0:
        weightSum = min(weightSum, STN.getEdgeWeight(i,k)-p_alpha[(j,k)])

      if j==0 and STN.getEdge(i,k).isContingent() and STN.getEdgeWeight(i,k)>0:
        weightSum = min(weightSum, STN.getEdgeWeight(k,0)-p_alpha[(k,i)])

      if STN.updateEdge(i, j, weightSum):
        updatedEdges.append([i, j])

    for edge in updatedEdges:
      adjTri = STN.getAdjacentTriangles(edge[0], edge[1])
      for tri in adjTri:
        triangleQueue.put(tri)
      sendEdge = STN.getEdge(edge[0], edge[1])
      updateQueue.put({'i':sendEdge.i, 'j':sendEdge.j, \
                       'Wij':sendEdge.Cij, 'Wji':sendEdge.Cji})

  return True

