##
# \file stnjsontools.py
#
# \brief Contains library functions for reading JSON files that represent STNs.
#


import json
from stn import STN
from stnreduce import stnreduce


##
# \fn loadSTNfromJSON
#
# \brief Wrapper for loadSTNfromJSONobj that loads an STN from a JSON string.
#
# \deprecated Notice, this function was quite horribly implemented the first
#   time around. A correction was made that breaks the entirely inconsistent
#   interface that was originally there.
#
# @param json_str   String representation of a full JSON object.
# @param reduction  Not sure what this does.
# @return           Returns the dictionary of the form...
#                   {'stn' : STN, 'agent_count' : numAgents}
def loadSTNfromJSON(json_str, reduction = True, using_PSTN=True):
  jsonSTN = json.loads(json_str)
  return loadSTNfromJSONobj(jsonSTN, reduction=reduction, using_PSTN=using_PSTN)


##
# \fn loadSTNfromJSONfile
# \brief Wrapper function for loadSTNfromJSON to allow reading from files.
#
# @param filepath Path of file to read in
# @param reduction Triangulate the STN, I guess?
# @return Returns a dictionary that has the format
#   {'stn' : STN, 'agent_count' : numAgents}
def loadSTNfromJSONfile(filepath, reduction=True, using_PSTN=True):
  with open(filepath,'r') as f:
    output_dict = loadSTNfromJSON(f.read(),reduction=reduction, using_PSTN=using_PSTN)
  return output_dict


##
# \fn loadSTNfromJSONobj
# \brief Returns a dictionary that gives us an STN from a JSON file, with
#   a few extra details.
#
# @param jsonSTN json object to use that represents an STN.
# @param reduction Triangulate the STN, I guess?
# @return Returns a dictionary that has the format
#   {'stn' : STN, 'agent_count' : numAgents}
def loadSTNfromJSONobj(jsonSTN, reduction=True, using_PSTN=True):
  stn = STN()

  # Add the root vertex and put it in the T_x set
  stn.addVertex(0, 0, None)
  # TODO: wtf? Why are we executing a point outside of a simulation in the
  # first place?
  stn.execute(0)
  agents = []

  # Add the vertices
  for v in jsonSTN['nodes']:
    # Accumulate a list of all the owners to retrieve the agents.
    if not v['owner_id'] in agents:
      agents.append(v['owner_id'])

    # We don't necessarily need a location, just set it to None if we don't.
    if not ('location' in v):
      v['location'] = None

    stn.addVertex(v['node_id'], v['local_id'], v['owner_id'], v['location'])

    # TODO: Change edge adding to allow integers to refer to vertecies,
    # rather than actually connecting integers together. (silly legacy support)
    stn.addEdge(0, v['node_id'], float(v['min_domain']), float(v['max_domain']))
    if 'executed' in v:
      if v['executed']:
        stn.execute(v['node_id'])

  # Add the edges
  for e in jsonSTN['constraints']:
    if 'distribution' in e and using_PSTN:
        stn.addEdge(e['first_node'],
                           e['second_node'],
                           float(e['min_duration']),
                           float(e['max_duration']),
                           e['distribution']['name'])
    else:
      stn.addEdge(e['first_node'],
                         e['second_node'],
                         float(e['min_duration']),
                         float(e['max_duration']))

  #numAgents = jsonSTN['num_agents']

  stn.numAgents = len(agents) # NOTE: deprecated, agents replaces numAgents.
  stn.agents = agents

  if reduction:
    # Triangulate the STN
    stnCopy = stn.copy()
    agentWait = []
    # Set up the wait timer for agent load balancing
    for a in stn.agents:
      agentWait.append(0)

    # Perform the reduction
    while len(stnCopy.verts) > 1:
      for a in stn.agents:
        if agentWait[a] == 0:
          created = stnreduce(stnCopy, a, stn)
          agentWait[a] = created
        else:
          agentWait[a] -= 1

  # Return back an dictionary for easy labelling of return types.
  # FIXME: Remove the deprecated numAgents value
  output_dict = {'stn' : stn, 'agent_count': stn.numAgents}

  # Ideally, this would return a class, however, this is already butchering
  # quite a bit of legacy support, and a class would end up being even more
  # damaging.
  return output_dict

