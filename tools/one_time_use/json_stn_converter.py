#!/usr/bin/env python

import sys, json
import argparse

## \file json_stn_converter.py
#
#  \brief Converts an old format JSON STN to a new format JSON STN
#
#  \details Everything should now be moved over to the new format,
#           but this is here in case it needs to be used again.

## \fn convertJSON_old2new(inSTN)
#  \brief Does the conversion between old JSON (as a dictionary) and new JSON
def convertJSON_old2new(inSTN):
  outSTN = dict()

  outSTN['num_agents'] = inSTN['numAgents']

  if 'reset_pts' in inSTN.keys():
    outSTN['reset_pts'] = dict()
    for pt in inSTN['reset_pts']:
      outSTN['reset_pts'][pt] = inSTN['reset_pts'][pt]

  outSTN['nodes'] = []
  for v in inSTN['verts']:
    v2 = dict()
    v2['node_id'] = v['nodeID']
    v2['owner_id'] = v['ownerID']
    v2['min_domain'] = int(v['min'])
    v2['max_domain'] = int(v['max'])

    # these are for backwards compatibility
    v2['local_id'] = v['localID']
    v2['start'] = v['start']
    v2['location'] = v['location']

    outSTN['nodes'].append(v2)


  outSTN['constraints'] = []
  for e in inSTN['edges']:
    e2 = dict()
    e2['first_node'] = e['i']
    e2['second_node'] = e['j']

    if e['max'] == 'Infinity':
      e2['max_duration'] = 'inf'
    else:
      e2['max_duration'] = e['max']

    if e['min'] == '-Infinity':
      e2['min_duration'] = '-inf'
    else:
      e2['min_duration'] = e['min']

    if 'dist' in e.keys():
      e2['distribution'] = dict()
      e2['distribution']['type'] = 'Empirical'
      e2['distribution']['name'] = e['dist']

    outSTN['constraints'].append(e2)

  return outSTN

# This tells doxygen not to document all of the variables in the script below
## \cond SKIP
if __name__ == "__main__":

  infile = sys.argv[1]
  outfile = sys.argv[2]
  with open(infile, 'r') as f:
    inSTN = json.load(f)

  outSTN = convertJSON_old2new(inSTN)

  with open(outfile, 'w') as f:
    json.dump(outSTN, f, indent=2, separators=(',', ':'))
## \endcond
