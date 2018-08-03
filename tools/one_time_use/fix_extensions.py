import argparse
import json, sys

## \file fixExtensions.py
#
#  \brief gets rid of periods at the end of file names for distributions in an STN
#


## \fn deletePeriods(stnFile)
#
#  \brief gets rid of periods at end of distribution names
def deletePeriods(stnFile):

  # goes through each edge and checks if it is contingent
  for edge in stnFile['constraints']:
    if 'distribution' in edge:
      if edge['distribution']['name'][-1] == '.':
        edge['distribution']['name'] = edge['distribution']['name'][:-1]

  return stnFile


def main():

  # handles arguments on command line
  parser = argparse.ArgumentParser(description='Get rid of the periods at the end of distribution names')

  parser.add_argument('jsonFile', metavar='JSON', type=str,
                      help='STN json file')
  args = parser.parse_args()


  with open(args.jsonFile) as f:
    data = json.load(f)

  # get rid of the periods
  modified_stn = deletePeriods(data)


  fileName = args.jsonFile.split("/")
  a,b = fileName[-1].split(".")

  # creates a file with the fake STN with the indicated change at the end of its filename
  f = open(args.jsonFile, 'w')
  f.write(json.dumps(modified_stn, sort_keys=False, indent = 2, separators=(',', ': ')))
  f.close()

if __name__ == '__main__':
  main()