#!/usr/bin/python
from subprocess import Popen, PIPE
import re, os, sys, json, copy
from json_stn_converter import convertJSON_old2new
from schedulase_stntools import *
from stntools import *
from guaranteedRobustness import guaranteedRobustness

## \file schedulaseSweep_stntools.py
#
#  \brief This file runs schedulase.py over a range of makespans
#
#  \details All of the data is piped to stderr, so to save the data you need to
#  pipe stderr to a file. Run the code as follows:
#  \code{.unparsed}
#  python schedulaseSweep.py JSONFILE MAKESPAN_LOW MAKESPAN_HIGH STEP [oldConstraints] [--intervalPlot PLOTFILE] 2> DATAFILE
#  \endcode
#  This file also runs the binary search on alpha, using the inverse CDFs
#  produced by `loadInverseCDFs.py`

# This tells doxygen not to document all of the variables in the script below
## \cond SKIP

dataRegex = re.compile(ur'result:\n([0-9\.]+)')

if __name__ == '__main__':
  usage = "USAGE: python schedulaseSweep.py JSONFILE MAKESPAN_LOW MAKESPAN_HIGH STEP [oldConstraints] [--intervalPlot PLOTFILE] [--guarantee PLOTFILE]"
  try:
    inSTN = loadSTNfromJSONfile(sys.argv[1])['stn']
    minVal = float(sys.argv[2])
    maxVal = float(sys.argv[3])
    step = float(sys.argv[4])
  except:
    print usage
    sys.exit(0)

  oldConstraints = False
  intervalPlot = False
  guarantee = False
  #check the rest of the command line arguments
  for i,arg in enumerate(sys.argv[5:]):

    if arg == "oldConstraints":
      print "using old constraints"
      oldConstraints = True

    elif arg == "--intervalPlot" and len(sys.argv) >= i:
      intervalPlotFile = sys.argv[5+i+1]
      try:
        with open(intervalPlotFile,"w+") as f:
          intervalPlot = True
      except:
        print "invalid file name: {0}".format(intervalPlotFile)
        print usage
        sys.exit(0)

    elif arg == "--guarantee" and len(sys.argv) >= i:
      guaranteeFile = sys.argv[5+i+1]
      try:
        with open(guaranteeFile,"w+") as f:
          guarantee = True
      except:
        print "invalid file name: {0}".format(guaranteeFile)
        print usage
        sys.exit(0)

  makespan = minVal


  while makespan <= maxVal:
    print '\n----------\nmakespan = {}'.format(makespan)

    # reset our copy of the STN
    mutSTN = inSTN.copy()
    mutSTN.setMakespan(makespan)

    result = schedulaseBinarySearch(mutSTN,debug=True,decouple=False)

    if not result:
      print "Could not decouple. Skip."
      print >> sys.stderr, '{0}: {1}% {2}'.format(makespan, 0, 1.0)
      makespan += step
      continue
    else:
      alpha,outSTN = result


    # The simulator requires file input
    outFile = os.path.dirname(os.path.abspath(__file__))+'/json/decouple_temp.json'
    with open(outFile,'w') as f:
      json.dump(outSTN.forJSON(), f, indent = 2, separators=(',',':'))


    if guarantee:
      stn_dict = loadSTNfromJSONfile(outFile)['stn']

      with open(guaranteeFile,"a+") as f:
        f.write("{} {}\n".format(makespan,guaranteedRobustness(STN)))

    # Find the robustness of the decoupling
    simulation = Popen(['../stpsimulator/target/release/simulator_stp',
                                  '--samples', str(10000),
                                  '--sample_directory', '../stpsimulator/samples/',
                                  outFile],
                                  env=os.environ, stdout=PIPE)
    simulation.wait()


    """
    # the old robusness algorithm, in case you want it
    simulation = Popen(['python','../robustness_simulator/sampling_sim.py',
                                  'json/decouple_temp.json',
                                  '10000'],
                                  env=os.environ, stdout=PIPE)

    simulation.wait()
    """

    # extract the robustness from simulator output
    simData = simulation.stdout.read()
    match = dataRegex.search(simData)

    if match is None:
      print 'NO REGEX MATCH. This was the output'
      print simData
      sys.exit(1)

    print >> sys.stderr, '{0}: {1:6.2f}% {2}'.format(int(makespan), float(match.group(1)), alpha)
    makespan += step

## \endcond
