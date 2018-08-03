#!/usr/bin/python
from optparse import OptionParser, OptionGroup
from dc_stn import *
#from loadInverseCDFs import invCDF

# TODO: Change the organisation of the code such that this path append is not
# needed.
import sys,os
sys.path.append(os.getenv("HOME")+"/robotbrunch/robustness_simulator/")
import harvest_distributions
from sampling_sim import runSimulation
from math import ceil, floor
import harvest_distributions

INVCDF = harvest_distributions.invcdf_map(os.getenv('HOME')+'/robotbrunch/data/calibration/normal')
##
# \file dc_binary_search.py
#
# \brief Finds the optimal risk level alpha such that the equivalent STNU is DC
# \note Very similar to the algorithm presented by Tsamardinos et. al. 2002
#
def sc_binary_search(STN,debug=False,lb = 0.0,ub = 1.0):
  # dictionary of alphas for binary search
  alphas = {i:i/1000.0 for i in range(1001)}

  # bounds for binary search
  lower = ceil(lb * 1000) - 1
  upper = floor(ub * 1000) + 1

  best_alpha = None
  best_SC_STN = None

  # Run binary search on alpha
  while upper - lower > 1:
    #set up the STN at risk level alpha
    alpha = alphas[(upper + lower) // 2]
    if debug:
      print 'trying alpha = {}'.format(alpha)
    set_contingent_bounds(STN,alpha)

    #check SC
    is_sc,sc_STN = STN.isStronglyControllable(returnSTN=True)
    if is_sc:
      best_alpha = alpha
      best_SC_STN = sc_STN
      upper = (upper + lower) // 2
      if debug:
        print "STNU was SC"
    else:
      lower = (upper + lower) // 2
      if debug:
        print "STNU was not SC"

  return best_SC_STN,best_alpha


def dc_binary_search(STN,debug=False,lb = 0.0,ub = 1.0):
  # dictionary of alphas for binary search
  alphas = {i:i/1000.0 for i in range(1001)}

  # bounds for binary search
  lower = ceil(lb * 1000) - 1
  upper = floor(ub * 1000) + 1

  best_alpha = 1
  best_DC_STN = None
  robustness = 0

  # Run binary search on alpha
  while upper - lower > 1:
    #set up the STN at risk level alpha
    alpha = alphas[(upper + lower) // 2]
    if debug:
      print 'trying alpha = {}'.format(alpha)
    set_contingent_bounds(STN,alpha)
    DC_STN = STNtoDCSTN(STN)

    #check DC
    if DC_STN.is_DC():
      best_alpha = alpha
      best_DC_STN = DC_STN
      upper = (upper + lower) // 2
      if debug:
        print "STNU was DC"
    else:
      lower = (upper + lower) // 2
      if debug:
        print "STNU was not DC"

  robustness = 100*(1-best_alpha)**len(STN.contingentEdges)
  return best_DC_STN,best_alpha,robustness


## \fn set_contingent_bounds(inputSTN,alpha)
#  \brief Sets all of the contingent edges in the input STN to a risk level of
#   alpha
def set_contingent_bounds(inputSTN,alpha):
  lower_bound = str(round(float(alpha)/2,3))
  upper_bound = str(round(1-float(alpha)/2,3))

  for (i,j),edge in inputSTN.contingentEdges.items():
    edge.Cij =  ceil(1000*INVCDF[edge.distribution][upper_bound])
    edge.Cji = ceil(-1000*INVCDF[edge.distribution][lower_bound])


##
#
# \fn getRobustness(STN)
# \brief Calls the rust simulator to compute the robustness of the input STN
#
def getRobustness(STN):
  tempSTN = 'json/temp.json'

  with open(tempSTN, 'w+') as f:
    json.dump(STN.forJSON(), f, indent=2, separators=(',', ':'))

  # Find the robustness of the decoupling
  simulation = Popen(['../stpsimulator/target/release/simulator_stp',
                                '--samples', str(10000),
                                '--threads', '4',
                                '--sample_directory', '../stpsimulator/samples/',
                                tempSTN],
                                stdout=PIPE, stderr=PIPE)
  simulation.wait()

  dataRegex = re.compile(ur'result:\n([0-9\.]+)')

  # extract the robustness from simulator output
  simData = simulation.stdout.read()
  match = dataRegex.search(simData)

  # simulator failed -- 0 robustness
  if match is None:
    return 0

  return float(match.group(1))


## \cond SKIP
if __name__ == '__main__':

  # handle command line input
  usage = "usage: python %prog [options] JSONFILE"
  parser = OptionParser(usage)

  group = OptionGroup(parser,"Simulation Options")
  group.add_option("-d","--debug", action = "store_true", dest = "debug",
          help = "Print debugging info.")
  parser.add_option_group(group)

  group = OptionGroup(parser, "Schedule Options")
  group.add_option("-m", "--makespan", type = "int", dest = "makespan",
          help = "Set the makespan for the schedule")
  parser.add_option_group(group)

  options,args = parser.parse_args()

  if len(args) != 1:
    parser.error("Incorrect number of arguments")

  #Load the file
  if options.debug:
    print "Loading STN"
  STN = loadSTNfromJSONfile(args[0])['stn']
  if options.makespan != None:
    STN.setMakespan(options.makespan)
    if options.debug:
      print "Set makespan to {}".format(options.makespan)

  DC = True

  if DC:
    DC_STN,alpha,robustness = dc_binary_search(STN,options.debug)
    if DC_STN == None:
      print "No dynamically controllable STNU found"
    else:
      print DC_STN
      print "The input PSTN was DC with a minimum risk level of {}".format(alpha)
      print "This corresponds to a robustness of %4.2f"%robustness
  else:
    SC_STN,alpha = sc_binary_search(STN,options.debug)
    if SC_STN == None:
      print "No strongly controllable STNU found"
    else:
      SC_STN.minimize()
      print SC_STN
      print "The input PSTN was SC with a minimum risk level of {}".format(alpha)
      robustness = 100*(1-alpha)**len(STN.contingentEdges)
      print "This corresponds to a robustness of %4.2f"%robustness
      # TODO: Change how samp_map is selected.
      # It probably should be a command line argument, or use some global path.
      samp_map = harvest_distributions.sample_map(os.getenv("HOME") \
        + "/robotbrunch/data/calibration/run_2")
      new_STN,sim_robustness = runSimulation(SC_STN, samp_map)
      print "The simulated robustness is %4.2f"%sim_robustness

## \endcond
