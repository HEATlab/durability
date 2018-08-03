from stntools import *
import sys,os,re
from subprocess import Popen, PIPE

def main():
  usage = 'USAGE: {} JSON_FILE'.format(sys.argv[0])
  try:
    mySTN = loadSTNfromJSONfile(sys.argv[1])['stn']
  except:
    print 'Could not load mySTN from JSON file: {}'.format(sys.argv[1])
    print usage
    raise

  bestSTN = mySTN.minimize_binary_search()

  if bestSTN == None:
    print "minimization failed"
  else:
    print bestSTN
    return
    # The simulator requires file input
    outFile = os.path.dirname(os.path.abspath(__file__))+'/json/decouple_temp.json'
    with open(outFile,'w') as f:
      json.dump(bestSTN.forJSON(), f, indent = 2, separators=(',',':'))
    # Find the robustness of the decoupling
    simulation = Popen(['../stpsimulator/target/release/simulator_stp',
                                  '--samples', str(10000),
                                  '--sample_directory', '../stpsimulator/samples/',
                                  outFile],
                                  env=os.environ, stdout=PIPE)
    simulation.wait()
    # extract the robustness from simulator output
    simData = simulation.stdout.read()
    dataRegex = re.compile(ur'result:\n([0-9\.]+)')
    match = dataRegex.search(simData)

    print "Robustness: {:6.2f}%".format(float(match.group(1)))

if __name__ == "__main__":
  main()
