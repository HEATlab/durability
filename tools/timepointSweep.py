from subprocess import Popen, PIPE
import re, os, sys, json

## \file flexSweep.py
#
#  \brief This file runs simpleLinearDecoupling.py over a range of makespans
#
#  \details All of the data is piped to stderr, so to save the data you need to
#  pipe stderr to a file. Run the code as follows:
#  \code{.unparsed}
#  python flexSweep.py JSONFILE MAKESPAN_LOW MAKESPAN_HIGH STEP 2> DATAFILE
#  \endcode
#

# This tells doxygen not to document all of the variables in the script below
## \cond SKIP

dataRegex = re.compile(ur'result:\n([0-9\.]+)')

if __name__ == '__main__':
  json_file = sys.argv[1]
  minVal = int(sys.argv[2])
  maxVal = int(sys.argv[3])
  step = int(sys.argv[4])



  makespan = minVal

  while makespan <= maxVal:
    print '\nmax = {}'.format(makespan)

    # Find the robustness of the decoupling
    simulation = Popen(['../stpsimulator/target/release/timepoint_sim',
                                  '--makespan', str(makespan),
                                  '--samples', '10000',
                                  '--sample_directory', '../stpsimulator/samples/',
                                  json_file],
                                  env=os.environ, stderr=sys.stderr)
    simulation.wait()

    makespan += step

## \endcond
