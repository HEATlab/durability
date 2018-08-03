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
  minVal = float(sys.argv[2])
  maxVal = float(sys.argv[3])
  step = float(sys.argv[4])

  # load json STN to save it
  with open(sys.argv[1], 'r') as f:
    initData = json.load(f)

  with open(sys.argv[1], 'r') as f:
    flexData = json.load(f)


  makespan = minVal

  while makespan <= maxVal:
    print '\nmax = {}'.format(makespan)

    # reset our copy of the STN just to be safe
    with open(sys.argv[1], 'r') as f:
      flexData = json.load(f)


    # update our copy of the STN with the makespan
    for v in initData['nodes']:
      v['max_domain'] = makespan

    # dump the STN so our scripts can read it
    with open('json/flex_temp.json', 'w') as f:
      json.dump(initData, f)

    # First run whatever decoupling algorithm we want
    flex = Popen(['python', 'modifiedWilsonDecoupling.py', 'json/flex_temp.json'], env=os.environ, stdout=PIPE)
    flex.wait()

    flexComp = flex.stdout.read()
    #print flexComp
    flexMatch = re.search(ur'flexibility is: ([0-9\.]+)', flexComp)
    if flexMatch != None:
      print 'flex is {}'.format(flexMatch.group(1))


    # update our copy of the STN with the decoupling
    for v in flexData['nodes']:
      regMin = 't'+str(v['node_id'])+ur'Lo = ([0-9\.-]+)'
      regMax = 't'+str(v['node_id'])+ur'Up = ([0-9\.-]+)'
      minMatch = re.search(regMin, flexComp)
      maxMatch = re.search(regMax, flexComp)
      #print 'matched {0} at {1}'.format(regMin,minMatch.group(1))
      #print 'matched {0} at {1}'.format(regMax,maxMatch.group(1))

      v['max_domain'] = float(maxMatch.group(1))
      v['min_domain'] = float(minMatch.group(1))

    with open('json/decouple_temp.json','w') as f:
      json.dump(flexData, f)

    # Find the robustness of the decoupling
    simulation = Popen(['../stpsimulator/target/release/simulator_stp',
                                  '--samples', str(10000),
                                  '--sample_directory', '../stpsimulator/samples/',
                                  'json/decouple_temp.json'],
                                  env=os.environ, stdout=PIPE)
    simulation.wait()


    simData = simulation.stdout.read()
    match = dataRegex.search(simData)

    if match is None:
      print 'NO MATCH. This was the output'
      print simData
      sys.exit(1)

    print >> sys.stderr, '{0}: {1}%'.format(int(makespan), match.group(1))


    makespan += step

## \endcond
