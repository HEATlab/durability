#!/usr/bin/env python

##
# \file harvestcli.py
# \brief Command line interface for doing a harvest manually.
# \details This program is symlinked to the robotbrunch/bin directory (not
#   the user directory, mind you). Run
#   \code{.sh}
#   rb-harvest-distributions --help
#   \endcode
#   for usage.

import argparse
from harvest_distributions import harvest, SAMPLE_NUM, BANDWIDTH_DEFAULT


##
# \fn main
# \brief Main function which is run when harvest distribution is run as a
#   standalone program.
def main():
  args = parse_args()

  # Actually run generation
  harvest(run_dir=args.run_dir,
          kde_factor=args.bandwidth,
          offset=args.off,
          silent=args.silent,
          normal=args.normal,
          suffix=args.suffix,
          sample_num=args.samples,
          offset_is_stddev=True)


##
# \fn parse_args
# \brief Parse command line arguments.
def parse_args():
  # Handle argument parsing
  parser = argparse.ArgumentParser(description='Harvest And Generate New'
    + ' Distributions')
  parser.add_argument('-n','--normal',
                      action='store_true',
                      help='Generate normal distributions instead of kernel'+
                      ' density distributions.')
  parser.add_argument('-s','--suffix',
                      metavar='SUFF',
                      default='',
                      type=str,
                      help='Suffix to append to the generated file names.'
                      +' Common values include "_shift" and "_norm"')
  parser.add_argument('-w','--bandwidth',
                      metavar='BANDWIDTH',
                      default=str(BANDWIDTH_DEFAULT),
                      type=float,
                      help=('KDE bandwidth factor. Default is {}. Not used if'
                      +' distribution is normal.').format(BANDWIDTH_DEFAULT))
  parser.add_argument('--off',
                      metavar='OFFSET',
                      default='0.0',
                      type=float,
                      help='Offset any generated distributions by some amount.')
  parser.add_argument('-c', '--samples',
                      metavar='COUNT',
                      default=SAMPLE_NUM,
                      type=int,
                      help='Number of samples to generate with. Default is'
                        +' {}.'.format(SAMPLE_NUM))
  parser.add_argument('--silent',
                      action='store_true',
                      help='Prevent any print statements from running.')
  parser.add_argument('run_dir', metavar='RUNDIR', type=str,
                      help='Directory containing samples directory to generate'\
                            +' into.')

  return parser.parse_args()


if __name__ == "__main__":
  main()

