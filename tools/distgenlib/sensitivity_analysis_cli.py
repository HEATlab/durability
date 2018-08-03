#!/usr/bin/env python


##
# \file sensitivity_analysis_cli.py
# \brief Command line interface for testing sensitivity analysis.
# \details This program is symlinked to the robotbrunch/bin directory (not
#   the user directory, mind you). Run
#   \code{.sh}
#   rb-sensitivity-harvest --help
#   \endcode
#   for usage.


import sys, os
import argparse
import shutil
from sensitivity_harvest import harvest_shifted, harvest_normalized
from harvest_distributions import SAMPLE_NUM, BANDWIDTH_DEFAULT

PROG_NAME = "rb-sensitivity-harvest"


##
# \fn main
# \brief Main function which is run when harvest distribution is run as a
#   standalone program.
def main():
  args = parse_args()
  # To understand how this line works see:
  # https://docs.python.org/3/library/argparse.html
  args.func(args)


## Perform a shifted harvest
def handle_shift(args):
  harvest_shifted(run_dir         = args.run_dir,
                  shift_min       = args.shift_min,
                  shift_max       = args.shift_max,
                  shift_step      = args.step_size,
                  shift_is_stddev = args.stddev,
                  kde_factor      = args.bandwidth,
                  sample_num      = args.samples,
                  silent          = args.silent)


## Perform a normalised harvest with variable standard deviations.
def handle_norm(args):
  harvest_normalized(run_dir             = args.run_dir,
                     stddev_min          = args.sd_min,
                     stddev_max          = args.sd_max,
                     stddev_step         = args.step_size,
                     stddev_factor_type  = args.factor_type,
                     sample_num          = args.samples,
                     silent              = args.silent)


##
# \fn parse_args
# \brief Parse command line arguments.
def parse_args():

  # This is more complicated parsing than normal. Check out this tutorial
  # for details: https://docs.python.org/3/library/argparse.html

  parser = argparse.ArgumentParser(prog=PROG_NAME, description='Sensitivity'
    +' experiment harvesting program.')
  subparsers = parser.add_subparsers(help='Use one of the following commands.'
    +' Each has its own help text if combined with --help.')

  # Shift command -------------------------------------------------------------
  parser_shift = subparsers.add_parser('shift',
                                       help='Generate shifted distributions')
  parser_shift.add_argument(
    'shift_min',
    metavar='SHIFT_MIN',
    type=float,
    help='Minimum to shift, as a float.')
  parser_shift.add_argument(
    'shift_max',
    metavar='SHIFT_MAX',
    type=float,
    help='Maximum to shift, as a float, inclusive.')
  parser_shift.add_argument(
    'step_size',
    type=float,
    metavar='STEP_SIZE',
    help='Step size to determine the density of shift values we want to test.')
  parser_shift.add_argument(
    'run_dir',
    type=str,
    metavar='RUNDIR',
    help="The robot run directory to harvest.")
  parser_shift.add_argument(
    '-w',
    '--bandwidth',
    metavar='BANDWIDTH',
    default=str(BANDWIDTH_DEFAULT),
    type=float,
    help=('KDE bandwidth factor. Default is {0}.').format(BANDWIDTH_DEFAULT))
  parser_shift.add_argument(
    '-s',
    '--stddev',
    action='store_true',
    help='Set shift to be in standard deviations instead of actual values.')
  parser_shift.add_argument(
    '--samples',
    metavar='COUNT',
    default=SAMPLE_NUM,
    type=int,
    help=('Set the number of samples to take per harvest. Default is {0}. See'
      + ' rb-harvest-distributions for more details.').format(SAMPLE_NUM))
  parser_shift.add_argument(
    '--silent',
    action='store_true',
    help="Don't print debug messages.")
  parser_shift.set_defaults(func=handle_shift)

  # Norm command --------------------------------------------------------------
  parser_norm = subparsers.add_parser('norm',
                                      help='Generate normal distributions')
  parser_norm.add_argument(
    'sd_min',
    type=float,
    metavar='SD_MIN',
    help='Minimum standard deviation factor to test, as a float.')
  parser_norm.add_argument(
    'sd_max',
    type=float,
    metavar='SD_MAX',
    help='Maximum standard deviation factor to test, as a float. Inclusive.')
  parser_norm.add_argument(
    'step_size',
    type=float,
    metavar='STEP_SIZE',
    help='Step size to determine the density of standard deviation values we'
      +' want to test.')
  parser_norm.add_argument(
    'factor_type',
    type=str,
    metavar='FACTOR_TYPE',
    help="How are the standard deviation factors applied to the normal"
      +" distributions? Must be either 'SET' or 'MULTIPLY'")
  parser_norm.add_argument(
    'run_dir',
    type=str,
    metavar='RUNDIR',
    help="The robot run directory to harvest.")
  parser_norm.add_argument(
    '--samples',
    metavar='COUNT',
    default=SAMPLE_NUM,
    type=int,
    help=('Set the number of samples to take per harvest. Default is {0}. See'
      + ' rb-harvest-distributions for more details.').format(SAMPLE_NUM))
  parser_norm.add_argument(
    '--silent',
    action='store_true',
    help="Don't print debug messages.")
  parser_norm.set_defaults(func=handle_norm)

  return parser.parse_args()


if __name__ == "__main__":
  main()


