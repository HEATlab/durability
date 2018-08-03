#!/usr/bin/env python

##
# \file Extendable plotting command line interface
# \details This program allows for an extendable plotting interface
#   to allow the developers of robot brunch to plot any data with ease.
#   Using
#
# \note
#   How to extend this program:
#   1. Add a subparser to the parse_args function that describes what the
#      command interface would look like.
#   2. Make sure you have the set_defaults(func=foobar) call.
#   3. Make a function with the same name as foobar.
#   4. Have this f function do all the plotty things (e.g. call Gnuplot)
#

import os, sys
import argparse
import subprocess

import linepointplot

# use CURDIR to find tools directory, which we then add to the python path.
CURDIR = os.path.dirname(os.path.realpath(__file__))
PROG_NAME = 'rb-plot'

HIST_PROG = 'hist.plt'
LINE_PROG = 'line.plt'


def main():
  args = parse_args()
  args.func(args)


##
# \fn handle_histogram
# \brief Handle histogram plotting through a subprocess.
def handle_histogram(args):

  # Get the absolute file path of the thing we want to plot.
  fp = os.path.abspath(args.fp)
  gnu_args = (
    "-p"
    +" -e \"binwidth=\'"+args.binwidth+"\'\""
    +" -e \"filename=\'"+fp+"\'\""
    +" "+CURDIR+"/"+HIST_PROG
  )
  subprocess.Popen('gnuplot '+gnu_args,shell=True)


##
# \fn handle_multiline
# \brief Plot multiple scatter plot lines
def handle_linepoint(args):
  linepointplot.generate_script(filepaths = args.fps,
                                xlabel = args.xlab,
                                ylabel = args.ylab,
                                errorbars = args.errorbars)

##
# \fn parse_args
# \brief Parse command line arguments.
def parse_args():

  # This is more complicated parsing than normal. Check out this tutorial
  # for details: https://docs.python.org/3/library/argparse.html

  parser = argparse.ArgumentParser(prog=PROG_NAME, description='Plotting'
   + ' utility.')
  subparsers = parser.add_subparsers(help='Use one of the following commands.'
    +' Each has its own help text if combined with --help.')

  # histogram command ---------------------------------------------------------
  parser_hist = subparsers.add_parser('histogram',
                                       help='Plot a single histogram.')
  parser_hist.add_argument(
    'fp',
    metavar='FILE',
    type=str,
    help='File to plot histogram of.')
  parser_hist.add_argument(
    '-w',
    '--binwidth',
    metavar='WIDTH',
    type=str,
    default="0.5",
    help='Binwidth to use, defaulted to 0.5.')
  parser_hist.set_defaults(func=handle_histogram)

  # multiline command ---------------------------------------------------------
  parser_line = subparsers.add_parser('linepoint',
                                       help='Plot multiple scatter plots given'
                                       +' xy data, which are connected via'
                                       +' lines.')
  parser_line.add_argument(
    '-x','--xlab',
    metavar='X_LAB',
    type=str,
    default="x",
    help='X Axis label')
  parser_line.add_argument(
    '-y','--ylab',
    metavar='Y_LAB',
    type=str,
    default="y",
    help='Y Axis label')
  parser_line.add_argument(
    '-e','--errorbars',
    metavar="(x|y|xy)",
    type=str,
    default=None,
    help="Error bar settings. This should be either one or two characters. It"
    +" indicates what error bar columns exist. Allowable settings include"
    +" 'x', 'y', or 'xy'."
  )
  parser_line.add_argument(
    'fps',
    metavar='FILE',
    type=str,
    nargs='+',
    help='Files to plot.')
  parser_line.set_defaults(func=handle_linepoint)

  return parser.parse_args()


if __name__ == "__main__":
  main()
