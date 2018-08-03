#!/usr/bin/env python

##
# \file linepointplot.py
# \brief Plotting program that generates a gnuplot of a set of lines.
#

import os
import subprocess

# =============================================================================

TEMP_STORAGE = "/tmp/"
SCRIPT_NAME = "multiconnectplot{0}.plt"

## A HERE doc for writing the gnuplot script
SCRIPT = (
"""
#!/bin/gnuplot
set datafile separator ","
set grid
{0}
{1}
{2}
""").strip()

## The gnuplot variable assignment command
ASSIGNMENT_COM = "v{0}='{1}'"
## The gnuplot plotting command
PLOT_BASE = "plot"
PLOT_ADD = (" v{0} with linespoints linetype 1 pointtype 7"
            +" linecolor rgb '{1}' title '{2}'{3}")

## List of colours to plot lines with. If we run out of colours, you'll have to
#  add some more. :c
COLORS = ('#892fce',
          '#1e8c16',
          '#c64545',
          '#3f6c9b',
          '#664d1b',
          '#000000')

# =============================================================================

##
# \fn generate_script
# \brief Program that generates a temporary gnuplot script that will generate a
#   pretty plot given a list of CSVs.
def generate_script(filepaths, xlabel='x', ylabel='y', errorbars=None):
  sets = ""
  assigns = ""
  plots = PLOT_BASE

  if xlabel != '':
   sets += "set xlabel '{}'\n".format(xlabel)
  if ylabel != '':
   sets += "set ylabel '{}'\n".format(ylabel)

  for i, path in enumerate(filepaths):
    # Create the variable assignment command.
    assign_com = ASSIGNMENT_COM.format(str(i),path)
    # Get error bars
    errorbar_com = ""

    if errorbars == "x":
      errorbar_com = ", '' with xerrorbars notitle rgb '{0}'"
    elif errorbars == "y":
      errorbar_com = ", '' with errorbars notitle rgb '{0}'"
    elif errorbars == "xy":
      errorbar_com = ", '' with xyerrorbars notitle rgb '{0}'"

    # Assign colours of error bars.
    if errorbar_com != "":
      errorbar_com = errorbar_com.format(COLORS[i])

    name = os.path.splitext(os.path.basename(path))[0]
    # Create the base for the plot command.
    plots += " " + PLOT_ADD.format(str(i),COLORS[i],name,errorbar_com)


    assigns += assign_com + "\n"

    # Add a comma if this is not the last line to graph.
    if i != len(filepaths) - 1:
      plots += ","
  # End for loop -------------------------------------------------------------

  new_script = SCRIPT.format(sets,assigns,plots)
  # Since we want every temp file to be unique, append the unique process ID to
  # the file name. This is standard practice, ask Geoff Kuenning.
  script_loc = TEMP_STORAGE + SCRIPT_NAME.format(os.getpid())

  # Write the script in a temporary location.
  with open(script_loc, 'w+') as s:
    s.write(new_script)

  # Add the -p to make the graph persist and halt the python program.
  gnu_args = ("-p "+script_loc)
  subprocess.call('gnuplot '+gnu_args,shell=True)

  # We just made a new file, delete it after the process is done.
  #os.remove(script_loc)
