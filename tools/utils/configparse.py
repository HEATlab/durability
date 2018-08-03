#!/usr/bin/env python

##
# \fn configparse
# \brief Utility function to read config files.
#
# @param filepath The config filepath to read.
# @return Returns a list of strings where each element is a line in the
#         config file, given that line is not a comment
#         e.g. (starts with a #)
def configparse(filepath):
  # filepath must be a string.
  if type(filepath) != str:
    raise TypeError

  output_list = []

  with open(filepath) as f:
    for l in f:
      if l[0] != "#" and l[0] != "\n":
        # Go to the second to last char,
        # because we don't want to include
        # the new line.
        output_list.append(l[:-1])

  return output_list
