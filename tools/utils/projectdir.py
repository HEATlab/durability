#!/usr/bin/env python

import os

##
# \file brunchpath.py
#
# \brief Utility function to get the root directory of robotbrunch.
#
# \details Anytime we want to find a file relative to the robotbrunch directory
#   this is how we should get it.
#
#

# Name of the main directory.
MAIN_DIR = 'robotbrunch'
# Directory delimiters
DIR_DELIM = '/'
# File system root directory. (on linux this is / )
SYSTEM_ROOT = '/'

##
# \fn projectdir
#
# \brief Allows converting
#
# @param internal_path Internal path to follow below the main directory.
#   Can be omitted to get the main directory. Does not require a delimiter
#   first. Also allows a list, in which case each element of the list will be
#   used as an input.
# @return Returns a converted path relative to the main directory to an absolute
#   path. If the input was a list, then a list of relative paths.
def projectdir(internal_path=''):
  # Allow polymorphism by just applying this function recursively
  # to each element of the list.
  if type(internal_path) == list:
    return [projectdir(path) for path in internal_path]

  # Get the location of this file
  script_path = os.path.realpath(__file__)
  # Assuming we are on a POSIX system, as opposed to say, Windows.
  dir_list = script_path.split(DIR_DELIM)

  parent_dir_count = 0

  for i in range(len(dir_list)):

    # If we haven't raised high enough in the directory tree,
    if dir_list[-(i+1)] == MAIN_DIR:
      break
    else:
      parent_dir_count += 1

    if i == (len(dir_list)-1):
      raise OSError("Could not find {} main directory.".format(MAIN_DIR))

  # Slice off the correct path.
  parent_dir_list = dir_list[1:-parent_dir_count]
  brunch_dir = SYSTEM_ROOT + DIR_DELIM.join(parent_dir_list) + DIR_DELIM

  # Return the wanted directory
  return brunch_dir + internal_path
