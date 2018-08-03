#!/usr/bin/env python

import os
import sys
import argparse

from dist_generator import dist_generator

##
# \file harvest_distributions.py
#
# \brief Collects and generates PDF and inverse CDF distributions for a run.
#
# \details
# Similar to generate_distribution.py, this program will generate PDF and
# inverse CDF data files given a
#
# However, this program takes in a full run directory, and harvests a full set
# of run data in one go. See `harvest` documentation for details on the structure
# of the run directory.
#
# Notice however, that this program provides the extensions ".samples" and
# ".invcdf" to the generated files, which is different behaviour than previous
# Mathematica scripts.
#


## Default number of samples to use during generation
SAMPLE_NUM = 10000

## The default KDE bandwidth factor.
BANDWIDTH_DEFAULT = 0.1

## label for the raw subdirectory
RAW_SUBDIR = "raw"

## label for the samples subdirectory
SAMPLE_SUBDIR = "samples"

## label for the inverse CDF subdirectory
INVCDF_SUBDIR = "inverseCDFs"

## Extension for samples files.
SAMPLE_EXT = ".samples"

## Extension for inverse cdf files.
INVCDF_EXT = ".invcdf"

## Extension for raw data files.
# Note, this is solely for reading raw data files.
# Writing raw data files happens in robot_main code.
RAW_EXT = ".dat"


##
# \fn harvest
#
# \brief Generate all PDF and inverse CDF files given a run directory.
#
# \details
# The run directory should have the following structure:
#
# \verbatim
#   run_dir/
#   |
#   |- raw/
#   |  |- ff.dat
#   |  |- ft.dat
#   |  |- etc...
#   |- samples/
#   |- inverseCDFs/
# \endverbatim
#
# However, the directories samples/ and inverseCDFs/ need not exist for this
# script to function. If they don't exist, they will be generated.
#
# @param run_dir Run directory.
# @param kde_factor Kernel Density Bandwidth Factor, default is 0.1.
# @param offset Amount to shift distributions. Default is zero.
# @param silent If true, don't print any messages (Default True)
# @param normal If set to true, all generated files will be normal
#   distributions instead of KDE distributions.
# @param suffix String added before the file extension
# @param sample_num Number of samples to generate. Default 10000.
# @param offset_is_stddev Is the offset in std deviations? Default is False.
# @post Generated files in samples directory and inverseCDFs directory.
def harvest(run_dir, kde_factor=BANDWIDTH_DEFAULT, offset=0, silent=True,
            normal=False, suffix='', sample_num=SAMPLE_NUM,
            offset_is_stddev=False, normal_args={}):

  # Initialise generator
  dist_gen = dist_generator()

  # Get the relative path to the directories
  raw_path = os.path.join(run_dir, RAW_SUBDIR)
  sample_path = os.path.join(run_dir, SAMPLE_SUBDIR)
  invcdf_path = os.path.join(run_dir, INVCDF_SUBDIR)

  # Retrieve all the files from the raw directory
  raw_files = [ \
    f for f in os.listdir(raw_path) \
    if os.path.isfile(os.path.join(os.path.join(raw_path,f))) \
  ]

  for f in raw_files:
    # Get file path and extension
    path = os.path.join(raw_path,f)
    name, file_ext = os.path.splitext(path)
    name = os.path.basename(name)

    # Make sure that we only grab files with a .dat extension
    if (file_ext != RAW_EXT):
      continue

    # Run the generator -----------------------------------
    dist_gen.read_dist(path)

    # Apply offsets to the the distributions
    if offset != 0:
      if offset_is_stddev:
        # If the inputted offset is in terms of standard deviations,
        # use that as the offset "unit"
        shift = dist_gen.distsd() * offset
      else:
        shift = offset
      dist_gen.shift(shift) # Apply offset if set

    if normal:
      # Some times we want to generate a normal distribution
      # with non-standard standard deviations (heh).
      # normal_args allows a set of keys and values that
      # can be used as a way of handling these cases.
      if not 'stddev_factor_type' in normal_args:
        # We're not using normal_args, likely. Act normally then (heheh).
        # Generate internal normal distribution, with the empirical data's
        # standard deviation and mean.
        dist_gen.gen_normal()
      else:
        # We're actually using normal_args now, check how we want to use it's
        # values.
        if normal_args['stddev_factor_type'] == 'MULTIPLY':
          sd = normal_args['stddev_factor'] * dist_gen.distsd()
        elif normal_args['stddev_factor_type'] == 'SET':
          sd = normal_args['stddev_factor']
        else:
          raise ValueError(
            ('normal_args stddev_factor_type was set to {0}'
              +', which is not a valid type.').format(
                normal_args['stddev_factor_type'])
          )

        # Now generate the normal distribution with a different standard dev.
        dist_gen.gen_normal(sd=sd)
      # End normal case.
    else:
      dist_gen.gen_pdf(kde_factor=kde_factor) # Generate internal PDF using KDE

    dist_gen.resample(sample_num)
    new_distribution = dist_gen.pdf_samples # retrieve resampled distribution
    max_value = max(new_distribution)

    dist_gen.gen_cdf(0,max_value,sample_num)
    dist_gen.invert_discrete_cdf(sample_num)
    # -----------------------------------------------------

    # Get new file names for generated files
    sample_file_name = name + suffix + SAMPLE_EXT # Add extensions.
    invcdf_file_name = name + suffix + INVCDF_EXT
    sample_file_path = os.path.join(sample_path, sample_file_name)
    invcdf_file_path = os.path.join(invcdf_path, invcdf_file_name)

    ensure_dir_exists(sample_path)
    ensure_dir_exists(invcdf_path)

    # Write generated PDF and CDF
    dist_gen.write_sample_file(sample_file_path)
    dist_gen.write_invcdf_file(invcdf_file_path)

    if not silent:
      print("Gathered {}".format(path))
      print("  samples -> {}".format(sample_file_path))
      print("  invcdf -> {}".format(invcdf_file_path))


## Ensure a directory exists
def ensure_dir_exists(path):
  if not os.path.exists(path):
    os.makedirs(path)


##
# \fn sample_map
# \brief Retrieves data for contingent edge sampling from a run directory.
# \details
#
# \code{.py}
# samp_map = sample_map(some_directory)
# # ... or ...
# samp_map = sample_map([some_directory_1, some_directory_2, ...])
# \endcode
#
# @param run_dir A directory that represents a robot run. This directory must
#   have a subdirectory with the name "samples". Alternatively, it will accept
#   a list of such directories (as strings), where it then acts like each of
#   these are paths to run directories.
# @return Return a dictionary of the form {file_name : generated_sample_list}
def sample_map(run_dirs):

  # Initialise the sample_map.
  sample_map = {}

  # -----------------------------------
  # We accept both a string and a list as a parameter.
  # Eh, polymorphism is weird, but it helps us in the long run.
  # As long as it's properly documented.
  if type(run_dirs) == str:
    dirlist = [run_dirs]
  elif type(run_dirs) == list:
    dirlist = run_dirs
  else:
    raise TypeError("run_dirs was not a string nor a list.")
  # -----------------------------------

  for run_dir in dirlist:
    sample_path = os.path.join(run_dir,SAMPLE_SUBDIR)
    sample_files = os.listdir(sample_path)\


    for name in sample_files:
      # Get file path and extension
      path = os.path.join(sample_path,name)
      _, file_ext = os.path.splitext(path)
      name_only, _ = os.path.splitext(name)

      # Note, .dat is the legacy file extension.
      if file_ext != SAMPLE_EXT and file_ext != ".dat":
        continue

      # Open the files, and assign the key values to lists of samples.
      with open(path, 'r') as f:
        # Convert each
        sample_map[name_only] = [ float(line) for line in f.read().split("\n") \
                                  if line != ""]

  return sample_map


##
# \fn invcdf_map
# \brief
#
# @param run_dirs Either a string representing a directory path
#   or a list of strings representing directory paths that lead to a
#   "run directory", which are created when a robot logs a run.
#   Each run directory must contain a folder called "invCDFs" to be properly
#
# @return Returns a dictionary of the form [file_name : invcdf]
#   where invcdf itself is another dictionary, of the form
#   invcdf = {str(probability) : time}
#   Oh god why world.
def invcdf_map(run_dirs):
  # Initialise the map.
  invcdf_map = {}

  # -----------------------------------
  # We accept both a string and a list as a parameter.
  # Eh, polymorphism is weird, but it helps us in the long run.
  # As long as it's properly documented.
  if type(run_dirs) == str:
    dirlist = [run_dirs]
  elif type(run_dirs) == list:
    dirlist = run_dirs
  else:
    raise TypeError("run_dirs was not a string nor a list.")
  # -----------------------------------

  for run_dir in dirlist:
    invcdf_path = os.path.join(run_dir, INVCDF_SUBDIR)
    invcdf_files = os.listdir(invcdf_path)

    for name in invcdf_files:
      # Get file path and extension
      path = os.path.join(invcdf_path,name)
      _, file_ext = os.path.splitext(path)
      name_only, _ = os.path.splitext(name)

      # Only grab sample files
      if file_ext != INVCDF_EXT and file_ext != ".dat":
        continue

      # Open the files, and assign the key values to lists of samples.
      with open(path, 'r') as f:
        invcdf_map[name_only] = {}
        for l in f.readlines():

          # Split the lines of the inverseCDF CSV file.
          prob = float(l.split(',')[0])
          time = float(l.split(',')[1].strip('\n'))

          # Due to terrible legacy support, the keys of this dictionary are
          # the counted distribution probabilities from inverse CDFs as rounded
          # to 3 decimal places as strings.
          key = str(round(prob,3))

          # Yes, we are assigning a dictionary value to another dictionary.
          invcdf_map[name_only][key] = time

  return invcdf_map
