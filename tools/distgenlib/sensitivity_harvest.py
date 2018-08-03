#!/usr/bin/env python

##
#
# \file sensitivity_harvest.py
# \brief Wrapper for harvest_distributions to allow mass distribution alterations.
#

import os, sys
from harvest_distributions import harvest, BANDWIDTH_DEFAULT, SAMPLE_NUM

# Suffixes for filenames
SHIFT_SUFF = "--shift_"
SHIFT_SD_SUFF = "--shiftsd_"
NORMM_SUFF = "--normm_"
NORMS_SUFF = "--norms_"

##
# \fn harvest_shifted
# \brief Systematically generates a set of shifted distributions.
#
# @param run_dir
# @param shift_min The minimum offset for the distribution shifts.
# @param shift_max The maximum offset for the distributions, inclusive.
# @param shift_step The step size for the offsets to generate.
# @param kde_factor Bandwidth factor for kernel density estimations.
# @param sample_num Number of samples to generate.
def harvest_shifted(run_dir,
                    shift_min,
                    shift_max,
                    shift_step,
                    shift_is_stddev=True,
                    kde_factor=BANDWIDTH_DEFAULT,
                    sample_num=SAMPLE_NUM,
                    silent = True):

  current_shift = shift_min
  while current_shift <= shift_max and current_shift >= shift_min:
    # Denote whether or not we used the standard deviation.
    if shift_is_stddev:
      suffix = SHIFT_SD_SUFF
    else:
      suffix = SHIFT_SUFF

    # Append the type of shift to the suffix.
    suffix += str(current_shift).replace(".","_")

    harvest(run_dir          = run_dir,
            kde_factor       = kde_factor,
            offset           = current_shift,
            silent           = silent,
            normal           = False,
            suffix           = suffix,
            sample_num       = sample_num,
            offset_is_stddev = shift_is_stddev)

    current_shift += shift_step
  # while loop end


##
# \fn harvest_shifted
# \brief Systematically shif
# @param run_dir
# @param stddev_min The minimum standard deviation for the normal dists.
# @param stddev_max The maximum standard deviation for the normal dists,
#   inclusive.
# @param stddev_step The step size for standard dev files.
# @param stddev_factor_type
# @param sample_num Number of samples to generate.
def harvest_normalized(run_dir,
                       stddev_min,
                       stddev_max,
                       stddev_step,
                       stddev_factor_type,
                       sample_num=SAMPLE_NUM,
                       silent = True):

    current_sd = stddev_min
    while current_sd <= stddev_max and current_sd >= stddev_min:
      if stddev_factor_type.upper() == "MULTIPLY":
        suffix_prepend = NORMM_SUFF
      elif stddev_factor_type.upper() == "SET":
        suffix_prepend = NORMS_SUFF
      else:
        raise ValueError("stddev_factor_type was not recognised.")

      suffix = (suffix_prepend + str(current_sd)).replace(".","_")
      norm_args = {
                   'stddev_factor_type': stddev_factor_type,
                   'stddev_factor': current_sd
                  }
      harvest(run_dir          = run_dir,
              silent           = silent,
              normal           = True,
              suffix           = suffix,
              sample_num       = sample_num,
              normal_args        = norm_args)

      current_sd += stddev_step
    # While loop end
