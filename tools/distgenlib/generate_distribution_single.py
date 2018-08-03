#!/usr/bin/env python

import sys
import argparse
from dist_generator import dist_generator

##
# \file generate_distribution_single.py
#
# \brief Python command line script which generates probability density and
# cumulative density samples.
#
# \details
# This python script is very similar to harvest_distributions.py. The main
# difference is that this script is designed for exactly one input .dat file.
# whereas harvest_distributions.py is meant to be run on a full run directory.
#
# Usage of this script:
# \code{.bash}
# $ ./generate_distributions_single.py READ_FILE PDF_OUTPUT CDF_OUTPUT
# \endcode
#

# Default sample amount
SAMPLE_NUM = 10000

## Implements the generation
def main():

  # Handle argument parsing
  parser = argparse.ArgumentParser(description='Generate a single distribution')
  parser.add_argument('-w','--bandwidth',
                      metavar='BANDWIDTH',
                      default='0.1',
                      type=float,
                      help='KDE bandwidth factor. Default is 0.1. Not used if'
                      +' distribution is normal.')
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
                        +' 10000.')
  parser.add_argument('readfile', metavar='READFILE', type=str,
                      help='Empirical raw sample files. (usually of type dat)')
  parser.add_argument('pdfout', metavar='PDFOUT', type=str,
                      help='PDF file output location.')
  parser.add_argument('invcdfout', metavar='INVCDFOUT', type=str,
                      help='CDF file output location.')

  args = parser.parse_args()

  read_file_name = args.readfile
  # Probability density function file output
  pdf_file_name = args.pdfout
  # Inverse cumulative density function file output
  invcdf_file_name = args.invcdfout
  bandwidth = args.bandwidth
  offset = args.off
  sample_num = args.samples


  # ---------------
  # Generate the distributions

  generator = dist_generator()
  generator.read_dist(read_file_name)

  if offset != 0.0:
    generator.shift(offset)

  generator.gen_pdf(bandwidth)
  generator.resample(sample_num)

  resampledMax = max(generator.pdf_samples)

  generator.gen_cdf(0,resampledMax,sample_num)
  generator.invert_discrete_cdf(sample_num)

  # ---------------
  # Write files
  generator.write_sample_file(pdf_file_name)
  generator.write_invcdf_file(invcdf_file_name)

  # Quit
  sys.exit(0)

if __name__ == "__main__":
  main()
