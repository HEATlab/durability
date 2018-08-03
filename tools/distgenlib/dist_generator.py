#!/usr/bin/env python

from scipy import stats # Note, this script requires scipy to run
from scipy.stats import norm

##
# \file dist_generator.py
#
# \brief Library class for generating PDFs and CDFs from sample data.
#
# \details
# Example usage in python scripts:
#
# \code{.py}
# #...
# gen = dist_generator()
# gen.read_dist(file_name)
# gen.gen_pdf(kde_factor=0.1)
# gen.resample(resample_num=10000)
# dist = gen.pdf_samples
#
# gen.gen_cdf(0,max(dist),sample_num=10000)
# cdf = gen.cdf_samples
#
# gen.invert_discrete_cdf(resample_num=10000)
# invcdf = gen.invcdf_samples
# #...
# \endcode
#
# See harvest_distributions.py for more details on usage.

class dist_generator:
  def __init__(self):
    ## Collected distribution from read files.
    self.dist = []
    ## Generated Probability Density Function
    self.pdf = None
    self.pdf_type = '' # String indicating the type of PDF.
    self.pdf_samples = []
    self.cdf_samples = []
    self.invcdf_samples = []


  ##
  # \fn read_dist
  #
  # \brief Reads a single column file of floats, and converts it to a list.
  #
  # @param read_file_name Name of file to read
  # @return A python List of all the entries in the file
  def read_dist(self, read_file_name):
    read_file = open(read_file_name,'r')
    self.dist = [float(line.strip()) for line in read_file if line != "\n"]
    read_file.close()


  ##
  # \fn shift
  #
  # \brief Adds a set offset to each of the underlying distribution data points.
  #
  # @param offset Double indicating how much we should shift the data.
  # @post Underlying distribution is modified. Invalidates generated PDFs and
  # CDFs.
  def shift(self, offset):
    for i in range(len(self.dist)):
      self.dist[i] += offset


  ##
  # \fn gen_pdf
  #
  # \brief Generate a PDF from a sample file
  #
  # @param kde_factor Bandwidth factor of the Gaussian Kernel Density model.
  # @post assigns the pdf property to a gaussian_kde object generated.
  def gen_pdf(self, kde_factor):
    self.pdf = stats.gaussian_kde(self.dist, kde_factor)
    self.pdf_type = 'kde'


  ##
  # \fn gen_normal
  #
  # \brief Generate a normal distribution that approximates the given
  #   underlying, read-in distribution.
  #
  # @param mean Mean of the normal distribution. If None, will take
  #   empirical distribution data.
  # @param sd Standard deviation of the normal distribution. If None, will
  #   take empirical distribution data.
  # @post Internally generates the PDF, and sets pdf_type to 'norm'.
  def gen_normal(self, mean=None,sd=None):
    # If mean or sd aren't assigned, generate the normal distribution from
    # the raw distribution.
    if mean == None:
      # Set mean
      mean = self.distmean()
    if sd == None:
      # Set standard deviation
      sd = self.distsd()
    self.pdf = {'mean' : mean, 'sd' : sd}
    self.pdf_type = 'norm'


  ## Return the mean of the raw distribution
  def distmean(self):
    return float(reduce(lambda x,y: x+y, self.dist))/len(self.dist)


  ## Return the standard deviation of the raw distribution
  def distsd(self):
    mean = self.distmean()
    n = len(self.dist)
    elems = list(map(lambda i: ((i-mean)**2), self.dist))
    sd = (float(reduce(lambda x,y: x+y, elems))/n)**0.5
    return sd


  ##
  # \fn resample
  #
  # \brief Generate a number of samples
  #
  # @param resample_num
  # @post assigns the pdf_samples property with a resampling from the internal
  #   pdf with resample_num elements
  def resample(self, resample_num):
    if self.pdf_type == 'kde':
      # Sample from the KDE object
      self.pdf_samples = self.pdf.resample(resample_num)[0]
    elif self.pdf_type == 'norm':
      # Randomly sample from Normal distribution
      self.pdf_samples = norm.rvs(loc=self.pdf['mean'], scale=self.pdf['sd'],
        size=resample_num)

    # Prevent any samples from being less than 0
    for i in range(len(self.pdf_samples)):
      if self.pdf_samples[i] < 0.0:
        self.pdf_samples[i] = 0


  ##
  # \fn gen_cdf
  #
  # \brief Internally generate a cumulative distribution function.
  #
  # @param pdf 1-Dim Probability distribution function to integrate.
  # @param sample_num Number of discrete floating point outputs
  # @post Sets cdf_samples to a list of lists with the format
  #   [times, cdf values], where the cdf values represent a cumulative
  #   distribution function.
  def gen_cdf(self, lower_lim, upper_lim, sample_num):

    if self.pdf == None:
      raise StandardError("No PDF model found. Did gen_pdf() or gen_normal()"
        +" run?")

    # Evaluate step size
    step_size = (upper_lim - lower_lim) / float(sample_num)
    # Get the time points
    times = [ step_size*count + lower_lim for count in range(sample_num) ]

    # --------------------------------------------------------------------------

    if self.pdf_type == 'norm': # FOR NORMAL
      cdflamb = lambda x: norm.cdf(x, self.pdf['mean'], self.pdf['sd'])
      cdfvals = [cdflamb(t) for t in times]

    elif self.pdf_type == 'kde': # FOR KDE
      # Use a list comprehension to integrate over the pdf to create a cdf
      cdfvals = [ self.pdf.integrate_box_1d(lower_lim, tPoint)
        for tPoint in times ]

      # Normalise the cdf so that the full integral is a cummulative probability
      # of 1.0

      # How far beyond the upper bound do we consider?
      f_integral_factor = 2

      full_integral = self.pdf.integrate_box_1d(lower_lim,
        f_integral_factor*upper_lim)
      cdfvals = [float(cdfval)/full_integral for cdfval in cdfvals]
    else:
      raise StandardError("Internal PDF Type Not Found: " + str(self.pdf_type))

    # --------------------------------------------------------------------------

    # Actual output in the form of [times, cdfvals] for easy evaluation
    self.cdf_samples = [times,cdfvals]


  ##
  # \fn invert_discrete_cdf
  #
  # \brief Invert a set of discrete CDF samples using linear interpolation.
  #
  # This function reads in discrete samples from a CDF, inverts them, and then
  # resamples the inversion using linear interpolation between the points
  # such that the resampling on the probability axis is evenly spaced.
  #
  # @param resample_num Amount of new, evenly spaced samples to take.
  # @return A list representing discrete points in an inverse CDF
  def invert_discrete_cdf(self, resample_num):
    if (self.cdf_samples != []):
      # Shorthand variable
      cdf_samples = self.cdf_samples
    else:
      raise StandardError("No CDF samples found, nothing to invert.")
      sys.exit(1)
    # Calculate the step
    step_size = (cdf_samples[1][-1] - cdf_samples[1][0])/resample_num
    # Initialise the invCDF list
    invCDF = [[0],[0]]

    sample_count = 1 # Samples we've taken so far
    next_sample = step_size # Next sample x value

    # Iterate over every element in the cdf_samples
    i = 1 # starting index
    last_threshold_index = 0
    while i < len(cdf_samples[1]):

      if cdf_samples[1][i] >= next_sample:
        # We've exceeded the next_sample threshold, now interpolate

        slope = (cdf_samples[0][i] - cdf_samples[0][i-1]) \
          / float(cdf_samples[1][i] - cdf_samples[1][i-1])


        invCDF[0].append(next_sample)
        invCDF[1].append(slope*(cdf_samples[1][i] - cdf_samples[1][i-1]) \
          + cdf_samples[0][i-1])

        # We've taken a new sample, increment the count
        sample_count += 1
        next_sample += step_size

        # Push back i to the last time we exceeded the next_sample
        temp_index = i # temporary variable to store i
        i = last_threshold_index # push back i
        last_threshold_index = temp_index # mark current index as exceeding

      else:
        i += 1

    self.invcdf_samples = invCDF

    return invCDF


  ##
  # \fn write_sample_file
  #
  # \brief Writes a sample file given from the generated resampling
  #
  # @param write_file_name The write file name.
  def write_sample_file(self, write_file_name):

    write_file = open(write_file_name, 'w')
    self.write_list(write_file, self.pdf_samples)
    write_file.close()


  ##
  # \fn write_invCDF_file
  #
  # \brief Writes an inverse CDF file for a given generated CDF
  #
  # @param inv_cdf_file_name The write file name
  def write_invcdf_file(self, inv_cdf_file_name):
    # Inverse Cumulative Density file to write to
    inv_cdf_file = open(inv_cdf_file_name,'w')
    zipped_list = zip(self.invcdf_samples[0],self.invcdf_samples[1])
    # Integrate the pdf and write to file
    self.write_list(inv_cdf_file, zipped_list, True)
    inv_cdf_file.close()


  ##
  # \fn write_list
  #
  # \brief Utility function. Prints a list to a given file as a column
  #
  # @param write_file File object to write to
  # @param l List to write
  # @param pairs boolean indicating whether this is a nested list/tuple.
  def write_list(self, write_file, l, pairs=False):
    for entry in l:
      if not pairs:
        write_file.write(str(entry)+'\n')
      else:
        write_file.write(str(entry[0])+','+str(entry[1])+'\n')


# end class dist_generator
