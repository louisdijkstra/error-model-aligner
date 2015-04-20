#!/usr/bin/env python

"""
Copyright (C) 2015 Louis Dijkstra

This file is part of error-model-aligner

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

from __future__ import print_function, division
from optparse import OptionParser
import os
import sys
from scipy.stats import norm
import numpy as np
from scipy.optimize import minimize
import math 

__author__ = "Louis Dijkstra"

usage = """%prog [options] <.insert-sizes>

	<.insert-sizes>	File containing the insert size observations
			when there is no indel
		
Outputs the mean and standard deviation of the null model (i.e., a discrete 
approximation of a Normal distribution that does not allow for negative values)

The file .insert-sizes must be orginazed in two columns (tab seperated): 

	x_1	c_1
	x_2	c_2
	...	...
	x_n	c_n

where x_1 is the minimal insert size observed x_n is the maximum value found. (Note: x_{i+1} = x_i + 1). 
c_i is the count for x_i. 
"""

def normalizationFactor(mu, sigma): 
	"""Returns the normalization factor given mean mu and STD sigma"""
	return 1.0 / (1.0 - norm.cdf((-mu - 0.5)/sigma))


def f(isize, mu, sigma):
	p = norm.cdf((isize + 0.5 - mu)/sigma) - norm.cdf((isize - 0.5 - mu)/sigma)
	if p < sys.float_info.min:
		return sys.float_info.min
	return p

def loglikelihood(mu, sigma, isizes, counts, n): 
	"""Returns the loglikelihood of mu and sigma given the data (isizes, counts and n)"""
	l = n * math.log(normalizationFactor(mu, sigma))
	
	for isize, count in zip(isizes, counts):
		l += count * math.log(f(isize, mu, sigma))
		
	return l

def aux_loglikelihood(var, isizes, counts, n):
	mu = var[0]
	sigma = var[1]
	return -1.0 * loglikelihood(mu, sigma, isizes, counts, n)

def main():
	parser = OptionParser(usage=usage)
	parser.add_option("-f", action="store", dest="maxfun", default=1000, type=int, 
                      		help="Maximum number of function evaluations (Default = 1000) ")
	parser.add_option("-i", action="store", dest="maxiter", default=100, type=int, 
                      		help="Maximum number of iterations (Default = 100) ")
	parser.add_option("-m", action="store", dest="mu_init", default=100.0, type=float, 
                      		help="Initial guess for the mean (mu). (Default is 100) ")
	parser.add_option("-s", action="store", dest="sigma_init", default=10.0, type=float, 
                      		help="Initial guess for the standard deviation (sigma). (Default is 10) ")
	parser.add_option("-v", action="store_true", dest="verbose", default=False,
                      		help="Verbose. Output of the optimizer is printed. ")
	(options, args) = parser.parse_args()

	if (len(args)!=1):
		parser.print_help()
		return 1

	isizes = [] # insert sizes that were observed
	counts = [] # number of times these insert sizes were observed
	for line in open(args[0], 'r'):
		values = map(int, line.split())
		isizes.append(values[0])
		counts.append(values[1])
	isizes 	= np.array(isizes)
	counts 	= np.array(counts)
	n 	= np.sum(counts)

	res = minimize	(	aux_loglikelihood, 
				[options.mu_init, options.sigma_init], 
				args=[isizes, counts, n], 
				method="L-BFGS-B", 
				bounds=[(0, None), (0, None)], 
				options={'disp': options.verbose, 'maxfun': options.maxfun, 'maxiter': options.maxiter})

	print("\n*** RESULTS ***\n")
	print("estimated mean: %lf\t estimated STD: %lf\n"%(res.x[0], res.x[1]))
	print(res.message)

if __name__ == '__main__':
	sys.exit(main())

