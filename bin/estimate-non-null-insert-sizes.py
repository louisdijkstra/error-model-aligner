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

usage = """%prog [options] <result-dir> <mu-est> <sigma-est>

	<result-dir>	Directory containing all the observations for the 
			non-null cases. 
	<mu-est> 	Estimate of the mean 
	<sigma-est> 	Estimate of the standard deviation

Outputs an estimate of the 'error rate', epsilon, for the non-null arm of 
the model. 

See 'estimate-null-insert-sizes.py' for estimating mu and sigma. 
"""

def normalizationFactor(mu, sigma, epsilon, length): 
	"""Returns the normalization factor given an estimate for the mean mu and STD sigma,
	   the error rate epsilon and the length of the indel involved."""
	return 1.0 / (1.0 - (1.0 - epsilon) * norm.cdf((-mu - 0.5 - length) / sigma) - epsilon * norm.cdf((-mu - 0.5)/sigma))

def f(isize, mu, sigma, length):
	p = norm.cdf((isize + 0.5 - mu - length)/sigma) - norm.cdf((isize - 0.5 - mu - length)/sigma)
	if p < sys.float_info.min:
		return sys.float_info.min
	return p

def negative_loglikelihood(epsilon, mu_est, sigma_est, isizes_del, counts_del, n_del, isizes_ins, counts_ins, n_ins):
	"""returns -1 * loglikelihood"""
	l = 0.0
	# work through the deletions 
	for length in range(1,1001):
		l += n_del[length-1] * math.log(normalizationFactor(mu_est, sigma_est, epsilon, length))
		for isize, count in zip(isizes_del[length-1], counts_del[length-1]):
			l += count * math.log(epsilon * f(isize, mu_est, sigma_est, 0) + (1.0 - epsilon) * f(isize, mu_est, sigma_est, length))
		# print('updated log-likelihood for observations induced by deletion of length %d: %lf'%(length, l))
	# work through the insertions 
	for length in range(1,1001):
		l += n_ins[length-1] * math.log(normalizationFactor(mu_est, sigma_est, epsilon, -1*length))
		for isize, count in zip(isizes_ins[length-1], counts_ins[length-1]):
			l += count * math.log(epsilon * f(isize, mu_est, sigma_est, 0) + (1.0 - epsilon) * f(isize, mu_est, sigma_est, -1*length))
		# print('updated log-likelihood for observations induced by insertions of length %d: %lf'%(length, l))
	return -1.0 * l 


def readInFile (filename):
	"""Reads in a histogram file."""
	isizes, counts = [], []
	for line in open(filename, 'r'):
		values = map(int, line.split())
		isizes.append(values[0])
		counts.append(values[1])
	return isizes, counts

def main():
	parser = OptionParser(usage=usage)
	parser.add_option("-f", action="store", dest="maxfun", default=1000, type=int, 
                      		help="Maximum number of function evaluations (Default = 1000) ")
	parser.add_option("-i", action="store", dest="maxiter", default=100, type=int, 
                      		help="Maximum number of iterations (Default = 100) ")
	parser.add_option("-k", action="store", dest="epsilon_min", default=0.001, type=float, 
                      		help="Lower bound for epsilon. (Default is 0.001) ")
	parser.add_option("-l", action="store", dest="epsilon_init", default=0.05, type=float, 
                      		help="Initial guess for epsilon. (Default is 0.05) ")
	parser.add_option("-m", action="store", dest="epsilon_max", default=0.10, type=float, 
                      		help="Upper bound for epsilon. (Default is 0.10) ")
	parser.add_option("-v", action="store_true", dest="verbose", default=False,
                      		help="Verbose. Output of the optimizer is printed. ")
	(options, args) = parser.parse_args()

	if (len(args)!=3):
		parser.print_help()
		return 1

	result_dir = args[0]
	if result_dir[-1] != '/':
		result_dir += '/'
	
	mu_est 		= float(args[1])
	sigma_est 	= float(args[2])
	
	# READ IN THE RAW DATA (BOTH DELETIONS and INSERTIONS)

	# every variant type and every length will be represented by two lists: 
	# insert sizes (isizes) and # of observations for that insert size (counts)
	# n - # total number of observations 
	isizes_del, isizes_ins 	= [], [] # insert sizes that were observed
	counts_del, counts_ins 	= [], [] # number of times these insert sizes were observed
	n_del, n_ins      	= [], [] # total number of observations 

	# walk through all the files
	for length in range(1,1001): 
		if options.verbose: 
			print("Reading data for indels of length %d"%length) 
		isizes, counts = readInFile(result_dir + 'deletion.length' + str(length) + '.insert-sizes')
		isizes_del.append(isizes)
		counts_del.append(counts)
		n_del.append(sum(counts))
		isizes, counts = readInFile(result_dir + 'insertion.length' + str(length) + '.insert-sizes')
		isizes_ins.append(isizes)
		counts_ins.append(counts)
		n_ins.append(sum(counts))
		
	isizes_del 	= np.array(isizes_del)
	counts_del	= np.array(counts_del)
	n_del	 	= np.array(n_del)
	isizes_ins 	= np.array(isizes_ins)
	counts_ins	= np.array(counts_ins)
	n_ins	 	= np.array(n_ins)
	if options.verbose:
		print("DONE Reading in data...")
	
	res = minimize	(	negative_loglikelihood, 
				options.epsilon_init, 
				args=[mu_est, sigma_est, isizes_del, counts_del, n_del, isizes_ins, counts_ins, n_ins], 
				method="L-BFGS-B", 
				bounds=[(options.epsilon_min, options.epsilon_max)], 
				options={'disp': options.verbose, 'maxfun': options.maxfun, 'maxiter': options.maxiter})
	

	print("\n*** RESULTS ***\n")
	print("estimated epsilon ", res.x)
	print(res.message)

if __name__ == '__main__':
	sys.exit(main())

