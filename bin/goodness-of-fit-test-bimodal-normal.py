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
import numpy as np
from scipy.stats import norm

__author__ = "Louis Dijkstra"

usage = """%prog [options] <null-histogram> <non-null-dir> <mu> <sigma> <epsilon>

	<null-histogram>	Histogram data of the null-case 
	<non-null-dir>	Directory containing all the observations for the 
			non-null cases. 
	<mu>	 	Estimate of the mean 
	<sigma> 	Estimate of the standard deviation
	<epsilon>	Estimate of the "error rate"

Performs a goodness-of-fit test for the model given the three parameter estimates 
and the data.  

See 'estimate-null-insert-sizes.py' for estimating mu and sigma. 
See 'estimate-non-null-insert-sizes.py' for estimating epsilon. 
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
		print('deletion %d: %lf'%(length, l))
	# work through the insertions 
	for length in range(1,1001):
		l += n_ins[length-1] * math.log(normalizationFactor(mu_est, sigma_est, epsilon, -1*length))
		for isize, count in zip(isizes_ins[length-1], counts_ins[length-1]):
			l += count * math.log(epsilon * f(isize, mu_est, sigma_est, 0) + (1.0 - epsilon) * f(isize, mu_est, sigma_est, -1*length))
		print('insertion %d: %lf'%(length, l))
	return -1.0 * l 


def readInFile (filename):
	"""Reads in a histogram file."""
	isizes, counts = [], []
	for line in open(filename, 'r'):
		values = map(int, line.split())
		isizes.append(values[0])
		counts.append(values[1])
	return isizes, counts

def determineTestStatisticNullCase(isizes, counts, mu, sigma):
	""""Computes the chi-squared test statistic for the null case (no indel present)"""
	T = 0.0
	n = float(sum(counts))
	for isize, count in zip(isizes, counts):
		expected_value = n * (normalizationFactor(mu, sigma, 0, 0) * f(isize, mu, sigma, 0)) 
		print(isize, count, expected_value)
		T += (count - expected_value)**2 / expected_value
	return T 

def determineTestStatistic(isizes, counts, is_deletion, length, mu, sigma, epsilon):
	""""Computes the chi-squared test statistic for the non-null case (indel present)"""
	if not is_deletion: 
		length = -1 * length
	T = 0.0
	n = float(sum(counts))
	for isize, count in zip(isizes, counts):
		expected_value = n * (normalizationFactor(mu, sigma, epsilon, length) * (epsilon * f(isize, mu, sigma, 0) + (1.0 - epsilon) * f(isize, mu, sigma, length))) 
		T += (count - expected_value)**2 / expected_value
	return T 

def main():
	parser = OptionParser(usage=usage)
	parser.add_option("-v", action="store_true", dest="verbose", default=False,
                      		help="Verbose.")
	(options, args) = parser.parse_args()

	if (len(args)!=5):
		parser.print_help()
		return 1

	non_null_dir = args[1]
	if non_null_dir[-1] != '/':
		non_null_dir += '/'
	
	mu 	= float(args[2])
	sigma	= float(args[3])
	epsilon	= float(args[4])
	
	# READ IN THE RAW DATA (BOTH DELETIONS and INSERTIONS)

	isizes_null, counts_null = readInFile(args[0])
	n_null = sum(counts_null)
	isizes_null = np.array(isizes_null)
	counts_null = np.array(counts_null)

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
		isizes, counts = readInFile(non_null_dir + 'deletion.length' + str(length) + '.insert-sizes')
		isizes_del.append(isizes)
		counts_del.append(counts)
		n_del.append(sum(counts))
		isizes, counts = readInFile(non_null_dir + 'insertion.length' + str(length) + '.insert-sizes')
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

	T0 = determineTestStatisticNullCase(isizes_null, counts_null, mu, sigma)
	print("Null case: %f"%T0)
	
	#for length in range(1,1001):
	#	T = determineTestStatistic(isizes_del[length-1], counts_del[length-1], True, length, mu, sigma, epsilon)
	#	print("Deletion length %d: %f"%(length, T))
	#for length in range(1,1001):
	#	T = determineTestStatistic(isizes_ins[length-1], counts_ins[length-1], False, length, mu, sigma, epsilon)
#		print("Insertion length %d: %f"%(length, T))


if __name__ == '__main__':
	sys.exit(main())


