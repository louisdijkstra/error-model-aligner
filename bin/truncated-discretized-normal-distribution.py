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

__author__ = "Louis Dijkstra"

usage = """%prog <mean> <std> <outputfile>	
	
	<mean> 		the mean. 
	<std> 		the standard deviation . 
	<outputfile> File with two columns. First column contains 
			the values. The second column the probabilities. 

Creates a file for a truncated discrete approximation of a normal distribution. 
"""

def normalizationFactor (mean, std): 
	return 1.0 / (1.0 - norm.cdf(-mean/std)) 

def main():
	parser = OptionParser(usage=usage)
	parser.add_option("-k", action="store", dest="min_value", default=0, type=int,
				  		help="Minimal value (Default=0)") 
	parser.add_option("-l", action="store", dest="max_value", default=1000, type=int,
				  		help="Maximal value (Default=1000)") 
	parser.add_option("-v", action="store_true", dest="verbose", default=False, 
				  		help="Print also to command line") 
	
	(options, args) = parser.parse_args()
	
	if (len(args)!=3):
		parser.print_help()
		return 1

	mean = float(args[0])
	std = float(args[1])
	
	normalization_factor = normalizationFactor(mean, std)
	
	outputfile = open(args[2], 'w')

	for x in range(options.min_value, options.max_value + 1): 
		probability = normalization_factor * (norm.cdf((x + .5 - mean) / std) - norm.cdf((x - 0.5 - mean)/std))
		outputfile.write("%d\t%0.64f\n"%(x, probability))
		if options.verbose: 
			print("%d\t%0.64f"%(x, probability)) 
	outputfile.close() 

if __name__ == '__main__':
	sys.exit(main())
