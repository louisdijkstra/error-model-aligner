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

__author__ = "Louis Dijkstra"

usage = """%prog <result-dir> 

	<result-dir> 	Folder containing all the data
			for creating the animations and 
			the 2d histogram plots 

The directory <result-dir> should contain three folders: 
	1) non-null/, and
	2) null/, and
	3) animations/ 
	4) 2d-histograms/

Creates both the animations and the 2d histograms. The animations are stored in
animations/ and the histograms in 2d-histograms/. 

See the options if you want to generate only one of these visualizations.

NOTE: should be started from the main directory of the project, so: 

	python bin/%prog.py <result-dir>
"""

def main():
	parser = OptionParser(usage=usage)
	parser.add_option("--animations-only", action="store_true", dest="animations_only", default=False,  
						help="Only animations.")
	parser.add_option("--2d-only", action="store_true", dest="histo_only", default=False,  
						help="Only 2d histograms.")
	(options, args) = parser.parse_args()

	if (len(args)!=1):
		parser.print_help()
		return 1

	result_dir = args[0]
	if result_dir[-1] != '/':
		result_dir += '/'


	# cases where there is a null-distribution! 
	# structure is : <extension>: <description label>
	comparitive_cases = 	{ 
					'.insert-sizes': 'insert size',					  
					'.length-deletion-splits': 'length deletion splits',
					'.length-insertion-splits': 'length insertion splits' 
				}

	# cases WITHOUT a null-distribution! 
	non_comparitive_cases = {
					'.dist-centerpoints-del-splits': 'distance between centerpoints deletion split and variant',
					'.dist-centerpoints-ins-splits': 'distance between centerpoints insertion split and variant',
					'.n-deletion-splits': 'number of deletion splits', 
					'.n-insertion-splits': 'number of insertion splits',
					'.n-insert-size-obs': 'number of insert size observations',
					'.n-obs': 'number of observations',
					'.n-overlapping-obs': 'number of overlapping observations'
				}

	if (not options.animations_only) and (not options.histo_only): # then process both
		options.animations_only = True
		options.histo_only = True 

	
	if options.animations_only: # create the animations
		for extension, description in comparitive_cases.items():
			for variant_type in ['deletion', 'insertion']:
				outputfile = result_dir + 'animations/' + variant_type + extension
				null_case = result_dir + 'null/histogram-data' + extension
				non_null_dir = result_dir + 'non-null/'
				command = ("python bin/create-animation.py -o %s -n %s %s %s %s %s"%(outputfile, null_case,  non_null_dir, variant_type, extension, description)) 
				print("EXECUTING: %s"%command)
				os.system(command)
		for extension, description in comparitive_cases.items():
			for variant_type in ['deletion', 'insertion']:
				outputfile = result_dir + 'animations/' + variant_type + extension
				non_null_dir = result_dir + 'non-null/'
				command = ("python bin/create-animation.py -o %s %s %s %s %s"%(outputfile, non_null_dir,  variant_type, extension, description)) 
				print("EXECUTING: %s"%command)
				os.system(command)	

	if options.histo_only: # create the 2d histograms
		non_comparitive_cases.update(comparitive_cases) # merge them! 
		for extension, description in non_comparitive_cases.items():
			for variant_type in ['deletion', 'insertion']:
				outputfile = result_dir + '2d-histograms/' + variant_type + extension + ".pdf"
				non_null_dir = result_dir + 'non-null/'
				command = ("python bin/create-2d-histogram.py -o %s %s %s %s %s"%(outputfile, non_null_dir, variant_type, extension, description)) 
				print("EXECUTING: %s"%command)
				os.system(command)	



if __name__ == '__main__':
	sys.exit(main())

