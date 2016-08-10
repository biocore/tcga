#!/bin/env python
# ----------------------------------------------------------------------------
# Copyright (c) 2016--, Gregory Poore.
#
# Distributed under the terms of the Modified BSD License.
#
# ----------------------------------------------------------------------------

import numpy as np
import pandas as pd
import subprocess
import glob

## Call Bash script to transpose data
print "Transposing CDE file..."
subprocess.call("./text_transpose.sh", shell=True)

## Find intersection of headers

lines_list = open('column_names_CDEs.txt').read().splitlines()
split_lines_list = map(str.split, lines_list)
intersection_headers_set = set.intersection(*map(set,split_lines_list))

## Convert set to string

intersection_headers_string = list(intersection_headers_set)
print "The intersection of all CDE columns yields: \n"
print '\n'.join(intersection_headers_string)

## Extract these columns from CDE files into one dataFrame object

print "\n""Concatenating transposed CDE files..."
allTransposedFiles = glob.glob("./Transposed*.txt")
dfIntersect = pd.concat(
	pd.read_csv(eachFile, 
		sep='\t', 
		usecols = intersection_headers_string) 
	for eachFile in allTransposedFiles)

allTransposedFiles = glob.glob("./Transposed*.txt")
dfFullJoin = pd.concat(
	pd.read_csv(eachFile, 
		sep='\t') 
	for eachFile in allTransposedFiles)

## Write to csv file
dfIntersect.to_csv('Intersected_CDEs.csv')
dfFullJoin.to_csv('Full_Joined_CDEs.csv')

print "Complete! Saved in Intersected_CDEs.csv"