#!/bin/env python
# ----------------------------------------------------------------------------
# Copyright (c) 2016--, Gregory Poore.
#
# Distributed under the terms of the Modified BSD License.
#
# ----------------------------------------------------------------------------
''' Tips:
- inputMetadataFilePath should be the FULL PATH to "Combined_CDEs" folder that "Pull_metadata_from_Firehose.sh" creates
- outputFilePath should also be given as the FULL PATH to its directory
- Make sure to end these paths with / to ensure that files are read and written correctly
- The data transformation (i.e. transpose) is dependent on GNU datamash (http://www.gnu.org/software/datamash/)
- The progress bar requires the tqdm dependency (https://github.com/tqdm/tqdm); install with "pip install tqdm"
'''
import numpy as np
import pandas as pd
import subprocess, glob, sys, getopt, os
from tqdm import tqdm # progressbar

## Defining getopts module

# args = '-i -o'.split()
# optlist, args = getopt.getopt(args, )

def main(argv):
   global inputMetadataFilePath, outputFilePath
   inputMetadataFilePath = './'
   outputFilePath = './'
   try:
      opts, args = getopt.getopt(argv,"hi:o:",["inputMetadataFilePath=","outputFilePath"])
   except getopt.GetoptError:
      print 'make_merged_CDE_table.py -i <inputMetadataFilePath> -o <outputFilePath>'
      sys.exit(2)
   for opt, arg in opts:
      if opt == '-h':
         print 'make_merged_CDE_table.py -i <inputMetadataFilePath> -o <outputFilePath>'
         sys.exit()
      elif opt in ("-i", "--inputMetadataFilePath"):
         inputMetadataFilePath = arg
      elif opt in ("-o", "--outputFilePath"):
         outputFilePath = arg
   print "Input file path to metadata is %r" % inputMetadataFilePath
   print "Output file is %r" % outputFilePath

if __name__ == "__main__":
   main(sys.argv[1:])

## Run Bash commands to transpose data
print "Transposing CDE file..."

# First bash command --> data transposing
#--------------------------------------------#
bashArgTransposeBase = """ for f in All_CDEs*.txt
do
    echo "Transposing $f file... Saved as Transposed_$f"
    cat $f | datamash transpose > %rTransposed_$f
done
"""
bashArgTranspose = bashArgTransposeBase % outputFilePath
subprocess.Popen(bashArgTranspose, cwd=inputMetadataFilePath, shell=True, stdout = subprocess.PIPE, stderr = subprocess.PIPE).wait()
#--------------------------------------------#
print "Transposing complete!"
print "Now extracting column headers and saving in column_names_CDEs.txt ...", # The comma here is intentional to print on the same line

# Second bash command --> header extraction for downstream intersection
#--------------------------------------------#
bashArgHeaderExtraction = """ grep '^bcr_patient_barcode' -h Transposed_All_CDEs*.txt > column_names_CDEs.txt """
subprocess.Popen(bashArgHeaderExtraction, cwd=outputFilePath, shell=True, stdout = subprocess.PIPE, stderr = subprocess.PIPE).wait()
#--------------------------------------------#

print "Success!" # Note that this is linked to the previous print command to print on the same line

## Find intersection of headers

print "Now finding intersection of column names ...", # comma intentional
colFileName = 'column_names_CDEs.txt'
lines_list = open(outputFilePath+colFileName).read().splitlines()
split_lines_list = map(str.split, lines_list)
intersection_headers_set = set.intersection(*map(set,split_lines_list))
bashArgRemoveColFile = """ rm -f column_names_CDEs.txt """
subprocess.Popen(bashArgRemoveColFile, cwd=outputFilePath, shell=True, stdout = subprocess.PIPE, stderr = subprocess.PIPE).wait()
print "Complete!\n" # linked to previous print command

## Convert set to string

intersection_headers_string = list(intersection_headers_set)
print "The intersection of all CDE columns yields: \n"
print '\n'.join(intersection_headers_string)

## Extract these columns from CDE files into one dataFrame object

allTransposedFiles = glob.glob(os.path.join(outputFilePath,r'Transposed*.txt'))

print "\nConcatenating intersected CDE fields:"
dfIntersect = pd.concat( pd.read_csv(eachFile, sep='\t', usecols = intersection_headers_string) 
    for eachFile in tqdm(allTransposedFiles)) # tqdm gives a progress bar
print "Concatenating fully-joined CDE fields:"
dfFullJoin = pd.concat( pd.read_csv(eachFile, sep='\t') 
    for eachFile in tqdm(allTransposedFiles)) # tqdm gives a progress bar

## Write to csv file
print "Writing csv files..."
dfIntersect.to_csv(os.path.join(outputFilePath,r'Intersected_CDEs.csv'))
dfFullJoin.to_csv(os.path.join(outputFilePath,r'Full_Joined_CDEs.csv'))
print "Complete! Files saved as Intersected_CDEs.csv and Full_Joined_CDEs.csv in directory: %r\n" % outputFilePath