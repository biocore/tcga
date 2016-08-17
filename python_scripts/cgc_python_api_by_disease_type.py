#!/usr/bin/env python3
# ----------------------------------------------------------------------------
# Copyright (c) 2016--, Gregory Poore.
#
# Distributed under the terms of the Modified BSD License.
#
# ----------------------------------------------------------------------------
import numpy as np
import pandas as pd
import sevenbridges as sbg
import glob, sys, getopt, os, json, unicodedata, argparse
from pandas.io.json import json_normalize
from tqdm import tqdm # progressbar
from tqdm import trange

parser = argparse.ArgumentParser(description='Grab TCGA metadata using CGC Sevenbridges API.')
parser.add_argument('-o', '--outputFilePath', help = 'Enter output file path (ending in /)', default = './', 
	dest='outputFilePath')
parser.add_argument('-p', '--cgcProjectName', help = 'Enter name of CGC Sevenbridges project', default = 'tcga-microbiome',
	dest='cgcProjectName')
parser.add_argument('-u', '--cgcUserName', help = 'Enter name of CGC Sevenbridges user name', default = 'ROBKNIGHT',
	dest='cgcUserName')
parser.add_argument('-d', '--diseaseTypeList', help = 'Enter desired disease name from which to pull metadata (surround with "")',
 default = '', required=True, nargs='+', dest='diseaseTypeList', choices = ['Cholangiocarcinoma', 'Lymphoid Neoplasm Diffuse Large B-cell Lymphoma', 'Uterine Carcinosarcoma',
'Adrenocortical Carcinoma', 'Mesothelioma', 'Uveal Melanoma', 'Thymoma', 'Kidney Chromophobe', 'Testicular Germ Cell Tumors',
'Pancreatic Adenocarcinoma', 'Pheochromocytoma and Paraganglioma', 'Sarcoma', 'Rectum Adenocarcinoma', 'Glioblastoma Multiforme', 
'Kidney Renal Papillary Cell Carcinoma','Cervical Squamous Cell Carcinoma and Endocervical Adenocarcinoma', 'Acute Myeloid Leukemia',
'Esophageal Carcinoma', 'Liver Hepatocellular Carcinoma', 'Skin Cutaneous Melanoma', 'Bladder Urothelial Carcinoma', 
'Brain Lower Grade Glioma', 'Stomach Adenocarcinoma', 'Prostate Adenocarcinoma', 'Thyroid Carcinoma', 'Lung Squamous Cell Carcinoma',
'Lung Adenocarcinoma', 'Colon Adenocarcinoma', 'Head and Neck Squamous Cell Carcinoma', 'Uterine Corpus Endometrial Carcinoma',
'Kidney Renal Clear Cell Carcinoma', 'Ovarian Serous Cystadenocarcinoma', 'Breast Invasive Carcinoma'])
args = parser.parse_args()
outputFilePath = args.outputFilePath
cgcProjectName = args.cgcProjectName
cgcUserName = args.cgcUserName
diseaseTypeList = args.diseaseTypeList

# #['-o', '--outputFilePath', '-p','--cgcProjectName', '-u', '--cgcUserName', '-d', '--diseaseTypeList']

print("\n")
print("Output directory         :%r" % outputFilePath)
print("CGC project name         :%r" % cgcProjectName)
print("CGC user name            :%r" % cgcUserName)
print("Selected disease type(s) :%r" % diseaseTypeList)

for diseaseTypeLabel in diseaseTypeList:

	### Define disease type to pull metadata from ###
	d_type_label = {'disease_type': diseaseTypeLabel}

	# Use config file for login. See here if you have questions: 
	# https://github.com/sbg/okAPI/blob/aa8a097c6f0be24170b0a0c800460c6defd0d6c9/Recipes/SBPLAT/Setup_API_environment.ipynb
	config_file = sbg.Config(profile='cgc')
	api = sbg.Api(config=config_file)

	print("\n")
	print("Grabbing list of files, their names, and IDs for "+diseaseTypeLabel+"...")
	# Grab list of files, their names and IDs
	file_list = api.files.query(
	    project=cgcUserName+'/'+cgcProjectName,
	    #limit = 1000, # Based on API calling limit of 1000 per hour
	    metadata = d_type_label)
	print("Extracting file names...")
	f_names = [f.name for f in file_list.all()] # Grab filenames
	print("Extracting file IDs...")
	f_id = [f.id for f in file_list.all()] # Grab file IDs

	# Grab and store filenames; convert from Unicode to Python strings
	filenames = []
	for ind, val in enumerate(f_names):
	    filenames.append(unicodedata.normalize('NFKD', f_names[ind]).encode('latin-1','ignore').decode("utf-8"))

	# Store filename strings in DF
	filenamesDF = pd.DataFrame(filenames, columns = ['filename'])

	print("Grabbing files' metadata...")
	# Grab metadata for files of desired disease type
	metaDF = pd.concat( json_normalize(api.files.get(id=f_id[k]).metadata) for k in trange(0,len(f_id)) )

	# Reindex metadata DF to join with filenames DF
	metaDF['index'] = np.arange(len(metaDF))
	metaDF = metaDF.set_index('index')

	# Merge with filenames DF
	joinedMetaDF = pd.concat([filenamesDF, metaDF], axis=1)

	print("Writing metadata to output path...")
	# Write to output path
	outputFileName = 'joinedMetaDF' + '_' + d_type_label["disease_type"].replace(" ","") + '.csv'
	outputFilePathName = outputFilePath + 'joinedMetaDF' + '_' + d_type_label["disease_type"].replace(" ","") + '.csv'
	joinedMetaDF.to_csv(path_or_buf = outputFilePathName)
	print("Complete! Saved as %r in the output directory. \n" % outputFileName)
