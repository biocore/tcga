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
import glob, sys, getopt, os, json, unicodedata, argparse, unittest
from pandas.io.json import json_normalize
from tqdm import tqdm # progressbar
from tqdm import trange
from shutil import rmtree
from tempfile import mkdtemp
from os.path import join
from skbio.util import remove_files

parser = argparse.ArgumentParser(description='Test results of TCGA metadata pulled from CGC Sevenbridges API.')
parser.add_argument('-i', '--inputFilePath', 
	help = 'Enter directory of previously created .csv file(s); note, this ends in /', 
	dest='inputFilePath', 
	required=True, 
	nargs='+',)
args = parser.parse_args()
inputFilePath = args.inputFilePath

class cgcMetadataAPITests(unittest.TestCase):
	""" Tests for cgc_python_api_by_disease_type.py """

	def setUp(self):
		""" Set up working directory and test files"""
		
		# Test output is written to this temporary directory
		self.working_dir = mkdtemp()

		# Test data output
		self.cgc_output_by_disease = join(
			self.working_dir, "")




    def testFoo(self):
        self.failUnless(False)

def main():
    unittest.main()

column_headers = ['', 
'filename',
'age_at_diagnosis',
'aliquot_id',
'aliquot_uuid',
'case_id',
'case_uuid',
'data_format',
'data_subtype',
'data_type',
'days_to_death',
'disease_type',
'ethnicity',
'experimental_strategy',
'gender',
'investigation',
'platform',
'primary_site',
'race',
'reference_genome',
'sample_id',
'sample_type',
'sample_uuid',
'vital_status'
]

table_lengths = {
	'Cholangiocarcinoma' : 45,
	'Uterine Carcinosarcoma' : 57,
	'Lymphoid Neoplasm Diffuse Large B-cell Lymphoma' : 62,
	'Adrenocortical Carcinoma' : 79,
	'Mesothelioma' : 87,
	'Thymoma' : 122,
	'Testicular Germ Cell Tumors' : 156,
	'Uveal Melanoma' : 182,
	'Pancreatic Adenocarcinoma' : 183,
	'Pheochromocytoma and Paraganglioma' : 187,
	'Kidney Chromophobe' : 191,
	'Acute Myeloid Leukemia' : 335,
	'Esophageal Carcinoma' : 340,
	'Sarcoma' : 347,
	'Rectum Adenocarcinoma' : 375,
	'Kidney Renal Papillary Cell Carcinoma' : 400,
	'Cervical Squamous Cell Carcinoma and Endocervical Adenocarcinoma' : 451,
	'Glioblastoma Multiforme' : 520,
	'Liver Hepatocellular Carcinoma' : 532,
	'Lung Squamous Cell Carcinoma' : 653,
	'Bladder Urothelial Carcinoma' : 729,
	'Brain Lower Grade Glioma' : 733,
	'Skin Cutaneous Melanoma' : 793,
	'Prostate Adenocarcinoma' : 830,
	'Thyroid Carcinoma' : 880,
	'Head and Neck Squamous Cell Carcinoma' : 907,
	'Lung Adenocarcinoma' : 951,
	'Colon Adenocarcinoma' : 1018,
	'Ovarian Serous Cystadenocarcinoma' : 1039,
	'Stomach Adenocarcinoma' : 1092,
	'Kidney Renal Clear Cell Carcinoma' : 1142,
	'Uterine Corpus Endometrial Carcinoma' : 1239,
	'Breast Invasive Carcinoma' : 1497,
}

if __name__ == '__main__':
    main()