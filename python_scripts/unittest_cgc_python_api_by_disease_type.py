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
import glob, sys, getopt, os, json, unicodedata, argparse, unittest, csv, random, subprocess
from pandas.io.json import json_normalize
from tqdm import tqdm # progressbar
from tqdm import trange
from shutil import rmtree
from tempfile import mkdtemp
from os.path import join
# from skbio.util import remove_files

class cgcMetadataAPITests(unittest.TestCase):
    """ Tests for results of cgc_python_api_by_disease_type.py """

    @classmethod
    def setUpClass(self):

        global selected_disease
        """ Set up working directory and test files"""
        cwd = os.getcwd()
        self.working_dir = cwd+'/tmpTestingDirec'
        if not os.path.isdir(self.working_dir):
            os.makedirs(self.working_dir)

        # Select random disease type for testing
        selected_disease = random.choice(disease_choices_lt_1000_files)

        print("\n")
        print("%r disease randomly selected for testing...\n" %selected_disease)

        # Run Python script and place output in temp direc
        subprocCall = 'python3 '+inputFilePath+' -o '+self.working_dir+'/'+' -p '+cgcProjectName+' -u '+cgcUserName+' -d '+'"'+selected_disease+'"'
        print("Using the following call for testing     : %r\n" %subprocCall)
        print("\n---------------------------Running Python Script---------------------------\n")
        subprocess.Popen(subprocCall, shell=True).wait()
        print("\n----------------------------------Complete---------------------------------\n")

    def testForColumnHeaders(self):
        tmpFile = glob.glob(self.working_dir+"/*.csv")[0]
        with open(tmpFile,'r') as f:
            reader = csv.DictReader(f)
            header = reader.fieldnames 
        self.assertCountEqual(header, true_column_headers) # Compares lists for equality; order does not matter  

    def testForTableLength(self):
        tmpFile = glob.glob(self.working_dir+"/*.csv")[0]
        with open(tmpFile,'r') as f:
            reader = csv.reader(f)
            data = list(reader)
            row_count = len(data)
        self.assertEqual( (row_count-1) , table_lengths[selected_disease]) # Subtract 1 from row_count b/c of header 

    @classmethod
    def tearDownClass(self):
        # remove_files(self.files_to_remove)
        rmtree(self.working_dir, ignore_errors=True)

true_column_headers = ['', 
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
    'Breast Invasive Carcinoma' : 1497
}

# Note: API calling limit is 1000 per 5 min
disease_choices_lt_1000_files = ['Cholangiocarcinoma', 'Lymphoid Neoplasm Diffuse Large B-cell Lymphoma', 'Uterine Carcinosarcoma',
'Adrenocortical Carcinoma', 'Mesothelioma', 'Uveal Melanoma', 'Thymoma', 'Kidney Chromophobe', 'Testicular Germ Cell Tumors',
'Pancreatic Adenocarcinoma', 'Pheochromocytoma and Paraganglioma', 'Sarcoma', 'Rectum Adenocarcinoma', 'Glioblastoma Multiforme', 
'Kidney Renal Papillary Cell Carcinoma','Cervical Squamous Cell Carcinoma and Endocervical Adenocarcinoma', 'Acute Myeloid Leukemia',
'Esophageal Carcinoma', 'Liver Hepatocellular Carcinoma', 'Skin Cutaneous Melanoma', 'Bladder Urothelial Carcinoma', 
'Brain Lower Grade Glioma', 'Prostate Adenocarcinoma', 'Thyroid Carcinoma', 'Lung Squamous Cell Carcinoma',
'Lung Adenocarcinoma', 'Head and Neck Squamous Cell Carcinoma']

if __name__ == '__main__':

    global inputFilePath, cgcProjectName, cgcUserName
    cwd = os.getcwd() # For input file
    default_inputFilePath = cwd+'/cgc_python_api_by_disease_type.py'

    parser = argparse.ArgumentParser(description='Testing script options.')
    parser.add_argument('-i', '--inputFilePath', help = 'Enter file path to location of Python script (ending in .py)', 
        default = default_inputFilePath, dest='inputFilePath')
    parser.add_argument('-p', '--cgcProjectName', help = 'Enter name of CGC Sevenbridges project', default = 'tcga-microbiome',
        dest='cgcProjectName')
    parser.add_argument('-u', '--cgcUserName', help = 'Enter name of CGC Sevenbridges user name', default = 'ROBKNIGHT',
        dest='cgcUserName')
    args = parser.parse_args()
    inputFilePath = args.inputFilePath
    cgcProjectName = args.cgcProjectName
    cgcUserName = args.cgcUserName

    print("\n")
    print("CGC project name                            :%r" % cgcProjectName)
    print("CGC user name                               :%r" % cgcUserName)
    print("Using this directory for Python script      :%r" % inputFilePath)

    unittest.main()