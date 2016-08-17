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

class cgcMetadataAPITests(unittest.TestCase):
	""" Tests for cgc_python_api_by_disease_type.py """
	
    def testFoo(self):
        self.failUnless(False)

def main():
    unittest.main()

if __name__ == '__main__':
    main()