#!/usr/bin/env python

# ----------------------------------------------------------------------------
# Copyright (c) 2016--, Evguenia Kopylova, Jad Kanbar
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from unittest import TestCase, main
from shutil import rmtree
from tempfile import mkdtemp
from os.path import join
import numpy as np

from biom.table import Table
from biom import load_table

from cgc_create_tcga_workflow_task import generate_mapping_file


class CreateTCGAWorkflowTask(TestCase):
    """Tests for script cgc_create_tcga_workflow_task.py
    """

    def setUp(self):
        """Set up working directory and test files
        """
        self.working_dir = mkdtemp()
        # Sample QIIME mapping file
        self.qiime_mapping_file_fp = join(
            self.working_dir, "qiime_mapping_file.txt")
        with open(self.qiime_mapping_file_fp, 'w') as tmp:
            tmp.write(qiime_mapping_file)
        disease_type = "Lymphoid Neoplasm Diffuse Large B-cell Lymphoma"
        file_list = list(
            api.files.query(
                project=config['project'],
                metadata = {'disease_type': disease_type}).all())

    def tearDown(self):
        """Remove tests directory
        """
        rmtree(self.working_dir)

    def test_generate_mapping_file(self):
        """Test functionality of generate_mapping_file().
        """
        all_files = []
        output_fp, sampleID_count = generate_mapping_file(
            self.qiime_mapping_file_fp,
            all_files,
            config,
            total_tasks_created,
            output_dp,
            sampleID_count)

qiime_mapping_file = """#SampleID\tBarcodeSequence\tLinkerPrimerSequence\tfile_name\tTreatment\tDOB\tDescription
PC.354\t\tf1.bam\tControl\t20061218\tControl_mouse__I.D._354
PC.355\t\tf2.bam\tControl\t20061218\tControl_mouse__I.D._355
PC.356\t\tf3.bam\tControl\t20061126\tControl_mouse__I.D._356
PC.481\t\tf4.bam\tControl\t20070314\tControl_mouse__I.D._481
PC.593\t\tf5.bam\tControl\t20071210\tControl_mouse__I.D._593
PC.607\t\tf6.bam\tFast\t20071112\tFasting_mouse__I.D._607
PC.634\t\tf7.bam\tFast\t20080116\tFasting_mouse__I.D._634
PC.635\t\tf8.bam\tFast\t20080116\tFasting_mouse__I.D._635
PC.636\t\tf9.bam\tFast\t20080116\tFasting_mouse__I.D._636
"""


if __name__ == '__main__':
    main()