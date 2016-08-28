#!/usr/bin/env python

# ----------------------------------------------------------------------------
# Copyright (c) 2016--, Evguenia Kopylova.
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

from parse_kraken_to_biom import (compute_biom_table, prepare_dataframe,
                                  write_biom_table)


class ParseKrakenToBIOMTests(TestCase):
    """Tests for function parse_kraken_to_biom.py
    """

    def setUp(self):
        """Set up working directory and test files
        """
        self.working_dir = mkdtemp()
        # Sample MPA report from Kraken
        self.kraken_translate_report_fp = join(
            self.working_dir, "kraken_mpa_report.txt")
        with open(self.kraken_translate_report_fp, 'w') as tmp:
            tmp.write(kraken_sample_report)
        # Sample MPA report from Kraken unsorted
        self.kraken_translate_report_unsorted_fp = join(
            self.working_dir, "kraken_mpa_report_unsorted.txt")
        with open(self.kraken_translate_report_unsorted_fp, 'w') as tmp:
            tmp.write(kraken_sample_report_unsorted)
        self.taxa_levels = {"domain": "d__",
                            "phylum": "|p__",
                            "class": "|c__",
                            "order": "|o__",
                            "family": "|f__",
                            "genus": "|g__",
                            "species": "|s__"}
        self.taxa_levels_idx = {"d__": 0, "|p__": 1, "|c__": 2,
                                "|o__": 3, "|f__": 4, "|g__": 5,
                                "|s__": 6, "6": "|s__", "5": "|g__",
                                "4": "|f__", "3": "|o__", "2": "|c__",
                                "1": "|p__", "0": "d__"}

    def tearDown(self):
        """Remove tests directory
        """
        rmtree(self.working_dir)

    def test_compute_biom_table(self):
        """Test functionality of compute_biom_table().
        """
        taxonomic_rank = "genus"
        columns = ["s1", "s2", "s3", "s4", "s5"]
        index = ["d__Bacteria|p__Actinobacteria|c__Actinobacteria|o__Actinomycetales|f__Propionibacteriaceae|g__Propionibacterium",
                 "d__Bacteria|p__Firmicutes|c__Bacilli|o__Bacillales|f__Staphylococcaceae|g__Staphylococcus",
                 "d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|o__Enterobacteriales|f__Enterobacteriaceae|g__Escherichia",
                 "d__Bacteria|p__Actinobacteria|c__Actinobacteria|o__Actinomycetales|f__Actinomycetaceae|g__Mobiluncus",
                 "d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|o__Xanthomonadales|f__Xanthomonadaceae|g__Stenotrophomonas",
                 "d__Bacteria|p__Actinobacteria|c__Actinobacteria|o__Actinomycetales|f__Corynebacteriaceae|g__Corynebacterium"]
        exp_table = Table(np.array([[1., 1., 0., 0., 1.],
                                    [1., 0., 0., 0., 0.],
                                    [0., 0., 1., 0., 1.],
                                    [0., 0., 0., 1., 0.],
                                    [0., 0., 0., 1., 0.],
                                    [0., 0., 0., 0., 1.]]),
                          ["k__Bacteria;p__Actinobacteria;c__Actinobacteria;o__Actinomycetales;f__Propionibacteriaceae;g__Propionibacterium",
                           "k__Bacteria;p__Firmicutes;c__Bacilli;o__Bacillales;f__Staphylococcaceae;g__Staphylococcus",
                           "k__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacteriales;f__Enterobacteriaceae;g__Escherichia",
                           "k__Bacteria;p__Actinobacteria;c__Actinobacteria;o__Actinomycetales;f__Actinomycetaceae;g__Mobiluncus",
                           "k__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Xanthomonadales;f__Xanthomonadaceae;g__Stenotrophomonas",
                           "k__Bacteria;p__Actinobacteria;c__Actinobacteria;o__Actinomycetales;f__Corynebacteriaceae;g__Corynebacterium"],
                          columns)
        obs_table = compute_biom_table(
            kraken_translate_report_fp=self.kraken_translate_report_fp,
            taxonomic_rank=taxonomic_rank,
            taxa_levels=self.taxa_levels,
            taxa_levels_idx=self.taxa_levels_idx,
            columns=columns,
            index=index)
        self.assertEqual(obs_table, exp_table)

    def test_compute_biom_table_species(self):
        """Test functionality of compute_biom_table() to species level.
        """
        taxonomic_rank = "species"
        columns = ["s1", "s2", "s3", "s4", "s5"]
        index = ["d__Bacteria|p__Actinobacteria|c__Actinobacteria|o__Actinomycetales|f__Propionibacteriaceae|g__Propionibacterium|s__Propionibacterium_acnes",
                 "d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|o__Enterobacteriales|f__Enterobacteriaceae|g__Escherichia|s__Escherichia_coli",
                 "d__Bacteria|p__Actinobacteria|c__Actinobacteria|o__Actinomycetales|f__Actinomycetaceae|g__Mobiluncus|s__Mobiluncus_curtisii",
                 "d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|o__Xanthomonadales|f__Xanthomonadaceae|g__Stenotrophomonas|s__Stenotrophomonas_maltophilia"]
        exp_table = Table(np.array([[1., 1., 0., 0., 1.],
                                    [0., 0., 1., 0., 1.],
                                    [0., 0., 0., 1., 0.],
                                    [0., 0., 0., 1., 0.]]),
                          ["k__Bacteria;p__Actinobacteria;c__Actinobacteria;o__Actinomycetales;f__Propionibacteriaceae;g__Propionibacterium;s__Propionibacterium_acnes",
                           "k__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacteriales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia_coli",
                           "k__Bacteria;p__Actinobacteria;c__Actinobacteria;o__Actinomycetales;f__Actinomycetaceae;g__Mobiluncus;s__Mobiluncus_curtisii",
                           "k__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Xanthomonadales;f__Xanthomonadaceae;g__Stenotrophomonas;s__Stenotrophomonas_maltophilia"],
                          columns)
        obs_table = compute_biom_table(
            kraken_translate_report_fp=self.kraken_translate_report_fp,
            taxonomic_rank=taxonomic_rank,
            taxa_levels=self.taxa_levels,
            taxa_levels_idx=self.taxa_levels_idx,
            columns=columns,
            index=index)
        self.assertEqual(obs_table, exp_table)

    def test_compute_biom_table_unsorted(self):
        """Test functionality of compute_biom_table() with unsorted input.
        """
        taxonomic_rank = "genus"
        columns = ["s1", "s2", "s5", "s4", "s3"]
        index = ["d__Bacteria|p__Actinobacteria|c__Actinobacteria|o__Actinomycetales|f__Propionibacteriaceae|g__Propionibacterium",
                 "d__Bacteria|p__Firmicutes|c__Bacilli|o__Bacillales|f__Staphylococcaceae|g__Staphylococcus",
                 "d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|o__Enterobacteriales|f__Enterobacteriaceae|g__Escherichia",
                 "d__Bacteria|p__Actinobacteria|c__Actinobacteria|o__Actinomycetales|f__Actinomycetaceae|g__Mobiluncus",
                 "d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|o__Xanthomonadales|f__Xanthomonadaceae|g__Stenotrophomonas",
                 "d__Bacteria|p__Actinobacteria|c__Actinobacteria|o__Actinomycetales|f__Corynebacteriaceae|g__Corynebacterium"]
        exp_table = Table(np.array([[1., 1., 1., 0., 0.],
                                    [1., 0., 0., 0., 0.],
                                    [0., 0., 1., 0., 1.],
                                    [0., 0., 0., 1., 0.],
                                    [0., 0., 0., 1., 0.],
                                    [0., 0., 1., 0., 0.]]),
                          ["k__Bacteria;p__Actinobacteria;c__Actinobacteria;o__Actinomycetales;f__Propionibacteriaceae;g__Propionibacterium",
                           "k__Bacteria;p__Firmicutes;c__Bacilli;o__Bacillales;f__Staphylococcaceae;g__Staphylococcus",
                           "k__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacteriales;f__Enterobacteriaceae;g__Escherichia",
                           "k__Bacteria;p__Actinobacteria;c__Actinobacteria;o__Actinomycetales;f__Actinomycetaceae;g__Mobiluncus",
                           "k__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Xanthomonadales;f__Xanthomonadaceae;g__Stenotrophomonas",
                           "k__Bacteria;p__Actinobacteria;c__Actinobacteria;o__Actinomycetales;f__Corynebacteriaceae;g__Corynebacterium"],
                          columns)
        obs_table = compute_biom_table(
            kraken_translate_report_fp=self.kraken_translate_report_fp,
            taxonomic_rank=taxonomic_rank,
            taxa_levels=self.taxa_levels,
            taxa_levels_idx=self.taxa_levels_idx,
            columns=columns,
            index=index)
        self.assertEqual(obs_table, exp_table)

    def test_prepare_dataframe(self):
        """Test functionality of prepare_dataframe().
        """
        taxonomic_rank = "genus"
        sample_ids_exp = ["s1", "s2", "s3", "s4", "s5"]
        taxonomies_exp = ["d__Bacteria|p__Actinobacteria|c__Actinobacteria|o__Actinomycetales|f__Propionibacteriaceae|g__Propionibacterium",
                          "d__Bacteria|p__Firmicutes|c__Bacilli|o__Bacillales|f__Staphylococcaceae|g__Staphylococcus",
                          "d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|o__Enterobacteriales|f__Enterobacteriaceae|g__Escherichia",
                          "d__Bacteria|p__Actinobacteria|c__Actinobacteria|o__Actinomycetales|f__Actinomycetaceae|g__Mobiluncus",
                          "d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|o__Xanthomonadales|f__Xanthomonadaceae|g__Stenotrophomonas",
                          "d__Bacteria|p__Actinobacteria|c__Actinobacteria|o__Actinomycetales|f__Corynebacteriaceae|g__Corynebacterium"]
        sample_ids_obs, taxonomies_obs = prepare_dataframe(
            kraken_translate_report_fp=self.kraken_translate_report_fp,
            taxonomic_rank=taxonomic_rank,
            taxa_levels=self.taxa_levels,
            taxa_levels_idx=self.taxa_levels_idx)
        self.assertListEqual(sorted(sample_ids_obs), sorted(sample_ids_exp))
        self.assertListEqual(sorted(taxonomies_obs), sorted(taxonomies_exp))

    def test_write_biom_table(self):
        """Test functionality of write_biom_table().
        """
        table_exp = Table(np.array([[1., 1., 1., 0., 0.],
                                    [1., 0., 0., 0., 0.],
                                    [0., 0., 1., 0., 1.],
                                    [0., 0., 0., 1., 0.],
                                    [0., 0., 0., 1., 0.],
                                    [0., 0., 1., 0., 0.]]),
                          ["k__Bacteria;p__Actinobacteria;c__Actinobacteria;o__Actinomycetales;f__Propionibacteriaceae;g__Propionibacterium",
                           "k__Bacteria;p__Firmicutes;c__Bacilli;o__Bacillales;f__Staphylococcaceae;g__Staphylococcus",
                           "k__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacteriales;f__Enterobacteriaceae;g__Escherichia",
                           "k__Bacteria;p__Actinobacteria;c__Actinobacteria;o__Actinomycetales;f__Actinomycetaceae;g__Mobiluncus",
                           "k__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Xanthomonadales;f__Xanthomonadaceae;g__Stenotrophomonas",
                           "k__Bacteria;p__Actinobacteria;c__Actinobacteria;o__Actinomycetales;f__Corynebacteriaceae;g__Corynebacterium"],
                          ["s1", "s2", "s3", "s4", "s5"])
        self.biom_output_fp = join(self.working_dir, "test_output_biom")
        write_biom_table(table_exp, self.biom_output_fp)
        table_obs = load_table(self.biom_output_fp)
        self.assertEqual(table_obs, table_exp)


kraken_sample_report = """s1_0\td__Bacteria|p__Actinobacteria|c__Actinobacteria|o__Actinomycetales|f__Propionibacteriaceae|g__Propionibacterium|s__Propionibacterium_acnes
s1_1\td__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|o__Enterobacteriales|f__Enterobacteriaceae
s1_2\td__Bacteria|p__Firmicutes|c__Bacilli|o__Bacillales|f__Staphylococcaceae|g__Staphylococcus
s2_1\td__Archaea
s2_2\td__Bacteria
s2_3\troot
s2_4\td__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|o__Pasteurellales|f__Pasteurellaceae
s2_5\td__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|o__Enterobacteriales|f__Enterobacteriaceae
s2_6\td__Bacteria|p__Actinobacteria|c__Actinobacteria|o__Actinomycetales|f__Propionibacteriaceae|g__Propionibacterium|s__Propionibacterium_acnes
s3_1\td__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|o__Enterobacteriales|f__Enterobacteriaceae
s3_2\td__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|o__Enterobacteriales|f__Enterobacteriaceae|g__Escherichia|s__Escherichia_coli
s4_1\td__Bacteria|p__Actinobacteria|c__Actinobacteria|o__Actinomycetales|f__Actinomycetaceae|g__Mobiluncus|s__Mobiluncus_curtisii
s4_2\td__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|o__Xanthomonadales|f__Xanthomonadaceae|g__Stenotrophomonas|s__Stenotrophomonas_maltophilia
s5_1\td__Bacteria|p__Actinobacteria|c__Actinobacteria|o__Actinomycetales|f__Propionibacteriaceae|g__Propionibacterium|s__Propionibacterium_acnes
s5_2\td__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|o__Enterobacteriales|f__Enterobacteriaceae|g__Escherichia|s__Escherichia_coli
s5_3\td__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|o__Enterobacteriales|f__Enterobacteriaceae
s5_4\td__Bacteria|p__Actinobacteria|c__Actinobacteria|o__Actinomycetales|f__Corynebacteriaceae|g__Corynebacterium
"""

kraken_sample_report_unsorted = """s1_0\td__Bacteria|p__Actinobacteria|c__Actinobacteria|o__Actinomycetales|f__Propionibacteriaceae|g__Propionibacterium|s__Propionibacterium_acnes
s1_1\td__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|o__Enterobacteriales|f__Enterobacteriaceae
s2_1\td__Archaea
s5_2\td__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|o__Enterobacteriales|f__Enterobacteriaceae|g__Escherichia|s__Escherichia_coli
s2_2\td__Bacteria
s2_3\troot
s2_4\td__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|o__Pasteurellales|f__Pasteurellaceae
s2_5\td__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|o__Enterobacteriales|f__Enterobacteriaceae
s1_2\td__Bacteria|p__Firmicutes|c__Bacilli|o__Bacillales|f__Staphylococcaceae|g__Staphylococcus
s4_2\td__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|o__Xanthomonadales|f__Xanthomonadaceae|g__Stenotrophomonas|s__Stenotrophomonas_maltophilia
s5_4\td__Bacteria|p__Actinobacteria|c__Actinobacteria|o__Actinomycetales|f__Corynebacteriaceae|g__Corynebacterium
s2_6\td__Bacteria|p__Actinobacteria|c__Actinobacteria|o__Actinomycetales|f__Propionibacteriaceae|g__Propionibacterium|s__Propionibacterium_acnes
s3_1\td__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|o__Enterobacteriales|f__Enterobacteriaceae
s3_2\td__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|o__Enterobacteriales|f__Enterobacteriaceae|g__Escherichia|s__Escherichia_coli
s4_1\td__Bacteria|p__Actinobacteria|c__Actinobacteria|o__Actinomycetales|f__Actinomycetaceae|g__Mobiluncus|s__Mobiluncus_curtisii
s5_1\td__Bacteria|p__Actinobacteria|c__Actinobacteria|o__Actinomycetales|f__Propionibacteriaceae|g__Propionibacterium|s__Propionibacterium_acnes
s5_3\td__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|o__Enterobacteriales|f__Enterobacteriaceae
"""

if __name__ == '__main__':
    main()
