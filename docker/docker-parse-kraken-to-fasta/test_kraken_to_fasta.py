#!/usr/bin/env python

# ----------------------------------------------------------------------------
# Copyright (c) 2017--, Evguenia Kopylova.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from unittest import TestCase, main
from shutil import rmtree
from tempfile import mkdtemp
from os.path import join


from kraken_to_fasta import get_kraken_report


class ParseKrakenToFastaTest(TestCase):
    """Tests for function kraken_to_fasta.py
    """

    def setUp(self):
        """Set up working directory and test files
        """
        self.working_dir = mkdtemp()
        # Sample Kraken report
        self.kraken_report_fp = join(
            self.working_dir, "kraken_report.txt")
        with open(self.kraken_report_fp, 'w') as tmp:
            tmp.write(kraken_sample_report)
        # Sample FASTA file
        self.sample_fasta_fp = join(
            self.working_dir, "sample.fasta")
        with open(self.sample_fasta_fp, 'w') as tmp:
            tmp.write(sample_fasta)

    def tearDown(self):
        """Remove tests directory
        """
        rmtree(self.working_dir)

    def test_get_kraken_report_classied(self):
        """Test functionality of get_kraken_report(), classified seqs
        """
        report_exp = set(['s6281_3', 's6281_57'])
        report_obs = get_kraken_report(kraken_report_fp=self.kraken_report_fp,
                                       classification_id='C')
        self.assertEqual(report_obs, report_exp)

    def test_get_kraken_report_unclassified(self):
        """Test functionality of get_kraken_report(), unclassified seqs
        """
        report_exp = set(['s6281_98058', 's6281_98059', 's6281_98066',
                          's6281_98067', 's6281_98073'])
        report_obs = get_kraken_report(kraken_report_fp=self.kraken_report_fp,
                                       classification_id='U')
        self.assertEqual(report_obs, report_exp)


kraken_sample_report = """U\ts6281_98058\t0\t50\t0:20
U\ts6281_98059\t0\t50\t0:20
U\ts6281_98066\t0\t50\t0:20
C\ts6281_3\t2\t50\t0:3 2:17
U\ts6281_98067\t0\t50\t0:20
U\ts6281_98073\t0\t50\t0:20
C\ts6281_57\t1773\t50\t0:1 1773:13 0:6
"""

sample_fasta = """>s6281_98059
GTACGTAGCTAGTCGCTGCTCGCTGATCGTAGCTCGCGCTAGCTAGTCGGCTAGTCGTAGCCGT
>s6281_98066
CGTCGCTGCTGCTAGCTGCTCGCCAAAAATGGCTGATCGTCGCGCTAGCTACGTCGTACGTCGC
>s6281_3
GTCGCTAGCTCGCTGATCGTCGTACGTACGACTGCTCGTACGTACGTACGTACGTAGCTCGTGA
>s6281_98067
CCCGATCGTCGCTGCTCGTACGTACGTCGTACGATCGTACGATCGATCGCTGACTACGTACGTC
>s6281_98073
GTCGTAGCTAGCTAGCTACGTCGTCGTACGTCGTACGTACGTCGTACGTACGTACGTCGCTAGC
>s6281_57
GTCGATGCTACGTCGCCGTACGAAATCGATGCTACGTCGTAGCTAGCTCGTCGATCGTACGTAG
"""


if __name__ == '__main__':
    main()
