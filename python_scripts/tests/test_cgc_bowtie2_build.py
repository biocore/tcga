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

from cgc_bowtie2_build import (load_kraken_mpa_report,
                               parse_repophlan_genome_ids,
                               collect_candidate_repophlan_genomes)


class CGCBowtie2Build(TestCase):
    """Tests for function cgc_bowtie2_build.py
    """

    def setUp(self):
        """Set up working directory and test files
        """
        self.working_dir = mkdtemp()
        # Create files for RepoPhlAn genomes
        self.genomes = ['G001487995', 'G001407635', 'G000306455',
                        'G000343285', 'G001569155', 'G001156185',
                        'G900034955', 'G000529175']
        for genome in self.genomes:
            f = open(join(self.working_dir, "%s.fna" % genome), 'w+')
        # Sample MPA report # 1 from Kraken
        self.kraken_mpa_report_task_1_fp = join(
            self.working_dir, "kraken_mpa_report_task_1.txt")
        with open(self.kraken_mpa_report_task_1_fp, 'w') as tmp:
            tmp.write(kraken_mpa_report_task_1_f)
        # Sample MPA report # 2 from Kraken
        self.kraken_mpa_report_task_2_fp = join(
            self.working_dir, "kraken_mpa_report_task_2.txt")
        with open(self.kraken_mpa_report_task_2_fp, 'w') as tmp:
            tmp.write(kraken_mpa_report_task_2_f)
        # Sample MPA report # 3 from Kraken
        self.kraken_mpa_report_task_3_fp = join(
            self.working_dir, "kraken_mpa_report_task_3.txt")
        with open(self.kraken_mpa_report_task_3_fp, 'w') as tmp:
            tmp.write(kraken_mpa_report_task_3_f)
        # RepoPhlAn genome IDs and taxonomy
        self.repophlan_genome_id_taxonomy_fp = join(
            self.working_dir, "repophlan_genome_id_taxonomy.txt")
        with open(self.repophlan_genome_id_taxonomy_fp, 'w') as tmp:
            tmp.write(repophlan_genome_ids_taxonomy_f)
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

    def test_load_kraken_mpa_report(self):
        """Test functionality of load_kraken_mpa_report()
        """
        kraken_mpa_report_fp = (self.kraken_mpa_report_task_1_fp,
                                self.kraken_mpa_report_task_2_fp,
                                self.kraken_mpa_report_task_3_fp)
        taxonomic_rank = "genus"
        taxa_levels_str = self.taxa_levels[taxonomic_rank]
        taxa_levels_idx_int = self.taxa_levels_idx[taxa_levels_str]
        if taxa_levels_idx_int < 6:
            split_on_level = self.taxa_levels_idx[str(
                taxa_levels_idx_int + 1)]
        else:
            split_on_level = '\t'
        taxonomic_set_obs, split_on_level = load_kraken_mpa_report(
            kraken_mpa_report_fp=kraken_mpa_report_fp,
            read_per_taxa=1,
            taxonomic_rank='genus')
        taxonomic_set_exp = set(['k__Bacteria|p__Proteobacteria|c__Epsilonproteobacteria|o__Campylobacterales|f__Campylobacteraceae|g__Campylobacter',
                                 'k__Bacteria|p__Spirochaetes|c__Spirochaetia|f__Leptospiraceae|g__Leptospira',
                                 'k__Bacteria|p__Firmicutes|c__Bacilli|o__Lactobacillales|f__Streptococcaceae|g__Streptococcus',
                                 'k__Bacteria|p__Firmicutes|c__Bacilli|o__Bacillales|f__Staphylococcaceae|g__Staphylococcus'])
        self.assertEquals(taxonomic_set_obs, taxonomic_set_exp)

    def test_parse_repophlan_genome_ids(self):
        """Test functionality of parse_repophlan_genome_ids()
        """
        repophlan_taxa_and_genomes_exp = {'k__Bacteria|p__Proteobacteria|c__Epsilonproteobacteria|o__Campylobacterales|f__Campylobacteraceae|g__Campylobacter':
                                          ['G001487995', 'G001407635',
                                           'G001493835', 'G001497035',
                                           'G001506965'],
                                          'k__Bacteria|p__Spirochaetes|c__Spirochaetia|f__Leptospiraceae|g__Leptospira':
                                          ['G000347235']}
        repophlan_taxa_and_genomes_obs = parse_repophlan_genome_ids(
            self.repophlan_genome_id_taxonomy_fp,
            split_on_level='|s__')
        self.assertDictEqual(repophlan_taxa_and_genomes_obs,
                             repophlan_taxa_and_genomes_exp)

    def test_collect_candidate_repophlan_genomes(self):
        """Test functionality of collect_candidate_repophlan_genomes()
        """
        # taxonomies from Kraken MPA report
        taxonomic_set = set(['k__Bacteria|p__Proteobacteria|c__Epsilonproteobacteria|o__Campylobacterales|f__Campylobacteraceae|g__Campylobacter',
                             'k__Bacteria|p__Spirochaetes|c__Spirochaetia|f__Leptospiraceae|g__Leptospira',
                             'k__Bacteria|p__Firmicutes|c__Bacilli|o__Lactobacillales|f__Streptococcaceae|g__Streptococcus',
                             'k__Bacteria|p__Firmicutes|c__Bacilli|o__Bacillales|f__Staphylococcaceae|g__Staphylococcus'])
        # repophlan taxonomies with list of associated genomes
        repophlan_taxa_and_genomes = {'k__Bacteria|p__Proteobacteria|c__Epsilonproteobacteria|o__Campylobacterales|f__Campylobacteraceae|g__Campylobacter':
                                      ['G001487995', 'G001407635'],
                                      'k__Bacteria|p__Spirochaetes|c__Spirochaetia|f__Leptospiraceae|g__Leptospira':
                                      ['G000306455', 'G000343285',
                                       'G001569155'],
                                      'k__Bacteria|p__Firmicutes|c__Bacilli|o__Lactobacillales|f__Streptococcaceae|g__Streptococcus':
                                      ['G001156185', 'G900034955'],
                                      'k__Bacteria|p__Firmicutes|c__Bacilli|o__Bacillales|f__Staphylococcaceae|g__Staphylococcus':
                                      ['G000529175']}
        genomes_to_align_to, genome_fps_str =\
            collect_candidate_repophlan_genomes(
                taxonomic_set=taxonomic_set,
                repophlan_taxa_and_genomes=repophlan_taxa_and_genomes,
                genome_dir=self.working_dir)


kraken_mpa_report_task_1_f = """s6980_6263\td__Bacteria|p__Firmicutes|c__Bacilli|o__Lactobacillales|f__Streptococcaceae|g__Streptococcus|s__Streptococcus_pneumoniae
s6980_6264\td__Bacteria|p__Firmicutes|c__Bacilli|o__Bacillales|f__Staphylococcaceae|g__Staphylococcus|s__Staphylococcus_aureus
s6980_6267\td__Bacteria
s6980_6269\td__Bacteria
s6980_6272\td__Bacteria|p__Firmicutes|c__Bacilli|o__Bacillales|f__Staphylococcaceae|g__Staphylococcus|s__Staphylococcus_aureus
s6980_6274\troot
s6980_6278\td__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|o__Alteromonadales|f__Pseudoalteromonadaceae|g__Pseudoalteromonas|s__Pseudoalteromonas_luteoviolacea
s6980_6281\td__Bacteria|p__Firmicutes|c__Bacilli|o__Bacillales|f__Staphylococcaceae|g__Staphylococcus|s__Staphylococcus_aureus
s6980_6285\troot
s6980_6287\td__Bacteria|p__Firmicutes|c__Bacilli|o__Bacillales|f__Staphylococcaceae|g__Staphylococcus|s__Staphylococcus_aureus
s6980_6288\td__Bacteria|p__Proteobacteria|c__Betaproteobacteria|o__Burkholderiales|f__Alcaligenaceae|g__Bordetella|s__Bordetella_bronchiseptica
s6980_6292\td__Bacteria|p__Firmicutes|c__Bacilli|o__Bacillales|f__Staphylococcaceae|g__Staphylococcus|s__Staphylococcus_aureus
s6980_6294\troot
"""

kraken_mpa_report_task_2_f = """s7022_370\td__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|o__Enterobacteriales|f__Enterobacteriaceae|g__Salmonella|s__Salmonella_enterica
s7022_473\td__Bacteria|p__Bacteroidetes|c__Bacteroidia|o__Bacteroidales|f__Porphyromonadaceae|g__Sanguibacteroides|s__Sanguibacteroides_justesenii
s7022_725\td__Bacteria|p__Planctomycetes|c__Planctomycetia|o__Planctomycetales|f__Planctomycetaceae|g__Gemmata|s__Gemmata_obscuriglobus
s7022_732\td__Bacteria|p__Spirochaetes|c__Spirochaetia|f__Leptospiraceae|g__Leptospira|s__Leptospira_interrogans
s7022_733\td__Bacteria|p__Spirochaetes|c__Spirochaetia|f__Leptospiraceae|g__Leptospira|s__Leptospira_interrogans
s7022_746\td__Bacteria
s7022_786\td__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|o__Vibrionales|f__Vibrionaceae|g__Grimontia|s__Grimontia_marina
s7022_789\td__Bacteria|p__Actinobacteria|c__Actinobacteria|o__Streptomycetales|f__Streptomycetaceae|g__Streptomyces|s__Streptomyces_clavuligerus
s7022_790\td__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|o__Enterobacteriales|f__Enterobacteriaceae|g__Escherichia|s__Escherichia_coli
s7022_792\td__Bacteria
s7022_793\td__Bacteria
"""

kraken_mpa_report_task_3_f = """s7076_40737\td__Bacteria|p__Firmicutes|c__Bacilli|o__Bacillales|f__Staphylococcaceae|g__Staphylococcus|s__Staphylococcus_aureus
s7076_40790\td__Bacteria|p__Firmicutes|c__Bacilli|o__Lactobacillales|f__Streptococcaceae|g__Streptococcus|s__Streptococcus_pneumoniae
s7076_40804\td__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|o__Enterobacteriales|f__Enterobacteriaceae|g__Yersinia|s__Yersinia_frederiksenii
s7076_40809\td__Bacteria|p__Proteobacteria|c__Gammaproteobacteria
s7076_40828\td__Bacteria|p__Proteobacteria|c__Gammaproteobacteria
s7076_40829\td__Bacteria|p__Proteobacteria|c__Epsilonproteobacteria|o__Campylobacterales|f__Campylobacteraceae|g__Campylobacter|s__Campylobacter_sp._MIT_97-5078
s7076_40833\td__Bacteria|p__Proteobacteria|c__Epsilonproteobacteria|o__Campylobacterales|f__Campylobacteraceae|g__Campylobacter|s__Campylobacter_sp._MIT_97-5078
s7076_40837\td__Bacteria
s7076_40890\td__Bacteria|p__Firmicutes|c__Bacilli|o__Lactobacillales|f__Streptococcaceae|g__Streptococcus|s__Streptococcus_pneumoniae
s7076_40892\td__Bacteria|p__Firmicutes|c__Bacilli|o__Bacillales|f__Staphylococcaceae|g__Staphylococcus|s__Staphylococcus_aureus
s7076_40893\td__Bacteria
s7076_40894\td__Bacteria
s7076_40895\td__Bacteria
s7076_40896\td__Bacteria
"""

repophlan_genome_ids_taxonomy_f = """G001487995\tk__Bacteria|p__Proteobacteria|c__Epsilonproteobacteria|o__Campylobacterales|f__Campylobacteraceae|g__Campylobacter|s__Campylobacter_coli
G001407635\tk__Bacteria|p__Proteobacteria|c__Epsilonproteobacteria|o__Campylobacterales|f__Campylobacteraceae|g__Campylobacter|s__Campylobacter_jejuni
G001493835\tk__Bacteria|p__Proteobacteria|c__Epsilonproteobacteria|o__Campylobacterales|f__Campylobacteraceae|g__Campylobacter|s__Campylobacter_coli
G001497035\tk__Bacteria|p__Proteobacteria|c__Epsilonproteobacteria|o__Campylobacterales|f__Campylobacteraceae|g__Campylobacter|s__Campylobacter_coli
G001506965\tk__Bacteria|p__Proteobacteria|c__Epsilonproteobacteria|o__Campylobacterales|f__Campylobacteraceae|g__Campylobacter|s__Campylobacter_jejuni
G000347235\tk__Bacteria|p__Spirochaetes|c__Spirochaetia|o__Spirochaetia_noname|f__Leptospiraceae|g__Leptospira|s__Leptospira_kirschneri|t__Leptospira_kirschneri_serovar_Bulgarica_str_Nikolaevo
"""

if __name__ == '__main__':
    main()
