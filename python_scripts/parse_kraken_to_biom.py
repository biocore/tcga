#!/usr/bin/env python

#-----------------------------------------------------------------------------
# Copyright (c) 2016--, Evguenia Kopylova.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

"""
Parse output of kraken-mpa-report and generate single BIOM table.
"""

import sys
import os
from os.path import isfile, join
import random
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
# requires matplot version 1.2.1 (or one that supports fontsize in the legend)
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import brewer2mpl


def compute_abundances_biom(reports_dp,
                            metadata_fp,
                            taxonomy_report_dp,
                            taxa_level,
                            taxa_levels,
                            taxa_levels_idx):
    """Absolute abundance of number of reads matching a defined taxa level.

    Parameters
    ----------
    reports_dp: string
        output directory path
    metadata_fp: string
        path to TCGA tsv metadata file
    taxonomy_report_dp: string
        path to directory with output(s) of kraken-mpa-report
    taxa_level: string
        taxonomy level (e.g., genus or species)
    taxa_levels: dictionary
        "domain": "d__",
        "phylum": "|p__",
        "class": "|c__",
        "order": "|o__",
        "family": "|f__",
        "genus": "|g__",
        "species": "|s__"
    taxa_levels_idx: dictionary
        "d__": 0, "|p__": 1, "|c__": 2,
        "|o__": 3, "|f__": 4, "|g__": 5,
        "|s__": 6, "6": "|s__", "5": "|g__",
        "4": "|f__", "3": "|o__", "2": "|c__",
        "1": "|p__", "0": "d__"

    Returns
    -------
    
    """
    total_levels = len(taxa_levels)

    # keys are file names
    metadata = {}
    num_samples = 1
    with open(metadata_fp, 'U') as metadata_f:
        for line in metadata_f:
            # skip the header
            if line.startswith('#'):
                continue
            line = line.strip().split("\t")
            file_n = line[6].split(".bam")[0]
            if file_n not in metadata:
                metadata[file_n] = [line[7], line[5], line[9], line[0]]
                num_samples += 1

    organs_list = sorted({metadata[file_n][1] for file_n in metadata if isfile(join(taxonomy_report_dp, "%s_kraken.report_mpa" % file_n))})
    print organs_list

    # species abundance
    abundances = {}
    labels = []
    total_organs = len(organs_list)
    taxa_l = taxa_levels[taxa_level]
    sample_id = 0
    for file_n in metadata:
        kraken_report_fp = join(taxonomy_report_dp, "%s_kraken.report_mpa" % file_n) 
        if isfile(kraken_report_fp):
            # record abundances
            with open(kraken_report_fp, 'U') as kraken_report_f:
                for line in kraken_report_f:
                    stop_search = False
                    #print "%s\t%s" % (taxa_l, line)
                    # get statistics only for specified taxa level
                    if taxa_l in line:
                        # cannot have deeper levels after specified level
                        for i in range(taxa_levels_idx[taxa_l]+1, total_levels):
                            if taxa_levels_idx[str(i)] in line:
                                stop_search = True
                                break
                    else:
                        stop_search = True
                    if stop_search:
                        continue
                    line = line.strip().split()
                    taxa = line[0]
                    abundance = np.float128(line[1])
                    # check if taxa exists
                    if taxa not in abundances:
                        abundances[taxa] = []
                        # new taxa not found previously in any organ,
                        # add to dictionary and set all abundances for
                        # other organs to 0.0
                        for j in range(0, num_samples):
                            abundances[taxa].append(0.0)
                    abundances[taxa][sample_id] += abundance
            sample_id += 1

    return (abundances, organs_list, total_organs, metadata)


def output_biom_table(abundances,
                      organs_list,
                      reports_dp,
                      taxa_level,
                      taxonomy_report_dp,
                      metadata):
    """
    """
    output_biom_fp = join(reports_dp, "tcga_output_%s.biom" % (taxa_level))
    with open(output_biom_fp, 'w') as output_f:
        # output BIOM tsv file
        output_f.write("# Constructed from biom file\n")
        output_f.write("# OTU ID\t")
        for file_n in metadata:
            if isfile(join(taxonomy_report_dp, "%s_kraken.report_mpa" % file_n)):
                output_f.write("%s\t" % metadata[file_n][3])
        output_f.write("\n")

        # print out sample IDs
        for taxa in abundances:
            taxa1 = taxa.replace("d__", "k__")
            taxa2 = taxa1.replace("|", ";")
            output_f.write("%s\t" % taxa2)
            for abundance in abundances[taxa]:
                output_f.write("%s\t" % abundance)
            output_f.write("\n")   


def main(argv):

    # kraken-mpa-report filepath
    report_fp = sys.argv[1]

    # metadata TSV file (QIIME format)
    metadata_fp = sys.argv[2]

    # taxonomy report directory path
    taxonomy_report_dp = sys.argv[3]

    # taxonomy level (ex. genus, species)
    taxa_level = sys.argv[4]

    taxa_levels = {"domain": "d__",
                   "phylum": "|p__",
                   "class": "|c__",
                   "order": "|o__",
                   "family": "|f__",
                   "genus": "|g__",
                   "species": "|s__"}

    taxa_levels_idx = {"d__": 0, "|p__": 1, "|c__": 2,
                       "|o__": 3, "|f__": 4, "|g__": 5,
                       "|s__": 6, "6": "|s__", "5": "|g__",
                       "4": "|f__", "3": "|o__", "2": "|c__",
                       "1": "|p__", "0": "d__"}

    abundances, organs_list, total_organs, metadata =\
        compute_abundances_biom(
            reports_dp=reports_dp,
            metadata_fp=metadata_fp,
            taxonomy_report_dp=taxonomy_report_dp,
            taxa_level=taxa_level,
            taxa_levels=taxa_levels,
            taxa_levels_idx=taxa_levels_idx)

    output_biom_table(abundances=abundances,
                      organs_list=organs_list,
                      reports_dp=reports_dp,
                      taxa_level=taxa_level,
                      taxonomy_report_dp=taxonomy_report_dp,
                      metadata=metadata)


if __name__ == "__main__":
    main(sys.argv[1:])
