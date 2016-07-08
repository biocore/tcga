#!/bin/python

#-----------------------------------------------------------------------------
# Copyright (c) 2016--, Evguenia Kopylova.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

# Edit the FASTA label to Kraken format:
#   >sequence16|kraken:taxid|32630  Adapter sequence
#   CAAGCAGAAGACGGCATACGAGATCTTCGAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA

import sys
import bz2
from skbio import Sequence
import skbio

from os import walk
from os.path import join, splitext, basename, isfile


def get_genome_paths(scores_average_fp, scores_repophlan_fp):
    """Return genome ID, tax_id and .fna.bz2 path
    """
    genomes = {}
    # Get quality filtered genome ID and filepath to .fna.bz2
    with open(scores_average_fp) as scores_average_f:
        # header
        line = scores_average_f.readline()
        line = line.strip().split('\t')
        fna_idx = line.index('fna_lname')
        for line in scores_average_f:
            line = line.strip().split('\t')
            genome_id = line[0]
            fna_fp = line[fna_idx]
            if genome_id not in genomes:
                genomes[genome_id] = [fna_fp]
            else:
                raise ValueError("Duplicate genome IDs %s" % genome_id)
    # Get tax_id
    genomes_without_tax_id = []
    with open(scores_repophlan_fp) as scores_repophlan_f:
        # header
        line = scores_repophlan_f.readline()
        line = line.strip().split('\t')
        tax_id_idx = line.index('taxid')
        for line in scores_repophlan_f:
            line = line.strip().split('\t')
            genome_id = line[0]
            # only want tax_ids for genomes passing quality filter
            if genome_id in genomes:
                tax_id = line[tax_id_idx]
                # tax_id must be an integer, if not check the field
                # immediately prior (bug in repophlan parsing)
                if not tax_id.isdigit():
                    tax_id = line[tax_id_idx-1]
                    if not tax_id.isdigit():
                        genomes_without_tax_id.append(genome_id)
                        continue
                genomes[genome_id].append(int(tax_id))
    return genomes, genomes_without_tax_id


def format_labels(genomes, repophlan_scores_filtered_genomes_dp):
    """Format FASTA labels in qualified genomes.
    """
    kraken_db_str = "|kraken:taxid|"
    for genome_id, info in genomes.items():
        genome_fp = info[0]
        taxid = info[1]
        genome_fp_name = basename(splitext(genome_fp)[0])
        # check FNA file exists for genome
        if genome_fp_name != "":
            output_fp = join(repophlan_scores_filtered_genomes_dp,
                             genome_fp_name)
            # skip files already modified (e.g. from previous run)
            if not isfile(output_fp):
                with open(output_fp, 'w') as output_f:
                    for seq in skbio.io.read(
                            genome_fp, format='fasta', compression='bz2'):
                        id_ = seq.metadata['id']
                        description = seq.metadata['description']
                        new_id_ = "%s%s%s" % (id_, kraken_db_str, taxid)
                        seq.metadata['id'] = new_id_
                        output_f.write(
                            ">%s%s%s %s\n%s\n" % (id_, kraken_db_str,
                                                  taxid, description,
                                                  str(seq)))


def main():
    """Edit qualified genomes' labels to Kraken format.
    """
    # .fna.bz2 genomes folder
    all_genomes_bz2_dp = sys.argv[1]
    # Output of RepoPhlAn's repophlan_get_microbes.py
    repophlan_scores_fp = sys.argv[2]
    # Output of Jon's score average script
    # https://github.com/tanaes/script_bin/blob/master/filter_repophlan.py
    repophlan_scores_average_fp = sys.argv[3]
    # Output directory path
    repophlan_scores_filtered_genomes_dp = sys.argv[4]

    genomes, genomes_without_tax_id = get_genome_paths(
        repophlan_scores_average_fp, repophlan_scores_fp)
    if len(genomes_without_tax_id) != 0:
        raise ValueError(
            "Some genomes do not have a taxid: %s" % genomes_without_tax_id)

    format_labels(genomes, repophlan_scores_filtered_genomes_dp)


if __name__ == '__main__':
    main()