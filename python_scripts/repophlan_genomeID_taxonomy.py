#!/bin/python

#-----------------------------------------------------------------------------
# Copyright (c) 2016--, Evguenia Kopylova.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

# Generate tsv file with RepoPhlAn quality genomes and their taxonomies.

import sys
import bz2
from skbio import Sequence
import skbio
import click

from os import walk
from os.path import join, splitext, basename, isfile


def get_genome_paths(scores_average_fp, scores_repophlan_fp):
    """Return genome ID, tax_id and .fna.bz2 path.

    Parameters
    ----------
    scores_average_fp: str
        RepoPhlAn screened report
    scores_repophlan_fp: str
        RepoPhlAn's output report from repophlan_get_microbes.py

    Returns
    -------
    genomes: dict
        keys are RepoPhlAn genome IDs and values are list of taxonomy str and
        FASTA filepath
    genomes_without_tax_id: list
        RepoPhlAn genome IDs without a listed taxonomy
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
    # Get taxonomy
    genomes_without_taxonomy = []
    with open(scores_repophlan_fp) as scores_repophlan_f:
        # header
        line = scores_repophlan_f.readline()
        line = line.strip().split('\t')
        tax_id_idx = line.index('taxonomy')
        for line in scores_repophlan_f:
            line = line.strip().split('\t')
            genome_id = line[0]
            # only want tax_ids for genomes passing quality filter
            if genome_id in genomes:
                taxonomy = line[tax_id_idx]
                # tax_id must be an integer, if not check the field
                # immediately prior (bug in repophlan parsing)
                if 'k__' not in taxonomy:
                    taxonomy = line[tax_id_idx-1]
                    if 'k__' not in taxonomy:
                        genomes_without_taxonomy.append(genome_id)
                        continue
                genomes[genome_id].append(taxonomy)
    return genomes, genomes_without_taxonomy


def output_results(genomes, output_fp):
    """Format FASTA labels in qualified genomes.

    Parameters
    ----------
    genomes: dict
        keys are RepoPhlAn genome IDs and values are list of taxonomy str and
        FASTA filepath
    output_fp: str
        Output filepath        
    """
    with open(output_fp, 'w') as output_f:
        for genome, info in genomes.items():
            output_f.write('%s\t%s\t%s\n' % (genome, info[0], info[1]))


@click.command()
@click.option('--repophlan-scores-fp', required=True,
              type=click.Path(resolve_path=True, readable=True, exists=True,
                              file_okay=True),
              help="RepoPhlAn's output report from repophlan_get_microbes.py")
# Output of Jon's score average script
# https://github.com/tanaes/script_bin/blob/master/filter_repophlan.py
@click.option('--repophlan-scores-average-fp', required=True,
              type=click.Path(resolve_path=True, readable=True, exists=True,
                              file_okay=True),
              help="RepoPhlAn's screened report")
@click.option('--output-fp', required=True,
              type=click.Path(resolve_path=True, readable=True, exists=False,
                              file_okay=True),
              help="Output filepath")
def main(repophlan_scores_fp,
         repophlan_scores_average_fp,
         output_fp):
    """Edit qualified genomes' labels to Kraken format.
    """
    genomes, genomes_without_taxonomy = get_genome_paths(
        repophlan_scores_average_fp, repophlan_scores_fp)
    if len(genomes_without_taxonomy) != 0:
        raise ValueError(
            "Some genomes do not have a taxid: %s" % genomes_without_taxonomy)

    output_results(genomes, output_fp)


if __name__ == '__main__':
    main()