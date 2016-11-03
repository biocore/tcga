#!/usr/bin/env python

#-----------------------------------------------------------------------------
# Copyright (c) 2016--, Evguenia Kopylova, Jad Kanbar
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

"""
Build Bowtie2 database on all reference genomes in Kraken report.
"""

import click
import collections

def load_kraken_mpa_report(kraken_mpa_report_fp,
                           taxa_levels,
                           read_per_taxa):
    """Absolute abundance of number of reads matching a defined taxa level.
    Parameters
    ----------
    kraken_translate_report_fp: str
        filepath to output of "kraken translate"
    taxa_levels: list
        list of two elements that includes the taxonomic rank at which
        to generate summary and rank below to split by
    read_per_taxa: int
        integer of number of minimum number of reads to keep per taxa

    Returns
    -------
    taxonomies: set
        set of taxonomies from kraken_translate_report_fp
    """
    taxonomic_abundances= []
    for report_fp in kraken_mpa_report_fp:
        with open(report_fp) as report_fp:
            for line in report_fp:
                label, taxonomy = line.strip().split('\t')
                if taxa_levels[0] in taxonomy:
                    # keep taxonomy string up to specified level
                    taxonomy_parse = taxonomy.split(taxa_levels[1])[0]
                    taxonomy_parse = taxonomy_parse.replace('d__','k__')
                    taxonomic_abundances.append(taxonomy_parse)

    taxonomies = set([k for k, v in
                      collections.Counter(taxonomic_abundances).iteritems()
                      if v > read_per_taxa])

    return taxonomies

def create_db_folder(repophlan_genomeid_taxonomy_fp,
                     genome_tax,
                     split_on_level,
                     output_filename):

    """Return sets for sample IDs and taxonomy strings.
    Parameters
    ----------
    repophlan_genomeid_taxonomy_fp: str
        filepath to output of repophlan file genome IDs and associated taxa
    genome_tax: set
        set of taxonomies from kraken_translate_report_fp
    split_on_level: str
        string that determines the level to split the taxonomy in genome_tax
    output_filename: str
        string with output filename

    Returns
    -------
    Writes a text file with genome IDs, directory to genomes, and associated
    taxa from repophlan_genomeid_taxonomy_fp containned with the set
    genome_tax
    """
    with open(output_filename, 'w') as f:
        with open(repophlan_genomeid_taxonomy_fp) as repophlan_genomeID_taxonomy_f:
            for line in repophlan_genomeID_taxonomy_f:
                tax_id_repo, directory_repo, taxonomy_repo =\
                    line.strip().split('\t')
                taxonomy_repo = taxonomy_repo.split(split_on_level)[0]
                if taxonomy_repo in genome_tax:
                    f.writelines(line)


@click.command()

@click.option('--kraken-mpa-report-fp', required=True, multiple=True,
              type=click.Path(resolve_path=True, readable=True, exists=True,
                              file_okay=True),
              help='Filepath to Kraken report')
@click.option('--repophlan-genomeID-taxonomy-fp', required=True,
              type=click.Path(resolve_path=True, readable=True, exists=False,
                              file_okay=True),
              help='Filepath to RepoPhlAn genome ID and taxonomy list')
@click.option('--taxonomic-rank', type=click.Choice(['genus', 'species',
                                                     'family', 'order',
                                                     'class', 'phylum',
                                                     'domain']),
              required=False, default=['genus'], show_default=True,
              help="Taxonomic rank at which to generate summary")
@click.option('--read-per-taxa', required=False,
              default=10, show_default=True,
              help="Minimum number of reads for each taxa")
@click.option('--output-filename', required=False,
              default='subset_repophlan_genomeID_taxonomy.good',
              show_default=True,
              help="Output filename for RepoPhlAn genome ID and taxonomy list")

def main(kraken_mpa_report_fp,
         repophlan_genomeid_taxonomy_fp,
         taxonomic_rank,
         read_per_taxa,
         output_filename):

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

    taxa_levels_str=taxa_levels[taxonomic_rank]
    taxa_levels_idx_int=taxa_levels_idx[taxa_levels_str]

    if taxa_levels_idx_int < 6:
        split_on_level = taxa_levels_idx[str(taxa_levels_idx_int + 1)]
    else:
        split_on_level = '\t'

    taxonomic_set =\
        load_kraken_mpa_report(
            kraken_mpa_report_fp=kraken_mpa_report_fp,
            taxa_levels=[taxa_levels_str,split_on_level],
            read_per_taxa=read_per_taxa)

    create_db_folder(
            repophlan_genomeid_taxonomy_fp=repophlan_genomeid_taxonomy_fp,
            genome_tax=taxonomic_set,
            split_on_level=split_on_level,
            output_filename=output_filename)

if __name__ == "__main__":
    main()
