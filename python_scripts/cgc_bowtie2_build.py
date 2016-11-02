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
                           taxonomic_rank,
                           taxa_levels,
                           taxa_levels_idx,
                           read_per_taxa):
    """Absolute abundance of number of reads matching a defined taxa level.
    Parameters
    ----------
    kraken_translate_report_fp: str
        filepath to output of "kraken translate"
    taxonomic_rank: str
        taxonomy level (e.g., genus or species)
    taxa_levels: dict
        keys are full name taxonomic ranks and values are abbreviated ranks
    taxa_levels_idx: dict
        2-way dict storing integer depths for abbreviated taxonomic ranks
    read_per_taxa: in
        integer of number of minimum number of reads to keep per taxa

    Returns
    -------
    taxonomies: set
        set of taxonomies from kraken_translate_report_fp
    split_on_level: str
        string that determines the level to split the taxonomy to be used in
        create_db_folder function
    repo_file_name: str
        string with filename suffix from kraken_translate_report_fp to be used
        in create_db_folder function

    """
    repo_file_name = kraken_mpa_report_fp.strip().split('/')[-1].split('mpa')[0]

    taxonomic_rank_level_str = taxa_levels[taxonomic_rank]
    taxonomic_rank_level_int = taxa_levels_idx[taxonomic_rank_level_str]

    if taxonomic_rank_level_int < 6:
        split_on_level = taxa_levels_idx[str(taxonomic_rank_level_int + 1)]
    else:
        # keep full string (to species level)
        split_on_level = '\n'
    taxonomic_abundances= []
    with open(kraken_mpa_report_fp) as kraken_translate_report_f:
        for line in kraken_translate_report_f:
            label, taxonomy = line.strip().split('\t')
            if taxonomic_rank_level_str in taxonomy:
                # keep taxonomy string up to specified level
                taxonomy_parse = taxonomy.split(split_on_level)[0]
                taxonomy_parse = taxonomy_parse.replace('d__','k__')
                taxonomic_abundances.append(taxonomy_parse)

    taxonomies = set([k for k, v in
                      collections.Counter(taxonomic_abundances).iteritems()
                      if v > read_per_taxa])

    return taxonomies, split_on_level, repo_file_name

def create_db_folder(repophlan_genomeid_taxonomy_fp,
                     genome_tax,
                     split_on_level,
                     repo_file_name):

    """Return sets for sample IDs and taxonomy strings.
    Parameters
    ----------
    repophlan_genomeid_taxonomy_fp: str
        filepath to output of repophlan file genome IDs and associated taxa
    genome_tax: set
        set of taxonomies from kraken_translate_report_fp
    split_on_level: str
        string that determines the level to split the taxonomy in genome_tax
    repo_file_name
        string with filename suffix from kraken_translate_report_fp

    Returns
    -------
    Writes a text file with genome IDs, directory to genomes, and associated
    taxa from repophlan_genomeid_taxonomy_fp containned with the set
    genome_tax
    """
    repo_file_name = repo_file_name+'repophlan_genomeID_taxonomy.good'

    with open(repo_file_name, 'w') as f:
        with open(repophlan_genomeid_taxonomy_fp) as repophlan_genomeID_taxonomy_f:
            for line in repophlan_genomeID_taxonomy_f:
                tax_id_repo, directory_repo, taxonomy_repo =\
                    line.strip().split('\t')
                taxonomy_repo = taxonomy_repo.split(split_on_level)[0]
                if taxonomy_repo in genome_tax:
                    f.writelines(line)


@click.command()
@click.option('--kraken-mpa-report-fp', required=True,
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
              help="Minimum number of reads needed for each taxa")


def main(kraken_mpa_report_fp,
         repophlan_genomeid_taxonomy_fp,
         taxonomic_rank,
         read_per_taxa):

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

    taxonomic_set, spit_on_level, repo_file_name =\
        load_kraken_mpa_report(
            kraken_mpa_report_fp=kraken_mpa_report_fp,
            taxonomic_rank=taxonomic_rank,
            taxa_levels=taxa_levels,
            taxa_levels_idx=taxa_levels_idx,
            read_per_taxa=read_per_taxa)

    create_db_folder(
            repophlan_genomeid_taxonomy_fp=repophlan_genomeid_taxonomy_fp,
            genome_tax=taxonomic_set,
            split_on_level=spit_on_level,
            repo_file_name=repo_file_name)

if __name__ == "__main__":
    main()
