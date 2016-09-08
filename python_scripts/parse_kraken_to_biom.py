#!/usr/bin/env python

#-----------------------------------------------------------------------------
# Copyright (c) 2016--, Evguenia Kopylova.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

"""
Parse output of "kraken translate" and generate single BIOM table.
"""

import click
import pandas as pd
import numpy as np
from biom.table import Table
from biom.util import biom_open


def compute_biom_table(kraken_translate_report_fp,
                       taxonomic_rank,
                       taxa_levels,
                       taxa_levels_idx,
                       columns,
                       index):
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
    columns: list
        sample IDs list
    index: list
        observation IDs (taxonomy strings) list

    Returns
    -------
    biom_table: biom.Table
        BIOM table
    """
    total_levels = len(taxa_levels)
    # columns are sample IDs and rows are taxonomy strings
    abundances = pd.DataFrame(columns=columns, index=index)    
    taxonomic_rank_level_str = taxa_levels[taxonomic_rank]
    taxonomic_rank_level_int = taxa_levels_idx[taxonomic_rank_level_str]
    if taxonomic_rank_level_int < 6:
        split_on_level = taxa_levels_idx[str(taxonomic_rank_level_int + 1)]
    else:
        # keep full string (to species level)
        split_on_level = '\t'
    with open(kraken_translate_report_fp) as kraken_translate_report_f:
        for line in kraken_translate_report_f:
            label, taxonomy = line.strip().split('\t')
            # record abundance
            if taxonomic_rank_level_str in taxonomy:
                # keep taxonomy string up to specified level
                taxonomy = taxonomy.split(split_on_level)[0]
                sample_id = label.split('_')[0]
                value = abundances.at[taxonomy, sample_id]
                if np.isnan(value):
                    abundances.set_value(taxonomy, sample_id, 1.)
                else:
                    abundances.set_value(taxonomy, sample_id, value+1.)
    obs_ids = abundances.index.values.tolist()
    for i in range(len(obs_ids)):
        obs_ids[i] = obs_ids[i].replace("d__", "k__")
        obs_ids[i] = obs_ids[i].replace("|", ";")
    return Table(abundances.fillna(0.).as_matrix(),
                 obs_ids,
                 abundances.columns.values.tolist())


def prepare_dataframe(kraken_translate_report_fp,
                      taxonomic_rank,
                      taxa_levels,
                      taxa_levels_idx):
    """Return sets for sample IDs and taxonomy strings.

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

    Returns
    -------
    sample_ids: list
        all unique sample IDs in file
    taxonomies: list
        all unique taxonomies in file
    """
    total_levels = len(taxa_levels)
    taxonomic_rank_level_str = taxa_levels[taxonomic_rank]
    taxonomic_rank_level_int = taxa_levels_idx[taxonomic_rank_level_str]
    if taxonomic_rank_level_int < 6:
        split_on_level = taxa_levels_idx[str(taxonomic_rank_level_int + 1)]
    else:
        # keep full string (to species level)
        split_on_level = '\t'
    sample_ids = set()
    taxonomies = set()
    with open(kraken_translate_report_fp) as kraken_translate_report_f:
        for line in kraken_translate_report_f:
            label, taxonomy = line.strip().split('\t')
            sample_id = label.split('_')[0]
            sample_ids.add(sample_id)
            # record abundance
            if taxonomic_rank_level_str in taxonomy:
                # keep taxonomy string up to specified level
                taxonomy = taxonomy.split(split_on_level)[0]
                taxonomies.add(taxonomy)
    return list(sample_ids), list(taxonomies)


def write_biom_table(table, biom_output_fp):
    """Write BIOM table to file.

    Parameters
    ----------
    table: biom.Table
        an instance of a BIOM table
    biom_output_fp: str
        filepath to output BIOM table
    """
    with biom_open(biom_output_fp, 'w') as f:
            table.to_hdf5(h5grp=f, generated_by="tcga-kraken-translate")


@click.command()
@click.option('--kraken-translate-report-fp', required=True,
              type=click.Path(resolve_path=True, readable=True, exists=True,
                              file_okay=True))
@click.option('--taxonomic-rank', type=click.Choice(['genus', 'species',
                                                     'family', 'order',
                                                     'class', 'phylum',
                                                     'domain']),
              required=False, default=['genus'], show_default=True,
              help="Taxonomic rank at which to generate summary")
@click.option('--biom-output-fp', required=True,
              type=click.Path(resolve_path=True, readable=True, exists=False,
                              file_okay=True),
              help="Filepath to output BIOM table")
def main(kraken_translate_report_fp,
         taxonomic_rank,
         biom_output_fp):
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
    columns, index = prepare_dataframe(
        kraken_translate_report_fp=kraken_translate_report_fp,
        taxonomic_rank=taxonomic_rank,
        taxa_levels=taxa_levels,
        taxa_levels_idx=taxa_levels_idx)
    biom_table =\
        compute_biom_table(
            kraken_translate_report_fp=kraken_translate_report_fp,
            taxonomic_rank=taxonomic_rank,
            taxa_levels=taxa_levels,
            taxa_levels_idx=taxa_levels_idx,
            columns=columns,
            index=index)
    # output table
    write_biom_table(biom_table, biom_output_fp)


if __name__ == "__main__":
    main()
