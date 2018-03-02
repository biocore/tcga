#!/usr/bin/env python

#-----------------------------------------------------------------------------
# Copyright (c) 2016--, Evguenia Kopylova.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

"""
Rename sample IDs in mapping file and BIOM table for combined analysis of
multiple disease types in CGC TCGA analysis.
"""

import click
from os.path import splitext, join
from biom.table import Table
from biom.util import biom_open
from biom import load_table


def rename_sample_ids(mapping_file_fp,
                      biom_fp,
                      count_start):
    """Rename sample IDs to join mapping files

    Parameters
    ----------
    mapping_file_fp: tuple
        Filepath to mapping file
    biom_fp: str
        Filepath to BIOM table
    count_start: int
        First new sample ID name (ascending in order)
    """
    output_mapping_file_fp = "%s_modified.txt" % splitext(mapping_file_fp)[0]
    id_map = {}
    modified_id = count_start
    with open(output_mapping_file_fp, 'w') as output_f:
        with open(mapping_file_fp, 'r') as mapping_f:
            for line in mapping_f:
                if line.startswith('#SampleID'):
                    output_f.write(line)
                    continue
                line = line.strip().split('\t')
                curr_sample_id = line[0]
                new_sample_id = "s%s" % modified_id
                id_map[curr_sample_id] = new_sample_id
                line[0] = new_sample_id
                output_f.write('\t'.join(map(str,line)))
                output_f.write('\n')
                modified_id += 1

    # update IDs in BIOM table to match modified mapping file
    output_biom_fp = "%s_modified.biom" % splitext(biom_fp)[0]
    table = load_table(biom_fp)
    table.update_ids(id_map, axis='sample')
    with biom_open(output_biom_fp, 'w') as f:
            table.to_hdf5(h5grp=f, generated_by="tcga-rename-sample-ids")


@click.command()
@click.option('--mapping-file-fp', required=False,
              type=click.Path(resolve_path=True, readable=True, exists=True,
                              file_okay=True),
              help='Filepath to mapping file')
@click.option('--biom-fp', required=True,
              type=click.Path(resolve_path=True, readable=True, exists=True,
                              file_okay=True),
              help='Filepath to BIOM table')
@click.option('--count-start', required=True, type=int,
              help='First new sample ID name (ascending in order)')
def main(mapping_file_fp,
         biom_fp,
         count_start):
    # rename sample ids in mapping file
    rename_sample_ids(mapping_file_fp,
                      biom_fp,
                      count_start)


if __name__ == "__main__":
    main()
