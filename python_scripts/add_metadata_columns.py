#!/usr/bin/env python

#-----------------------------------------------------------------------------
# Copyright (c) 2016--, Evguenia Kopylova.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

"""
Add metadata columns from super metadata file to reduced metadata file.
"""

import click
from collections import OrderedDict


def add_columns(qiime_mapping_file_fp,
                cgc_mapping_file_fp,
                metadata_column,
                output_fp):
    """Copy metadata columns from CGC metadata file to QIIME metadata file.

    Parameters
    ----------
    qiime_mapping_file_fp: str
        Filepath to QIIME mapping file
    cgc_mapping_file_fp: str
        Filepath to CGC mapping file
    metadata_column: tuple
        Metadata columns to copy from CGC to QIIME mapping file
    output_fp: str
        Filepath to updated QIIME mapping file
    """
    qiime_mapping_file = {}
    qiime_mapping_header = ""
    with open(qiime_mapping_file_fp) as qiime_f:
        for line in qiime_f:
            line = line.strip().split('\t')
            if line[0] == '#SampleID':
                qiime_mapping_header = line
                continue
            # use filename as key
            key = line[3]
            qiime_mapping_file[key] = line

    cgc_mapping_header = []
    cgc_metadata_columns = OrderedDict()
    # initialize list indexes to 0
    for name in metadata_column:
        cgc_metadata_columns[str(name)] = 0
    files_searched = set()
    with open(cgc_mapping_file_fp) as cgc_f:
        for line in cgc_f:
            line = line.strip().split('\t')
            if line[0] == "case_name":
                cgc_mapping_header = line
                # set correct index values for metadata columns in CGC header
                for key in cgc_metadata_columns:
                    cgc_metadata_columns[key] = cgc_mapping_header.index(key)+1
                continue
            # use filename as key
            key_file = line[22]
            if key_file in files_searched:
                continue
            if key_file in qiime_mapping_file:
                for key_column, value in cgc_metadata_columns.iteritems():
                    qiime_mapping_file[key_file].append(line[value])
                files_searched.add(key_file)

    # output updated QIIME mapping file
    num_columns_qiime = len(qiime_mapping_header)
    num_columns_cgc = len(cgc_metadata_columns)
    with open(output_fp, 'w') as output_f:
        output_f.write('\t'.join(map(str,qiime_mapping_header[:-1])))
        for key in cgc_metadata_columns:
            output_f.write('\t%s' % key)
        output_f.write('\t%s' % qiime_mapping_header[-1])
        output_f.write('\n')
        for key in qiime_mapping_file:
            # Output QIIME mapping file up to Description column
            output_f.write('\t'.join(map(str,qiime_mapping_file[key][:num_columns_qiime-1])))
            output_f.write('\t')
            # Output added keys from CGC mapping file
            output_f.write('\t'.join(map(str,qiime_mapping_file[key][num_columns_qiime:])))
            # Output Description column
            output_f.write('\t%s' % qiime_mapping_file[key][num_columns_qiime-1])
            output_f.write('\n')


@click.command()
@click.option('--qiime-mapping-file-fp', required=False,
              type=click.Path(resolve_path=True, readable=True, exists=True,
                              file_okay=True),
              help='Filepath to QIIME metadata file')
@click.option('--cgc-mapping-file-fp', required=True,
              type=click.Path(resolve_path=True, readable=True, exists=True,
                              file_okay=True),
              help='Filepath to CGC metadata file')
@click.option('--metadata-column', type=str, required=True, multiple=True,
              help="Metadata column to add from CGC to QIIME mapping file")
@click.option('--output-fp', required=True,
              type=click.Path(resolve_path=True, readable=True, exists=False,
                              file_okay=True),
              help='Filepath to updated QIIME mapping file')
def main(qiime_mapping_file_fp,
         cgc_mapping_file_fp,
         metadata_column,
         output_fp):
    # rename sample ids in mapping file
    add_columns(qiime_mapping_file_fp,
                cgc_mapping_file_fp,
                metadata_column,
                output_fp)


if __name__ == "__main__":
    main()