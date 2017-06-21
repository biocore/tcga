#!/usr/bin/env python

# ----------------------------------------------------------------------------
# Copyright (c) 2017--, Evguenia Kopylova
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

"""
Split FASTA file into two files with WGS and RNA-Seq reads.
"""

import click
from os.path import join, basename
import skbio

def read_metadata(metadata_file_fp):
    """Read metadata into dictionary.
    """
    metadata = {}
    with open(metadata_file_fp) as metadata_file_f:
        for line in metadata_file_f:
            line = line.strip().split('\t')
            sample_id = line[0]
            if sample_id not in metadata:
                metadata[sample_id] = line[1:]
            else:
                raise ValueError("Duplicate samples: %s" % sample_id)
    return metadata


def split_file(input_fasta_fp,
               output_dir_dp,
               metadata):
    """Split FASTA file (QIIME formatted labels) into WGS and RNA-Seq files.
    """
    print(input_fasta_fp)
    output_file_base = basename(input_fasta_fp)
    output_wgs_fp = join(output_dir_dp, "%s_wgs.fasta" % output_file_base)
    output_rnaseq_fp = join(output_dir_dp, "%s_rnaseq.fasta" % output_file_base)
    with open(output_wgs_fp, 'w') as output_wgs_f:
        with open(output_rnaseq_fp, 'w') as output_rnaseq_f:       
            for seq in skbio.io.read(input_fasta_fp, format='fasta'):
                seq_id = seq.metadata['id']
                seq_desc = seq.metadata['desc']
                sample_id = seq_id.split()[0].split('_')[0]
                experimental_strategy = metadata[sample_id][15]
                if experimental_strategy  == "WGS":
                    output_wgs_f.write(">%s %s\n%s\n" % (seq_id, seq_desc, str(seq)))
                elif experimental_strategy == "RNA-Seq":
                    output_rnaseq_f.write(">%s %s\n%s\n" % (seq_id, seq_desc, str(seq)))
                else:
                    raise ValueError("Unknown experimental_strategy: %s" % experimental_strategy )


@click.command()
@click.option('--input-fasta-fp', required=True,
              type=click.Path(resolve_path=True, readable=True, exists=True,
                              file_okay=True),
              help='FASTA file containing WGS and RNA-Seq reads')
@click.option('--metadata-file-fp', required=True,
              type=click.Path(resolve_path=True, readable=True, exists=False,
                              file_okay=True),
              help='Filepath to QIIME formatted metadata')
@click.option('--output-dir-dp', required=True,
              type=click.Path(resolve_path=True, readable=True, exists=False,
                              file_okay=True),
              help='Output directory path for WGS and RNA-Seq files')
def main(input_fasta_fp,
         metadata_file_fp,
         output_dir_dp):

    metadata = read_metadata(metadata_file_fp)
    split_file(input_fasta_fp=input_fasta_fp,
               output_dir_dp=output_dir_dp,
               metadata=metadata)


if __name__ == "__main__":
    main()
