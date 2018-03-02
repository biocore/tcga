#!/usr/bin/env python

# ----------------------------------------------------------------------------
# Copyright (c) 2017--, Evguenia Kopylova
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

"""
Read sequence IDs and copy corresponding FASTA record from one corrupt file to another clean file.
"""

import click
import skbio

def process_file(input_corrupt_fasta_fp, input_seqids_fp, output_fp):
    """Write non-corrupted FASTA seqs to file.
    """
    # Store seqids in set
    print("Read seqids ..")
    #seqids = set()
    num_str = 0
    with open(input_seqids_fp) as input_seqids_f:
        for line in input_seqids_f:
            #seqid = line.strip()
            #seqids.add(seqid)
            num_str += 1
            if num_str % 10000000 == 0:
                print("%s seqids processed" % num_str)
    print("Total sequences processed: %s" % num_str)
    #print("Total seqids: %s" % len(seqids))
    # Write
    print("Read through FASTA file ..")
    num_str = 0
    for seq in skbio.io.read(input_corrupt_fasta_fp, format='fasta'):
        seq_id = seq.metadata['id']
        seq_de = seq.metadata['description']
        #if seq_id in seqids:
        with open(output_fp, "a") as output_f:
            output_f.write('>%s %s\n%s\n' % (seq_id, seq_de, seq))
        num_str += 1
        if num_str % 10000000 == 0:
            print("%s reads processed" % num_str)


@click.command()
@click.option('--input-corrupt-fasta-fp', required=True,
              type=click.Path(resolve_path=True, readable=True, exists=True,
                              file_okay=True),
              help='Corrupt FASTA file')
@click.option('--input-seqids-fp', required=True,
              type=click.Path(resolve_path=True, readable=True, exists=True,
                              file_okay=True),
              help='Sequence IDs to keep')
@click.option('--output-fp', required=True,
              type=click.Path(resolve_path=True, readable=True, exists=False,
                              file_okay=True),
              help='Path to write FASTA seqs for seqids')
def main(input_corrupt_fasta_fp, input_seqids_fp, output_fp):

    process_file(input_corrupt_fasta_fp=input_corrupt_fasta_fp,
                 input_seqids_fp=input_seqids_fp,
                 output_fp=output_fp)


if __name__ == "__main__":
    main()