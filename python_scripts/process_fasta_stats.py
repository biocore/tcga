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

def process_file(input_fasta_fp):
    """Get file stats.
    """
    min_len = 10000000000000
    max_len = 0
    tot_len = 0
    num_str = 0
    num_M = 0
    for seq in skbio.io.read(input_fasta_fp, format='fasta'):
        num_str += 1
        if num_str % 1000000 == 0:
        	num_M += 1
        	print("%sM reads processed" % num_M)
        seq_id = seq.metadata['id']
        str_len = len(str(seq)) 
        if str_len < min_len:
        	min_len = str_len
        if str_len > max_len:
        	max_len = str_len
        tot_len = tot_len + str_len
    avg_len = tot_len/num_str
    return min_len, max_len, num_str, avg_len


@click.command()
@click.option('--input-fasta-fp', required=True,
              type=click.Path(resolve_path=True, readable=True, exists=True,
                              file_okay=True),
              help='FASTA file containing WGS and RNA-Seq reads')
def main(input_fasta_fp):

    min_len, max_len, num_str, avg_len = process_file(input_fasta_fp=input_fasta_fp)
    print("Total reads: %s" % num_str)
    print("Minimum length read: %s" % min_len)
    print("Maximum length read: %s" % max_len)
    print("Average length read: %s" % avg_len)


if __name__ == "__main__":
    main()