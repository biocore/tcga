#!/usr/bin/env python

# ----------------------------------------------------------------------------
# Copyright (c) 2017--, Evguenia Kopylova
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

"""
Read and copy lines from a corrupt file to a clean file.
"""

import click
import skbio

def process_file(input_corrupt_kraken_output_fp, output_fp):
    """Write non-corrupted Kraken output lines to file.
    """
    print("Read through Kraken output file ..")
    num_str = 0
    with open(input_corrupt_kraken_output_fp) as input_f:
        for line in input_f:
            with open(output_fp, "a") as output_f:
                output_f.write(line)
            num_str += 1
            if num_str % 10000000 == 0:
                print("%s lines processed" % num_str)


@click.command()
@click.option('--input-corrupt-kraken-output-fp', required=True,
              type=click.Path(resolve_path=True, readable=True, exists=True,
                              file_okay=True),
              help='Corrupt Kraken output file')
@click.option('--output-fp', required=True,
              type=click.Path(resolve_path=True, readable=True, exists=False,
                              file_okay=True),
              help='Path to write FASTA seqs for seqids')
def main(input_corrupt_kraken_output_fp, output_fp):

    process_file(input_corrupt_kraken_output_fp=input_corrupt_kraken_output_fp,
                 output_fp=output_fp)


if __name__ == "__main__":
    main()