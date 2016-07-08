#!/usr/bin/env python

#-----------------------------------------------------------------------------
# Copyright (c) 2016--, Evguenia Kopylova.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import sys

import skbio
from skbio import Sequence

def main():
    """Output FASTA file using sequences in Kraken output.
    """
    kraken_results_fp = sys.argv[1]
    fasta_fp = sys.argv[2]
    output_fp = sys.argv[3]
    fasta_d = {}
    # Load all FASTA seqs into dictionary
    for seq in skbio.io.read(fasta_fp, format='fasta'):
        label = seq.metadata['id']
        description = seq.metadata['description']
        if label not in fasta_d:
            fasta_d[label] = [description, str(seq)]
        else:
            raise ValueError("Duplicate FASTA labels %s" % label)
    # Output sequences to FASTA.
    with open(output_fp, 'w') as output_f:
        with open(kraken_results_fp) as kraken_results_f:
            for line in kraken_results_f:
                label = line.strip().split('\t')[1]
                if label in fasta_d:
                    output_f.write(">%s%s\n%s\n" % (label, fasta_d[label][0],
                                                    fasta_d[label][1]))


if __name__ == '__main__':
    main()