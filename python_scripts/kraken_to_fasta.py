#!/usr/bin/env python

# -----------------------------------------------------------------------------
# Copyright (c) 2017--, Evguenia Kopylova
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# -----------------------------------------------------------------------------

"""
Output FASTA file for classified and/or unclassified reads in Kraken output.
"""

import click
import skbio.io


def get_kraken_report(kraken_output_fp,
                      classification_id):
    """Parse Kraken report from output file.

    Parameters
    ----------
    kraken_output: str
        Filepath to Kraken output file
    classification_id: str
        Classification ID for which to output reads
    Returns
    -------
    kraken_report: set
        read IDs to output to FASTA
    """
    kraken_report = set()
    with open(kraken_output_fp) as f:
        for line in f:
            line = line.strip().split('\t')
            classification = line[0]
            if (classification != 'U') and (classification != 'C'):
                raise ValueError(
                    'Classification must be a U or C: %s' % classification)
            if classification == classification_id:
                read_id = line[1]
                if read_id in kraken_report:
                    raise ValueError('Duplicate read IDs: %s' % read_id)
                kraken_report.add(read_id)
    return kraken_report


def output_reads(kraken_report, input_fp, output_fp):
    """Output FASTA file for read IDs in Kraken report.

    Parameters
    ----------
    kraken_report: dict
        keys are read IDs and values are classification info
    input_fp: str
        Filepath to FASTA file from which to extract reads
    """
    with open(output_fp, 'w') as output_f:
        with open(input_fp) as input_f:
            for seq in skbio.io.read(input_f, format='fasta'):
                if seq.metadata['id'] in kraken_report:
                    seq.write(output_f)


@click.command()
@click.option('--kraken-output-fp', required=True,
              type=click.Path(resolve_path=True, readable=True, exists=False,
                              file_okay=True),
              help='Filepath to Kraken output')
@click.option('--input_fp', required=True,
              type=click.Path(resolve_path=True, readable=True, exists=False,
                              file_okay=True),
              help='Filepath to FASTA file from which to extract reads')
@click.option('--output-classified', required=True, type=bool, default=True,
              show_default=True, help='Output classified reads to FASTA')
@click.option('--output-unclassified', required=False, type=bool,
              default=False, show_default=False,
              help='Output unclassified reads to FASTA')
def main(kraken_output_fp,
         input_fp,
         output_classified,
         output_unclassified):
    classification_id = 'C'
    if output_classified and output_unclassified:
        raise ValueError(
            'Only one allowed: --output-classified or --output-unclassified')
    if output_unclassified:
        classification_id = 'U'
    kraken_report = get_kraken_report(kraken_output_fp, classification_id)
    output_reads(kraken_report, input_fp)


if __name__ == "__main__":
    main()
