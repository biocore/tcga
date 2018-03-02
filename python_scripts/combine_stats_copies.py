#!/usr/bin/env python

# ----------------------------------------------------------------------------
# Copyright (c) 2017--, Evguenia Kopylova
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import click
import os


def compute_stats(stats_dp, output_stats_unique_numbers, output_stats_unique_charact):
    for dp in stats_dp:
        for fn in os.listdir(dp):
            with open(os.path.join(dp, fn)) as f:
                



@click.command()
@click.option('--stats-dp', required=False, multiple=True,
              type=click.Path(resolve_path=True, readable=True, exists=False,
                              file_okay=True),
              help='Directory with stats files per disease')
@click.option('--output-stats-unique-numbers', required=True,
              type=click.Path(resolve_path=True, readable=True, exists=False,
                              file_okay=True),
              help='Directory with unique stats files per disease (numbers)')
@click.option('--output-stats-unique-charact', required=True,
              type=click.Path(resolve_path=True, readable=True, exists=False,
                              file_okay=True),
              help='Directory with unique stats files per disease (characters)')                     

def main(stats_dp, output_stats_unique_numbers, output_stats_unique_charact):
    combine_stats(stats_dp, output_stats_unique_numbers, output_stats_unique_charact)


if __name__ == "__main__":
    main()