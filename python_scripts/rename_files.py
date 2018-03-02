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


def rename_files(stats_dp, delimiter):
    for fn in os.listdir(stats_dp):
        if fn.startswith(delimiter):
            os.rename(fn, fn[len(delimiter):])


@click.command()
@click.option('--stats-dp', required=False,
              type=click.Path(resolve_path=True, readable=True, exists=False,
                              file_okay=True),
              help='Directory with stats files per disease')
@click.option('--delimiter', required=False, type=str, default='_1_',
              help='Beginning delimiter of file name')
def main(stats_dp, delimiter):
    rename_files(stats_dp, delimiter)


if __name__ == "__main__":
    main()