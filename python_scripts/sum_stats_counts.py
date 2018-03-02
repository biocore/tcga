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


def compute_stats(stats_dp):
    """ BAM file statistics.
    """
    total_reads=0
    total_reads_mapped=0
    total_reads_unmapped=0
    for fn in os.listdir(stats_dp):
        with open(os.path.join(stats_dp, fn)) as f:
            lines = f.read().splitlines()
            total_reads += int(lines[0])
            total_reads_mapped += int(lines[1])
            #total_true=0
            #mapped_true=0
            #for line in f:
            #    if "in total" in line:
            #        line = line.strip().split()[0]
            #        total_reads += int(line)
            #        total_true=1
            #    elif "mapped (" in line:
            #        line = line.strip().split()[0]
            #        total_reads_mapped += int(line)
            #        mapped_true=1
            #if total_true != 1:
            #    print("ERROR: in total doesn't exist %s" % fn)
            #elif mapped_true != 1:
            #    print("ERROR: mapped doesn't exist %s" % fn)
    return total_reads, total_reads_mapped


@click.command()
@click.option('--stats-dp', required=False,
              type=click.Path(resolve_path=True, readable=True, exists=False,
                              file_okay=True),
              help='Directory with stats files per disease')

def main(stats_dp):
    total_reads, total_reads_mapped = compute_stats(stats_dp)
    print("%s\t%s" % (total_reads, total_reads-total_reads_mapped))


if __name__ == "__main__":
    main()