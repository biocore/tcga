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


def validate_stats(stats_dp, data_dp):
    """ Validate same files used for both stats.
    """
    data_bam_files = []
    stats_txt_files = []
    for fn in os.listdir(stats_dp):
        fn_b = fn.replace("_stats.txt", "")
        stats_txt_files.append(fn_b)
    for fn in os.listdir(data_dp):
        if "mapping_file" in fn:
            with open(os.path.join(data_dp,fn)) as f:
                next(f)
                for line in f:
                    line = line.strip().split('\t')
                    line2 = [x for x in line if x] 
                    file = line2[1].replace(".bam", "")
                    data_bam_files.append(file)
    print("Total stats files: %s" % len(stats_txt_files))
    print("Total data files: %s" % len(data_bam_files))
    for item in stats_txt_files:
        if item not in data_bam_files:
            print("%s in stats but not in data" % item)
    for item in data_bam_files:
       if item not in stats_txt_files:
            print("%s in data but not in stats" % item)


@click.command()
@click.option('--stats-dp', required=False,
              type=click.Path(resolve_path=True, readable=True, exists=False,
                              file_okay=True),
              help='Directory with stats files per disease')
@click.option('--data-dp', required=False,
              type=click.Path(resolve_path=True, readable=True, exists=False,
                              file_okay=True),
              help='Directory with data files per disease')

def main(stats_dp,
         data_dp):
    validate_stats(stats_dp, data_dp)


if __name__ == "__main__":
    main()