#!/bin/python

# ----------------------------------------------------------------------------
# Copyright (c) 2016--, Evguenia Kopylova.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

# Move low quality genomes from input folder to output folder.

import click
from os.path import join, basename, exists
from os import listdir, makedirs
import shutil
import logging


def get_logger(log_fp):
    """ Get logging instance.
    """
    logger = logging.getLogger('log')
    log_handler = logging.FileHandler(log_fp)
    formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
    log_handler.setFormatter(formatter)
    logger.addHandler(log_handler)
    logger.setLevel(logging.DEBUG)

    return logger


def get_genomes_ids(repophlan_scores_fp):
    """Return genome IDs to keep.
    """
    genomes = set()
    # Get quality filtered genome ID and filepath to .fna.bz2
    with open(repophlan_scores_fp) as repophlan_scores_f:
        # header
        line = repophlan_scores_f.readline()
        line = line.strip().split('\t')
        fna_idx = line.index('faa_lname')
        for line in repophlan_scores_f:
            line = line.strip().split('\t')
            filename = basename(line[fna_idx])
            if filename not in genomes:
                genomes.add(filename)
            else:
                raise ValueError('Duplicate file name: %s' % filename)
    return genomes


def filter_genomes(all_genomes_bz2_dp,
                   genomes,
                   low_quality_genomes_dp,
                   logger):
    """Move low quality genomes into another directory.
    """
    files_moved = 0
    files_kept = 0
    for filename in listdir(all_genomes_bz2_dp):
        if filename.endswith('.faa.bz2'):
            if filename not in genomes:
                src = join(all_genomes_bz2_dp, filename)
                dst = join(low_quality_genomes_dp, filename)
                logger.info('Moving %s' % filename)
                shutil.move(src, dst)
                files_moved += 1
            else:
                files_kept += 1
        else:
            continue
    logger.info('Total files moved to %s: %s' % (
        low_quality_genomes_dp, files_moved))
    logger.info('Total files kept in %s: %s' % (
        all_genomes_bz2_dp, files_kept))


@click.command()
@click.option('--all-genomes-bz2-dp', required=True,
              type=click.Path(resolve_path=True, readable=True, exists=True,
                              file_okay=True),
              help='Directory path to all 54K genomes (.faa)')
@click.option('--repophlan-scores-fp', required=True,
              type=click.Path(resolve_path=True, readable=True, exists=False,
                              file_okay=True),
              help='Filepath to RepoPhlAn screened results')
@click.option('--low-quality-genomes-dp', required=True,
              type=click.Path(resolve_path=True, readable=True, exists=False,
                              file_okay=True),
              help='Output directory path to low quality genomes')
@click.option('--log-fp', required=True,
              type=click.Path(resolve_path=True, readable=True, exists=False,
                              file_okay=True),
              help='Filepath to log file')
def main(all_genomes_bz2_dp,
         repophlan_scores_fp,
         low_quality_genomes_dp,
         log_fp):
    logger = get_logger(log_fp)
    genomes = get_genomes_ids(repophlan_scores_fp)
    if not exists(low_quality_genomes_dp):
        makedirs(low_quality_genomes_dp)
    filter_genomes(all_genomes_bz2_dp,
                   genomes,
                   low_quality_genomes_dp,
                   logger)

if __name__ == '__main__':
    main()
