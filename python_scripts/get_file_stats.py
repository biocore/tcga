#!/usr/bin/env python

# ----------------------------------------------------------------------------
# Copyright (c) 2017--, Evguenia Kopylova
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import print_function
import logging, yaml
import click
import sevenbridges as sb
from sevenbridges.errors import SbgError
from os.path import join, splitext, basename, isfile, getsize
from os import listdir
from collections import OrderedDict
import re
import inspect
import sys
import time


def load_config(yaml_fp):
    """Load CGC API configuration file.

    Parameters
    ----------
    yaml_fp: str
        Filepath to CGC API configuration file

    Return
    ------
    logger: logger instance
        Log
    """
    try:
        fp = open(yaml_fp)
        config = yaml.load(fp)
    except:
        raise SbgError('%s file missing!' % yaml_fp)

    logger = logging.getLogger('log')
    log_handler = logging.FileHandler(config['log_file'])
    formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
    log_handler.setFormatter(formatter)
    logger.addHandler(log_handler)
    logger.setLevel(logging.DEBUG)

    return logger, config


def download_links(api,
                   logger,
                   config,
                   metadata):
    """Generate download links for specified files.
    """
    logger.info('Generate download links.')
    # BAM files for disease type
    bam_list = list(
      api.files.query(
        project=config['project'],
        metadata={'disease_type': config['disease'],
                  'data_format': ['BAM'],
                  'experimental_strategy': ['RNA-Seq', 'WGS'],
                  'data_type': ['Raw sequencing data']}).all())
    # Flagstat files for disease type
    stats_list = []
    print("Number of BAM files: %s" % len(bam_list))
    metadata_set = set()
    for bam_file in bam_list:
      metadata_set.add(bam_file.metadata[metadata])
    print(metadata_set)


@click.command()
@click.option('--yaml-fp', required=True,
              type=click.Path(resolve_path=True, readable=True, exists=False,
                              file_okay=True),
              help='Filepath to API yaml file')
@click.option('--metadata', required=False, type=str, default='')
def main(yaml_fp, metadata):
    logger, config = load_config(yaml_fp)
    api = sb.Api(url=config['api-url'], token=config['token'])

    download_links(api=api, logger=logger, config=config,
                   metadata=metadata)


if __name__ == "__main__":
    main()
