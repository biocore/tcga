#!/usr/bin/env python

# ----------------------------------------------------------------------------
# Copyright (c) 2017--, Evguenia Kopylova, Jad Kanbar
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
                   download_links_fp):
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
    # Files for all disease types
    all_list = list(
      api.files.query(project=config['project']).all())
    # Flagstat files for disease type
    stats_list = []
    print("Number of BAM files: %s" % len(bam_list))
    for bam_file in bam_list:
      bam_file_bn = splitext(bam_file.name)[0]
      for all_file in all_list:
        if "%s_stats.txt" % bam_file_bn in all_file.name:
          stats_list.append(all_file)
    print('Total stats files: %s' % len(stats_list))
    with open(download_links_fp, 'w') as download_f:
      for i, _file in enumerate(stats_list):
        download_f.write("%s\n" % _file.download_info().url)
        # Rate limit
        if i % 400 == 0:
          time.sleep(120)



@click.command()
@click.option('--yaml-fp', required=True,
              type=click.Path(resolve_path=True, readable=True, exists=False,
                              file_okay=True),
              help='Filepath to API yaml file')
@click.option('--download-links-fp', required=True,
              type=click.Path(resolve_path=True, readable=True, exists=False,
                              file_okay=True),
              help='Filepath to store download links')
def main(yaml_fp, download_links_fp):
    logger, config = load_config(yaml_fp)
    api = sb.Api(url=config['api-url'], token=config['token'])

    download_links(api=api, logger=logger, config=config,
                   download_links_fp=download_links_fp)


if __name__ == "__main__":
    main()
