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


def compare_files(api,
                  logger,
                  config,
                  dp,
                  download_links_fp):
    """Compare copied files to originals.
    """
    logger.info('Compare files.')
    file_list = list(api.files.query(project=config['project']).all())
    my_files = [_file for _file in file_list if u"Output: %s" % config['disease'] in _file.tags]
    print('Total files on CGC: %s' % len(my_files))
    with open(download_links_fp, 'w') as download_f:
      files_missing = 0
      for i, _file in enumerate(my_files):
        sys.stdout.write("%s) %s" % (i+1, _file.name))
        if isfile(join(dp, _file.name)):
          if _file.size == getsize(join(dp, _file.name)):
            sys.stdout.write("\t%s\t%s\n" % (_file.size, getsize(join(dp, _file.name))))
          else:
            sys.stdout.write("\t%s\t%s\tSIZE NOT EQUAL\n" % (_file.size, getsize(join(dp, _file.name))))
        else:
          sys.stdout.write("\tmissing locally\n")
          files_missing += 1
          download_f.write("%s\n" % _file.download_info().url)
      print('\n%s files missing' % files_missing)


@click.command()
@click.option('--yaml-fp', required=True,
              type=click.Path(resolve_path=True, readable=True, exists=False,
                              file_okay=True),
              help='Filepath to output CGC API yaml file')
@click.option('--dp', required=True,
              type=click.Path(resolve_path=True, readable=True, exists=False,
                              file_okay=True),
              help='Directory path where copied files are')
@click.option('--download-links-fp', required=True,
              type=click.Path(resolve_path=True, readable=True, exists=False,
                              file_okay=True),
              help='Filepath to store download links for missing files')
def main(yaml_fp, dp, download_links_fp):
    logger, config = load_config(yaml_fp)
    api = sb.Api(url=config['api-url'], token=config['token'])

    compare_files(api=api, logger=logger, config=config, dp=dp, download_links_fp=download_links_fp)


if __name__ == "__main__":
    main()
