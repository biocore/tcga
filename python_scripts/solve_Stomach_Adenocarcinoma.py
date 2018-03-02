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
import re
import inspect
import sys
import glob


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


def retrieve_files(api,
                   config):
    """Download info regarding WGS and RNA-Seq files for disease type.
    """
    print('Retrieve files for cancer type: %s' % config['disease'])
    # Retrieve BAM files associated with project, disease type,
    # data format experimental strategy and data type
    bam_inputs = list(
        api.files.query(
            project=config['project'],
            metadata={'disease_type': config['disease'],
                      'data_format': ['BAM'],
                      'experimental_strategy': ['RNA-Seq', 'WGS'],
                      'data_type': ['Raw sequencing data']}).all())
    # Loop through BAM files computing total size, create task if size within
    # lower and upper bounds
    num_wgs = 0
    num_rna_seq = 0
    for i, file in enumerate(bam_inputs):
        file_size_gb = file.size/float(1073741824)
        exp_str = file.metadata['experimental_strategy']
        if exp_str == 'WGS':
            num_wgs += 1
        elif exp_str == 'RNA-Seq':
            num_rna_seq += 1
        else:
            raise ValueError('%s is not supported' % exp_str)
    print("# WGS: %s" % num_wgs)
    print("# RNA-Seq: %s" % num_rna_seq)

    return bam_inputs


def sort_files(bam_inputs, data_dp):
    """Sort files.
    """
    # Group files.
    file_groups = {}
    for file in bam_inputs:
        file_basename = splitext(file.name)[0]
        all_files = glob.glob(join(data_dp, "*%s_*" % file_basename))
        if file_basename not in file_groups:
            file_groups[file_basename] = all_files
        else:
            file_groups[file_basename].extend(all_files)
    # Remove duplicates.
    for group in file_groups:
        group_set = set(file_groups[group])
        file_groups[group] = list(group_set)
    print("# Total groups: %s" % len(file_groups))
    # Split files.
    total_reads = 0
    for group in file_groups:
        file_stats = {}
        for file in file_groups[group]:
            with open(file) as f:
                first_line = f.readline()
                if "in total" in first_line:
                    total = first_line.strip().split()[0]
                    file_stats[basename(file)] = (int(total), 'all')
                    continue
                first_line = first_line.strip().split()
                if len(first_line) == 1:
                    total = first_line[0]
                    file_stats[basename(file)] = (int(total), "#")
        total_union = set()
        for stat in file_stats:
            total_union.add(file_stats[stat][0])
        if len(total_union) == 1:
            total_reads+=int(list(total_union)[0])
    print("# Total reads: %s" % total_reads)


@click.command()
@click.option('--yaml-fp', required=True,
              type=click.Path(resolve_path=True, readable=True, exists=False,
                              file_okay=True),
              help='Filepath to output CGC API yaml file')
@click.option('--data-dp', required=True,
              type=click.Path(resolve_path=True, readable=True, exists=False,
                              file_okay=True),
              help='Directory containing data files')
def main(yaml_fp,
         data_dp):
    logger, config = load_config(yaml_fp)
    api = sb.Api(url=config['api-url'], token=config['token'])

    bam_inputs = retrieve_files(api=api,
                                config=config)
    sort_files(bam_inputs=bam_inputs,
               data_dp=data_dp)


if __name__ == "__main__":
    main()