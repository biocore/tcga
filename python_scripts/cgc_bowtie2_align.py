#!/usr/bin/env python

#-----------------------------------------------------------------------------
# Copyright (c) 2016--, Evguenia Kopylova.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

"""
Remove false positive matches from Kraken output by aligning reads with Bowtie2.

Algorithm:
	1. Retrieve taxonomy IDs from Kraken output and use those to identify
		corresponding genomes in 54K reference database
	2. Pass comma-separated list of genomes to bowtie2-build
"""

import click

@click.command()
@click.option('--kraken-report-fp', required=True,
              type=click.Path(resolve_path=True, readable=True, exists=True,
                              file_okay=True),
              help='Filepath to Kraken report')
@click.option('--yaml-fp', required=True,
              type=click.Path(resolve_path=True, readable=True, exists=False,
                              file_okay=True),
              help='Filepath to output CGC API yaml file')

def main(kraken_mpa_report_fp):
    logger, config = load_config(yaml_fp)
    sb_config = sb.Config(url=config['api-url'], token=config['token'])
    api = sb.Api(config=sb_config)

    if create_draft_tasks:
        create_tasks(api, mapping_fp, logger, config, lower_bound_group_size,
                     upper_bound_group_size, output_dp, count_start)
    elif run_draft_tasks:
        run_tasks(api, logger, config)
    elif check_status:
        show_status(api)
    else:
        raise ValueError('Please select one of --create-draft-tasks, '
                         '--run-draft-tasks or --check-status')


if __name__ == "__main__":
    main()