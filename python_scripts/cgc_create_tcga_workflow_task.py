#!/usr/bin/env python

#-----------------------------------------------------------------------------
# Copyright (c) 2016--, Evguenia Kopylova, Jad Kanbar, Erik Lehnert.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

"""
Create tasks for tcga-workflow-fasta-input-full-kraken-test workflow.
"""

from __future__ import print_function
import logging, yaml
import click
import sevenbridges as sb
from sevenbridges.errors import SbgError


def load_config(yaml_fp):
    """Load CGC API configuration file.

    Parameters
    ----------
    yaml_fp: str
        Filepath to CGC API configuration file

    Return
    ------
    logger
    """
    try:
        fp = open(yaml_fp)
        config = yaml.load(fp)
    except:
        raise SbgError('config.yaml file missing!')

    logger = logging.getLogger('log')
    log_handler = logging.FileHandler(config['log_file'])
    formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
    log_handler.setFormatter(formatter)
    logger.addHandler(log_handler)
    logger.setLevel(logging.DEBUG)

    return logger, log_handler, config


def create_task_cgc(local_mapping_fp,
                    all_files,
                    total_size_gb,
                    total_tasks_created,
                    task_basename,
                    api,
                    config,
                    logger):
    """Create CGC task.

    Parameters
    ----------
    local_mapping_fp: str
        Filepath to master QIIME mapping file
    all_files: list
        TCGA file IDs
    total_size_gb: float
        Total size of all TCGA files
    total_tasks_created: int
        Number of task
    task_basename: str
        CGC task basename
    api:
    config:
    logger:

    Returns
    -------
    all_files: list
        TCGA file IDs
    total_size_gb: float
        Total size of all TCGA files
    """
    inputs = {"input_bam_file" : all_files,
              "bacterial_database_idx" : ,
              "bacterial_nodes_dmp": ,
              "bacterial_names_dmp": ,
              "bacterial_database_kdb": ,
              "viral_database_idx": ,
              "viral_names_dmp": ,
              "viral_nodes_dmp": ,
              "viral_database_kdb": ,
              "fasta_file_input": ,
              "qiime_mapping_file_1": }
    task_name = "%s_%s_task_%s_files_%s_Gb_%s" % (task_basename,
                                                  str(total_tasks_created),
                                                  str(len(all_files)),
                                                  str(total_size_gb),
                                                  config['disease'])
    task_name = name
    my_project = api.projects.get(id = config['project'])
    try:
        api.tasks.create(task_name,
                         my_project.id,
                         config['app'],
                         inputs=inputs,
                         description=task_name)
    except SbgError as e:
        logger.error("Draft task was not created!", exc_info=e)
        raise SbgError("Draft task was not created!")
    # Initialize files array and total size
    all_files = []
    total_size_gb = 0.0
    return all_files, total_size_gb


def generate_mapping_file(mapping_fp, all_files):
    """Create mini mapping file based on defined sample IDs.

    Parameters
    ----------
    mapping_fp: str
        Filepath to master QIIME mapping file
    all_files: list
        List of CGC file IDs for which to generate mini-mapping file
    """
    pass

def create_tasks(api,
                 mapping_fp,
                 task_basename,
                 logger,
                 log_handler,
                 config,
                 lower_bound_group_size,
                 upper_bound_group_size):
    """Create draft tasks for tcga-workflow-fasta-input-full-kraken-test
       workflow.

    Parameters
    ----------
    api:

    mapping_fp: str
        Filepath to master QIIME mapping file
    task_basename: str
        Basename to use for CGC task
    logger:

    log_handler:

    config:

    lower_bound_group_size: int
        Lower bound on total size of input files to pass to workflow
    upper_bound_group_size: int
        Upper bound on total size of input files to pass to workflow
    """
    logger.info('Creating draft tasks!')
    input_config = config['inputs']
    # Retrieve all files associated with project and disease type
    file_list = list(
        api.files.query(
            project=config['project'],
            metadata = {'disease_type': config['disease']}).all())
    # Keep only BAM files
    bam_inputs = [_file for _file in file_list if
                  _file.name.lower().endswith(input_config['input_bam_file'])]
    # Loop through BAM files computing total size, create task if size within
    # lower and upper bounds
    my_project  = 
    total_size_gb = 0.0
    all_files = []
    total_files_tasked = 0
    total_tasks_created = 0
    for i, file in enumerate(bam_inputs):
        file_size_gb = file.size/float(1073741824)
        # New file will cause total file size to exceed upper limit, create
        # task and add new file to next task
        if (total_size_gb + file_size_gb > upper_bound_group_size and
                len(all_files) > 1):
            total_files_tasked += len(all_files)
            local_mapping_fp = generate_mapping_file(mapping_fp, all_files)
            total_tasks_created += 1
            all_files, total_size_gb = create_task_cgc(
                local_mapping_fp, all_files, total_size_gb,
                total_tasks_created, task_basename, api, config, logger)
            # Add new file to next task
            all_files.append(file)
            total_size_gb += file_size_gb
            continue
        # Add new file to next task
        all_files.append(file)
        total_size_gb += file_size_gb
        # If:
        # (1) Single file larger than upper bound limit, or
        # (2) Group of files fall within defined limit, or
        # (3) Last file encountered and group size less than defined lower
        #     bound, then
        # Create task.
        if ( (len(all_files) == 1 and
                total_size_gb >= upper_bound_group_size) or
                (total_size_gb > lower_bound_group_size and
                total_size_gb < upper_bound_group_size) or
                (i+1 == len(bam_inputs) and
                total_size_gb <= lower_bound_group_size)):
            local_mapping_fp = generate_mapping_file(mapping_fp, all_files)
            total_tasks_created += 1
            all_files, total_size_gb = create_task_cgc(
                local_mapping_fp, all_files, total_size_gb,
                total_tasks_created, task_basename, api, config, logger)
            total_files_tasked += len(all_files)
    logger.info('Total tasks created: %s' % str(total_tasks))
    logger.info('Total files tasked: %s' % str(total_files_tasked))


def run_tasks(api):
    logger.info('Running tasks!')

    running_tasks = list(
        api.tasks.query(project=project, limit=100, status='RUNNING').all()
    )
    queued_tasks = list(
        api.tasks.query(project=project, limit=100, status='QUEUED').all()
    )
    if len(running_tasks) + len(queued_tasks) >= max_task_number:
        logger.info("Maximum number of active tasks reached!")
        raise SbgError(
            'Unable to run! You already have {active} active tasks. '
            'Please try later!'.format
            (active=len(running_tasks) + len(queued_tasks)))

    draft_tasks = list(
        api.tasks.query(project=project, limit=100, status='DRAFT').all()
    )
    if len(draft_tasks) == 0:
        print('No draft tasks left to be run!')
        return

    executable_tasks = draft_tasks[0:max_task_number - len(running_tasks)]
    for task in executable_tasks:
        try:
            task.run()
        except SbgError as e:
            logger.error("Task was not started! Error happened ", exc_info=e)
            raise SbgError('Task was not started! Error happened')
        if task.status == 'DRAFT':
            logger.error("Task was not started! Task state is DRAFT!")
            raise SbgError("Task was not started! Task state is DRAFT!")


def show_status(api):
    logger.info('Fetching task statuses!')
    queued = api.tasks.query(project=project, status='QUEUED').total
    running = api.tasks.query(project=project, status='RUNNING').total
    completed = api.tasks.query(project=project, status='COMPLETED').total
    draft = api.tasks.query(project=project, status='DRAFT').total
    failed = api.tasks.query(project=project, status='FAILED').total
    aborted = api.tasks.query(project=project, status='ABORTED').total
    print("Draft={}, Queued={}, Running={}, Completed={},"
          " Failed={}, Aborted={} ".format(draft, queued,
                                           running, completed,
                                           failed, aborted)
          )


@click.command()
@click.option('--mapping-fp', required=True,
              type=click.Path(resolve_path=True, readable=True, exists=True,
                              file_okay=True),
              help='Filepath to QIIME mapping file')
@click.option('--yaml-fp', required=True,
              type=click.Path(resolve_path=True, readable=True, exists=False,
                              file_okay=True),
              help='Filepath to output CGC API yaml file')
@click.option('--task-basename', required=True, type=str,
              help='CGC task basename')
@click.option('--create-draft-tasks', required=True, type=bool, default=True,
              show_default=True, help='Create CGC draft tasks')
@click.option('--run-draft-tasks', required=False, type=bool, default=False,
              show_default=False, help='Run CGC draft tasks')
@click.option('--check-status', required=False, type=bool, default=False,
              show_default=True, help='Show CGC task status')
@click.option('--lower-bound-group-size', required=False, type=int,
              default=400, show_default=True,
              help='Lower bound on total size of input files to pass to '
              'workflow')
@click.option('--upper-bound-group-size', required=False, type=int,
              default=700, show_default=True,
              help='Upper bound on total size of input files to pass to '
              'workflow')
def main(mapping_fp,
         yaml_fp,
         task_basename,
         create_draft_tasks,
         run_draft_tasks,
         check_status,
         lower_bound_group_size,
         upper_bound_group_size):
    logger, log_handler, config = load_config(yaml_fp)
    sb_config = sb.Config(url=config['api-url'], token=config['token'])
    api = sb.Api(config=sb_config)

    if create_draft_tasks:
        create_tasks(api, mapping_fp, task_basename, logger,
                     log_handler, config, lower_bound_group_size,
                     upper_bound_group_size)
    if run_draft_tasks:
        run_tasks(api)
    if check_status:
        show_status(api)


if __name__ == "__main__":
    main()