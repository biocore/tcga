#!/usr/bin/env python

#-----------------------------------------------------------------------------
# Copyright (c) 2016--, Evguenia Kopylova, Jad Kanbar, SevenBridges dev team.
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
from os.path import join


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
        raise SbgError('config.yaml file missing!')

    logger = logging.getLogger('log')
    log_handler = logging.FileHandler(config['log_file'])
    formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
    log_handler.setFormatter(formatter)
    logger.addHandler(log_handler)
    logger.setLevel(logging.DEBUG)

    return logger, config


def create_task_bam2fasta_cgc(all_files,
                              logger,
                              task_name,
                              config,
                              api):
    """Create CGC task for samtools-bam2fasta-worklfow workflow.

    Parameters
    ----------
    all_files: list
        TCGA file IDs
    logger: logger instance
        Log
    task_name: str
        CGC task name
    config: dict
        YAML configuration file
    api: SevenBridges Api instance
        Api       
    """
    inputs = {"input_bam_file" : all_files}
    task_name = "bam2fasta_%s" % task_name
    logger.info('\tName: %s' % task_name)
    my_project = api.projects.get(id = config['project'])
    #try:
    #    api.tasks.create(task_name,
    #                     my_project.id,
    #                     config['app-bam2fasta'],
    #                     inputs=inputs,
    #                     description=task_name)
    #except SbgError as e:
    #    logger.error("Draft task was not created!", exc_info=e)
    #    raise SbgError("Draft task was not created!")


def create_task_workflow_cgc(local_mapping_fp,
                             all_files,
                             task_name,
                             api,
                             config,
                             logger):
    """Create CGC task for tcga_fasta_input_disease_type_workflow workflow.

    Parameters
    ----------
    local_mapping_fp: str
        Filepath to master QIIME mapping file
    all_files: list
        TCGA file IDs
    task_name: str
        CGC task name
    api: SevenBridges Api instance
        Api
    config: dict
        YAML configuration file
    logger: logger instance
        Log

    Returns
    -------
    all_files: list
        TCGA file IDs
    total_size_gb: float
        Total size of all TCGA files
    """
    inputs = {"input_bam_file" : all_files,
              "bacterial_database_idx" : '',
              "bacterial_nodes_dmp": '',
              "bacterial_names_dmp": '',
              "bacterial_database_kdb": '',
              "viral_database_idx": '',
              "viral_names_dmp": '',
              "viral_nodes_dmp": '',
              "viral_database_kdb": '',
              "fasta_file_input": '',
              "qiime_mapping_file_1": ''}
    task_name = "workflow_%s" % task_name
    my_project = api.projects.get(id = config['project'])
    #try:
    #    api.tasks.create(task_name,
    #                     my_project.id,
    #                     config['app-workflow'],
    #                     inputs=inputs,
    #                     description=task_name)
    #except SbgError as e:
    #    logger.error("Draft task was not created!", exc_info=e)
    #    raise SbgError("Draft task was not created!")
    # Initialize files array and total size
    all_files = []
    total_size_gb = 0.0
    return all_files, total_size_gb


def generate_mapping_file(mapping_fp,
                          all_files,
                          config,
                          total_tasks_created,
                          output_dp,
                          sampleID_count):
    """Create mini mapping file based on defined sample IDs.

    Parameters
    ----------
    mapping_fp: str
        Filepath to master QIIME mapping file
    all_files: list
        List of CGC file IDs for which to generate mini-mapping file
    config: dict
        YAML configuration file
    total_tasks_created: int
        Number of task
    output_dp: str
        Output directory path
    sampleID_count: int
        Begin naming sample IDs from this integer

    Returns
    -------
    output_fp: str
        Filepath to mini-mapping file
    sampleID_count: int
        Updated sampleID count start
    """
    disease_type = config['disease'].split()
    filename = "%s_cgc_qiime_mapping_file_%s.txt" % (
        '_'.join(disease_type), total_tasks_created)
    output_fp = join(output_dp, filename)
    all_files_names = [file.name for file in all_files]
    with open(output_fp, 'w') as output_f:
        with open(mapping_fp) as mapping_f:
            for line in mapping_f:
                if line.startswith('#SampleID'):
                    output_f.write(line)
                else:
                    line = line.strip().split('\t')
                    # file name
                    filename = line[3]
                    if filename in all_files_names:
                        # update sampleID count
                        output_f.write('s%s\t' % sampleID_count)
                        sampleID_count += 1
                        output_f.write('\t'.join(line[1:]))
                        output_f.write('\n')
    return output_fp, sampleID_count


def create_tasks(api,
                 mapping_fp,
                 task_basename,
                 logger,
                 config,
                 lower_bound_group_size,
                 upper_bound_group_size,
                 output_dp,
                 count_start):
    """Create draft tasks for tcga-workflow-fasta-input-full-kraken-test
       workflow.

    Parameters
    ----------
    api: SevenBridges Api instance
        Api
    mapping_fp: str
        Filepath to master QIIME mapping file
    task_basename: str
        Basename to use for CGC task
    logger: logger instance
        Log
    config: dict
        YAML configuration file
    lower_bound_group_size: int
        Lower bound on total size of input files to pass to workflow
    upper_bound_group_size: int
        Upper bound on total size of input files to pass to workflow
    output_dp: str
        Directory path to output QIIME mini mapping files
    count_start: int
        Count from which to start SampleID generation
    """
    logger.info('Creating draft tasks.')
    input_config = config['inputs-bam2fasta']
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
    total_size_gb = 0.0
    all_files = []
    total_files_tasked = 0
    total_tasks_created = 0
    sampleID_count = count_start
    for i, file in enumerate(sorted(bam_inputs)):
        file_size_gb = file.size/float(1073741824)
        # New file will cause total file size to exceed upper limit, create
        # task and add new file to next task
        if (total_size_gb + file_size_gb > upper_bound_group_size and
                len(all_files) > 1):
            total_files_tasked += len(all_files)
            total_tasks_created += 1
            # Add info to logger
            logger.info('Task %s: %s files, %.2f Gb' % (total_tasks_created,
                                                        len(all_files),
                                                        total_size_gb))
            local_mapping_fp, sampleID_count = generate_mapping_file(
                mapping_fp, all_files, config, total_tasks_created, output_dp,
                sampleID_count)
            task_name = "%s_%s_task_%s_files_%.2f_Gb_%s" % (
                task_basename,
                str(total_tasks_created),
                str(len(all_files)),
                total_size_gb,
                '_'.join(config['disease']))
            # Create draft tasks for samtools-bam2fasta-workflow workflow
            create_task_bam2fasta_cgc(all_files, logger, task_name, config,
                                      api)
            # Create draft tasks for tcga_fasta_input_disease_type_workflow
            # workflow
            all_files, total_size_gb = create_task_workflow_cgc(
                local_mapping_fp, all_files, task_name, api, config, logger)
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
            total_files_tasked += len(all_files)
            total_tasks_created += 1
            # Add info to logger
            logger.info('Task %s: %s files, %.2f Gb' % (total_tasks_created,
                                                        len(all_files),
                                                        total_size_gb))
            local_mapping_fp, sampleID_count = generate_mapping_file(
                mapping_fp, all_files, config, total_tasks_created, output_dp,
                sampleID_count)
            task_name = "%s_%s_task_%s_files_%.2f_Gb_%s" % (
                task_basename,
                str(total_tasks_created),
                str(len(all_files)),
                total_size_gb,
                '_'.join(config['disease']))
            # Create draft tasks for samtools-bam2fasta-workflow workflow
            create_task_bam2fasta_cgc(all_files, logger, task_name, config,
                                      api)
            # Create draft tasks for tcga_fasta_input_disease_type_workflow
            # workflow
            all_files, total_size_gb = create_task_workflow_cgc(
                local_mapping_fp, all_files, task_name, api, config, logger)
    logger.info('Total tasks created: %s' % str(total_tasks_created))
    logger.info('Total files tasked: %s' % str(total_files_tasked))
    logger.info('Total files for disease type: %s' % str(len(bam_inputs)))


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
@click.option('--output-dp', required=True,
              type=click.Path(resolve_path=True, readable=True, exists=False,
                              file_okay=True),
              help='Directory path to output QIIME mini-mapping files')
@click.option('--count-start', required=True, type=int,
              help='Count from which to start SampleID generation')
def main(mapping_fp,
         yaml_fp,
         task_basename,
         create_draft_tasks,
         run_draft_tasks,
         check_status,
         lower_bound_group_size,
         upper_bound_group_size,
         output_dp,
         count_start):
    logger, config = load_config(yaml_fp)
    sb_config = sb.Config(url=config['api-url'], token=config['token'])
    api = sb.Api(config=sb_config)

    if create_draft_tasks:
        create_tasks(api, mapping_fp, task_basename, logger, config,
                     lower_bound_group_size, upper_bound_group_size,
                     output_dp, count_start)
    if run_draft_tasks:
        run_tasks(api)
    if check_status:
        show_status(api)


if __name__ == "__main__":
    main()