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
from os.path import join, splitext, basename
from collections import OrderedDict
import re


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
    api: SevenBridges API instance
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
    project = config['project']
    # Upload local mapping file to project
    try:
        api.files.upload(project=project, path=local_mapping_fp)
    # File already exists
    except SbgError as e:
        logger.error(
            "Could not upload file, trying to query for it", exc_info=e)
        pass
    # Retrieve File object for mapping file
    local_mapping_file = list(
        api.files.query(
            project=project, names=[basename(local_mapping_fp)]).all())
    if len(local_mapping_file) > 1:
        raise ValueError(
            'List contains >1 files: %s' % len(local_mapping_file))
    # Retrieve File objects for all bacterial and viral database files.
    # We're not calling files directly by their ID because this can change,
    # whereas file names are expected to stay the same.
    input_index_files = list(api.files.query(
        project=project,
        names=['bacterial_database.idx',
               'bacterial_nodes.dmp',
               'bacterial_names.dmp',
               'bacterial_database.kdb',
               'viral_database.idx',
               'viral_names.dmp',
               'viral_nodes.dmp',
               'viral_database.kdb']).all())
    inputs = {}
    inputs['qiime_mapping_file'] = local_mapping_file[0]
    inputs['fasta_file_input'] = all_files
    for _file in input_index_files:
        name = _file.name
        if name == 'bacterial_database.idx':
            inputs['bacterial_database_idx'] = _file
        elif name == 'bacterial_nodes.dmp':
            inputs['bacterial_nodes_dmp'] = _file
        elif name == 'bacterial_names.dmp':
            inputs['bacterial_names_dmp'] = _file
        elif name == 'bacterial_database.kdb':
            inputs['bacterial_database_kdb'] = _file
        elif name == 'viral_database.idx':
            inputs['viral_database_idx'] = _file
        elif name == 'viral_names.dmp':
            inputs['viral_names_dmp'] = _file
        elif name == 'viral_nodes.dmp':
            inputs['viral_nodes_dmp'] = _file
        elif name == 'viral_database.kdb':
            inputs['viral_database_kdb'] = _file
        else:
            raise ValueError(
                "File %s not assigned to any input argument." % name)
    task_name = "workflow_%s" % task_name
    my_project = api.projects.get(id = config['project'])
    #try:
    #    api.tasks.create(name=task_name,
    #                     project=my_project.id,
    #                     app=config['app-workflow'],
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
                          sampleID_count,
                          logger,
                          fasta_files_workflow):
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
    logger: logger instance
        Log
    fasta_files_workflow: list
        FASTA file names

    Returns
    -------
    output_fp: str
        Filepath to mini-mapping file
    sampleID_count: int
        Updated sampleID count start
    all_files: list
        List of updated CGC file IDs (duplicates removed)
    """
    disease_type = config['disease'].split()
    filename = "%s_cgc_qiime_mapping_file_%s.txt" % (
        '_'.join(disease_type), total_tasks_created)
    output_fp = join(output_dp, filename)
    all_files_names = [file.name for file in all_files]
    all_files_names_added = []
    with open(output_fp, 'w') as output_f:
        with open(mapping_fp) as mapping_f:
            for line in mapping_f:
                if line.startswith('#SampleID'):
                    output_f.write(line)
                else:
                    line = line.strip().split('\t')
                    # FASTA file name
                    filename = line[4]
                    if filename in all_files_names:
                        # update sampleID count
                        output_f.write('s%s\t' % sampleID_count)
                        sampleID_count += 1
                        output_f.write('\t'.join(line[1:]))
                        output_f.write('\n')
                        all_files_names_added.append(filename)
    files_not_added = set(all_files_names) - set(all_files_names_added)
    all_files_updated = list(all_files)
    # At least one FASTA file analyzed not found in mapping file
    if len(files_not_added) > 0:
        logger.error(
            'Following files missing in mapping file:\n')
        # Check missing files are duplicates of those that have been added
        files_accounted_for = 0
        for _file in files_not_added:
            # Remove prefix _*_ which signifies duplicate
            regex = re.compile('_._')
            prefix = _file[0:3]
            if re.match(regex, prefix):
                original_file_name = _file[3:]
                if original_file_name not in fasta_files_workflow:
                    logger.error('\t%s, [status] missing file')
                else:
                    files_accounted_for += 1
                    logger.info('\t%s, [status] duplicate' % _file)
            # File does not have prefix _*_ which signifies it is not a
            # duplicate
            else:
                logger.error('\t%s [status] missing file' % _file)
        # Reassure user all missing files are duplicates
        if files_accounted_for == len(files_not_added):
            logger.info('All missing files are duplicates and original files '
                        'exists in analysis.')
            # Remove duplicate FASTA files from analysis
            for __file in all_files:
                if __file.name in files_not_added:
                    all_files_updated.remove(__file)
        # Non-duplicate files do not exist in mapping file, raise error
        else:
            logger.error('Not all missing files have been accounted for.')
            raise ValueError('Not all missing files have been accounted for.')
    return output_fp, sampleID_count, all_files_updated


def create_tasks(api,
                 mapping_fp,
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
    # Retrieve all files associated with project and disease type
    file_list = list(
        api.files.query(
            project=config['project'],
            metadata = {'disease_type': config['disease']}).all())
    # BAM files
    bam_inputs = [_file.name for _file in file_list if
                  _file.name.lower().endswith('bam')]
    bam_inputs_derep = list(bam_inputs)
    # Remove duplicates from analysis
    regex = re.compile('_._')
    for _file in bam_inputs:
        prefix = _file.name[0:3]
        # File is duplicate, check original file exists
        if re.match(regex, prefix):
            original_file_name = _file[3:]
            if original_file_name in bam_inputs:
                bam_inputs_derep.remove(_file)
            else:
                logger.info('%s does not have the original file in list')
    if len(bam_inputs_dere) < len(bam_inputs):
        logger.info('%s duplicate files removed' % 
                    len(bam_inputs) - len(bam_inputs_derep))
    # FASTA files
    fasta_files = {}
    for _file in file_list:
        if _file.name.lower().endswith('fasta'):
            if _file.name not in fasta_files:
                fasta_files[_file.name] = _file
            else:
                raise ValueError('%s already exists' % _file.name)
    # Check BAM associated FASTA file exists
    for _file in bam_inputs:
        file_name, file_ext = splitext(_file)
        if "%s.fasta" % file_name not in fasta_files:
            raise ValueError(
                '%s.fasta is missing from FASTA files' % file_name)
    # Remove all non BAM associated FASTA files from further analysis
    fasta_files_workflow = OrderedDict(fasta_files)
    for key, value in fasta_files.iteritems():
        file_name, file_ext = splitext(key)
        if "%s.bam" % file_name not in bam_inputs:
            del fasta_files_workflow[key]
    # Check number of BAM files equals to number of bam2fasta FASTA files
    if len(fasta_files_workflow) != len(bam_inputs):
        raise ValueError('%s != %s' % (
            len(fasta_files_workflow), len(bam_inputs)))
    # Loop through FASTA files computing total size, create task if size
    # within lower and upper bounds
    total_size_gb = 0.0
    all_files = []
    total_files_tasked = 0
    total_tasks_created = 0
    sampleID_count = count_start
    for i, key in enumerate(fasta_files_workflow):
        file = fasta_files_workflow[key]
        file_size_gb = file.size/float(1073741824)
        # If:
        # (1) File will cause total file size to exceed upper limit, then
        # Create task and add file to next task
        if (total_size_gb + file_size_gb > upper_bound_group_size and
                len(all_files) > 1):
            # Create local mapping file
            local_mapping_fp, sampleID_count, all_files =\
                generate_mapping_file(
                    mapping_fp, all_files, config, total_tasks_created,
                    output_dp, sampleID_count, logger,
                    fasta_files_workflow.keys())
            total_files_tasked += len(all_files)
            total_tasks_created += 1
            task_name = "%s_%s_task_%s_files_%.2fGb" % (
                config['disease'],
                str(total_tasks_created),
                str(len(all_files)),
                total_size_gb)
            # Add info to logger
            logger.info('Task %s: %s files, %.2f Gb' % (total_tasks_created,
                                                        len(all_files),
                                                        total_size_gb))
            # Create draft tasks for tcga_fasta_input_disease_type_workflow
            # workflow
            all_files, total_size_gb = create_task_workflow_cgc(
                local_mapping_fp, all_files, task_name, api, config, logger)
        # Add new file to next task
        all_files.append(file)
        total_size_gb += file_size_gb
        # If:
        # (1) Single file larger than upper bound limit, or
        # (2) Group of files fall within defined limit, or
        # (3) Last file encountered, then
        # Create task.
        if ( (len(all_files) == 1 and
                total_size_gb >= upper_bound_group_size) or
                (total_size_gb > lower_bound_group_size and
                total_size_gb < upper_bound_group_size) or
                i+1 == len(bam_inputs) ):
            # Create local mapping file
            local_mapping_fp, sampleID_count, all_files =\
                generate_mapping_file(
                    mapping_fp, all_files, config, total_tasks_created,
                    output_dp, sampleID_count, logger,
                    fasta_files_workflow.keys())
            total_files_tasked += len(all_files)
            total_tasks_created += 1
            task_name = "%s_%s_task_%s_files_%.2fGb" % (
                config['disease'],
                str(total_tasks_created),
                str(len(all_files)),
                total_size_gb)
            # Add info to logger
            logger.info('Task %s: %s files, %.2f Gb' % (total_tasks_created,
                                                        len(all_files),
                                                        total_size_gb))
            # Create draft tasks for tcga_fasta_input_disease_type_workflow
            # workflow
            all_files, total_size_gb = create_task_workflow_cgc(
                local_mapping_fp, all_files, task_name, api, config, logger)
    logger.info('Total tasks created: %s' % str(total_tasks_created))
    logger.info('Total files tasked: %s' % str(total_files_tasked))
    logger.info('Total files for disease type: %s' % str(len(bam_inputs)))


def run_tasks(api,
              logger,
              config):
    """Run CGC tasks.

    Parameters
    ----------
    api: SevenBridges Api instance
        Api
    mapping_fp: str
        Filepath to master QIIME mapping file
    logger: logger instance
        Log    
    """
    logger.info('Running tasks!')
    project = config['project']
    max_task_number = config['task_max_per_run']
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
    """Show CGC status.

    Parameters
    ----------
    api: SevenBridges API instance
        Api
    """
    logger.info('Fetching task statuses!')
    project = config['project']
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
