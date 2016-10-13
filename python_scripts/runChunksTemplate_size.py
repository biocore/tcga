from __future__ import print_function

import argparse, unicodedata, logging, yaml
import sevenbridges as sb
from sevenbridges.errors import SbgError
import math
import pdb
import numpy as np

global config
maxInputsForScatter = 15
task_basename = "Breast_Invasive_Carcinoma_bam2fasta"

try:
    fp = open('yaml_test_config.yaml')
    config = yaml.load(fp)
except:
    raise SbgError('config.yaml file missing!')

logger = logging.getLogger('log')
log_handler = logging.FileHandler(config['log_file'])
formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
log_handler.setFormatter(formatter)
logger.addHandler(log_handler)
logger.setLevel(logging.DEBUG)

max_task_number = config['task_max_per_run']
project = config['project']
app = config['app']


def create_draft_tasks(api):
    logger.info('Creating draft tasks!')
    input_config = config['inputs']
    file_list = list(api.files.query(project=project, metadata = {'disease_type': config['disease']}).all())

    bam_inputs = [_file for _file in file_list if
                  _file.name.lower().endswith(input_config['input_bam_file'])]
    my_project  = api.projects.get(id = project)
    gb = []
    group = []
    gb_group_name = []
    bam_files_inputs = []
    counter = 0
    task_counter = 0

    for i in range(len(bam_inputs)):
        single_file = api.files.get(id = bam_inputs[i].id)

        gb = (single_file.size/float(1073741824))
        gb_group_name.append(bam_inputs[i].name)
        bam_files_inputs.append(api.files.query(project=my_project, \
                                        names =[bam_inputs[i].name])[0])
        group.append(gb)
        counter += 1

        #RUN if file size is bigger than 400GB
        if gb > 400:
            inputs = {
                        'input_bam_file': bam_files_inputs
                     }
            task_counter += 1
            name = task_basename + "_" + str(task_counter) + \
                        "_task_" + str(len(gb_group_name)) + "_files_" + str(np.sum(group))+"Gb" + \
                        "_" + config['disease']
            task_name = name

            try:
                api.tasks.create(task_name, project, app, inputs=inputs, description=task_name)
            except SbgError as e:
                logger.error("Draft task was not created!", exc_info=e)
                raise SbgError("Draft task was not created!")

            gb_group_name = []
            gb = []
            group = []
            bam_files_inputs = []

        #RUN chunks of files if sum of file sizes is greater than 400GB
        #and less thank 700GB
        if np.sum(group) > 400 and np.sum(group) < 700:
            inputs = {
                        'input_bam_file': bam_files_inputs
                     }
            task_counter += 1
            name = task_basename + "_" + str(task_counter) + \
                        "_task_" + str(len(gb_group_name)) + "_files_" + str(np.sum(group))+"Gb" + \
                        "_" + config['disease']
            task_name = name

            try:
                api.tasks.create(task_name, project, app, inputs=inputs, description=task_name)
            except SbgError as e:
                logger.error("Draft task was not created!", exc_info=e)
                raise SbgError("Draft task was not created!")

            gb_group_name = []
            gb = []
            group = []
            bam_files_inputs = []

        #RUN last remainig files that did not sum to greater than 400GB
        if counter == len(bam_inputs) and np.sum(group) > 0:

            inputs = {
                        "input_bam_file" : bam_files_inputs
                     }
            task_counter += 1
            name = task_basename + "_" + str(task_counter) + \
                        "_task_" + str(len(gb_group_name)) + "_files_" + str(np.sum(group))+"Gb" + \
                        "_" + config['disease']
            task_name = name
            try:
                api.tasks.create(task_name, my_project.id, app, inputs=inputs, description=task_name)
            except SbgError as e:
                logger.error("Draft task was not created!", exc_info=e)
                raise SbgError("Draft task was not created!")

            gb_group_name = []
            gb = []
            group = []
            bam_files_inputs = []

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


def status(api):
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


def compare_tasks(task, task2):
  if (task.app == task2.app) and task.inputs.equals(task2.inputs):
    return True
  else:
    return False

def taskInTaskList(task, taskList):
  for task2 in taskList:
      if compare_tasks(task, task2):
        return True
  return False


def rerun_failed(api):
	logger.info('Running failed tasks!')

	failed_tasks = list(api.tasks.query(project=project, limit=100, status='FAILED').all())
	completed_tasks = list(api.tasks.query(project=project, limit=100, status='COMPLETED').all())
	running_tasks = list(api.tasks.query(project=project, limit=100, status='RUNNING').all())
	draft_tasks = list(api.tasks.query(project=project, limit=100, status='DRAFT').all())
	queued_tasks = list(api.tasks.query(project=project, limit=100, status='QUEUED').all())

	for task in failed_tasks:
		if taskInTaskList(task, completed_tasks) or taskInTaskList(task, running_tasks) or taskInTaskList(task, draft_tasks) or taskInTaskList(task, queued_tasks):
			continue
		else:
			try:
				api.tasks.create(
				task.name, project, task.app, inputs=task.inputs, description=task.name
				)
			except SbgError as e:
				logger.error("Draft task was not created!", exc_info=e)
				raise SbgError("Draft task was not created!")



if __name__ == '__main__':
    print('got to main')

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "option", nargs="?", help="create|run|status|rerun_failed"
    )
    args = parser.parse_args()

    sb_config = sb.Config(url=config['api-url'], token=config['token'])
    api = sb.Api(config=sb_config)

    if args.option == 'create':
        create_draft_tasks(api)

    if args.option == 'run':
        run_tasks(api)

    if args.option == 'status':
        status(api)

    if args.option == 'rerun_failed':
        rerun_failed(api)
