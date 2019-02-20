import os
import sys
import yaml
import argparse
from pprint import pprint
from radical.entk import Pipeline, Stage, Task, AppManager

def get_pre_proc_stage(wcfg):
    
    stage_cfg = wcfg['pre_process']

    # Create a stage and give it a name -- names are useful for referencing
    s = Stage()
    s.name = 'pre-process-stage'

    # Create and add Task to Stage
    t = Task()
    t.name = 'pre-process-task'

    # Pre exec let you load modules, set environment before executing the workload
    t.pre_exec = [stage_cfg['pre_exec']]

    # Executable to use for the task
    t.executable = stage_cfg['executable']

    # List of arguments for the executable
    t.arguments = stage_cfg['arguments']

    # CPU requirements for this task
    t.cpu_threads = {
                        'processes': stage_cfg['cpu']['processes'],
                        'process-type': stage_cfg['cpu']['process-type'],
                        'threads-per-process': stage_cfg['cpu']['threads-per-process'],
                        'thread-type': stage_cfg['cpu']['thread-type'],
                    }

    # Upload data from your local machine to the remote machine
    # Note: Remote machine can be the local machine
    t.upload_input_data = stage_cfg['upload_input_data']

    # Add task to stage
    s.add_tasks(t)

    return s

def get_fit_stage(wcfg, pipe_name):
    
    stage_cfg = wcfg['fit']

    # Create a stage and give it a name -- names are useful for referencing
    s = Stage()
    s.name = 'fit-stage'

    # Create and add Task to Stage
    t = Task()
    t.name = 'fit-task'

    # Pre exec let you load modules, set environment before executing the workload
    t.pre_exec = [stage_cfg['pre_exec']]

    # Executable to use for the task
    t.executable = stage_cfg['executable']

    # List of arguments for the executable
    t.arguments = stage_cfg['arguments']

    # CPU requirements for this task
    t.cpu_threads = {
                        'processes': stage_cfg['cpu']['processes'],
                        'process-type': stage_cfg['cpu']['process-type'],
                        'threads-per-process': stage_cfg['cpu']['threads-per-process'],
                        'thread-type': stage_cfg['cpu']['thread-type'],
                    }

    # Upload data from your local machine to the remote machine
    # Note: Remote machine can be the local machine
    t.upload_input_data = stage_cfg['upload_input_data']

    # Copy output of preprocess stage to the folder where this task executes
    loc = '$Pipeline_%s_Stage_pre-process-stage_Task_pre-process-task'%pipe_name
    t.copy_input_data = ['%s/%s'%(loc, x) for x in stage_cfg['copy_input_data']]

    # Add task to stage
    s.add_tasks(t)

    return s

def get_projection_stage(wcfg, pipe_name):
    
    stage_cfg = wcfg['project']

    # Create a stage and give it a name -- names are useful for referencing
    s = Stage()
    s.name = 'project-stage'

    # Create and add Task to Stage
    t = Task()
    t.name = 'project-task'

    # Pre exec let you load modules, set environment before executing the workload
    t.pre_exec = [stage_cfg['pre_exec']]

    # Executable to use for the task
    t.executable = stage_cfg['executable']

    # List of arguments for the executable
    t.arguments = stage_cfg['arguments']

    # CPU requirements for this task
    t.cpu_threads = {
                        'processes': stage_cfg['cpu']['processes'],
                        'process-type': stage_cfg['cpu']['process-type'],
                        'threads-per-process': stage_cfg['cpu']['threads-per-process'],
                        'thread-type': stage_cfg['cpu']['thread-type'],
                    }

    # Upload data from your local machine to the remote machine
    # Note: Remote machine can be the local machine
    t.upload_input_data = stage_cfg['upload_input_data']

    # Copy output of fit and preprocess stage to the folder where this task executes
    pre_proc_loc = '$Pipeline_%s_Stage_pre-process-stage_Task_pre-process-task'%pipe_name
    fit_loc = '$Pipeline_%s_Stage_fit-stage_Task_fit-task'%pipe_name
    t.copy_input_data = ['%s/%s'%(fit_loc, stage_cfg['copy_input_data'][0]),
                         '%s/%s'%(pre_proc_loc, stage_cfg['copy_input_data'][1])]

    # Add task to stage
    s.add_tasks(t)

    return s

def get_post_proc_stage(wcfg):
    return None


def create_pipelines(wcfg):

    p = Pipeline()
    p.name = 'kopp'
    
    # Create and add pre-process, fit, project, and post-process stages if
    # they are defined


    pre_proc_stage = get_pre_proc_stage(wcfg)
    if pre_proc_stage:
        p.add_stages(pre_proc_stage)

    fit_stage = get_fit_stage(wcfg, p.name)
    if fit_stage:
        p.add_stages(fit_stage)

    projection_stage = get_projection_stage(wcfg, p.name)
    if projection_stage:
        p.add_stages(projection_stage)

    post_proc_stage = get_post_proc_stage(wcfg)
    if post_proc_stage:
        p.add_stages(post_proc_stage)

    return p

    

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Process some arguments to get resource and workflow cfgs')
    parser.add_argument('--wcfg', help='path to workflow cfg file', required=False, default='./workflow_cfg.yml')
    parser.add_argument('--rcfg', help='path to resource cfg file', required=False, default='./resource_cfg.yml')

    args = parser.parse_args()
    if not os.path.isfile(args.wcfg):
        print '%s does not exist' % args.wcfg
        sys.exit(1)

    if not os.path.isfile(args.rcfg):
        print '%s does not exist' % args.rcfg
        sys.exit(1)

    with open(args.rcfg, 'r') as fp:
        rcfg = yaml.load(fp)

    with open(args.wcfg, 'r') as fp:
        wcfg = yaml.load(fp)

    pipelines = create_pipelines(wcfg)
    if not isinstance(pipelines, list):
        pipelines = [pipelines]

    # If you don't want to run it yet, but just take a look at all the stages
    # in the pipeline and the tasks in all stages, set debug_mode in the
    # workflow_cfg.yml to True. False to run.
    if wcfg['debug_mode']:
        for p in pipelines:
            print 'Pipeline %s:'%p.name
            print '############################'
            pprint(p.to_dict())
            for s in p.stages:
                print 'Stage %s:'%s.name
                print '==============================='
                pprint(s.to_dict())
                for t in s.tasks:
                    print 'Task %s:'%t.name
                    print '--------------------------------'
                    pprint(t.to_dict())

        sys.exit(0)

    amgr = AppManager(hostname=rcfg['rabbitmq']['hostname'],
                      port=rcfg['rabbitmq']['port'])

    res_desc = {
        'resource': rcfg['resource-desc']['name'],
        'walltime': rcfg['resource-desc']['walltime'],
        'cpus': rcfg['resource-desc']['cpus'],
        'queue': rcfg['resource-desc']['queue'],
        'project': rcfg['resource-desc']['project']
    }

    amgr.resource_desc = res_desc    

    amgr.workflow = pipelines
    amgr.run()