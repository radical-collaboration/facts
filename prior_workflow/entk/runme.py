import os
import sys
import yaml
import argparse
from pprint import pprint
from radical.entk import Pipeline, Stage, Task, AppManager

def get_pre_proc_stage(wcfg):
    return None

def get_calibration_stage(wcfg):
    
    stage_cfg = wcfg['calibration']

    s = Stage()
    s.name = 'calibration-stage'
    t = Task()
    t.name = 'calibration-task'
    t.pre_exec = [stage_cfg['pre_exec']]
    t.executable = stage_cfg['executable']
    t.arguments = [stage_cfg['arguments']['filename']]
    t.cpu_threads = {
                        'processes': stage_cfg['cpu']['processes'],
                        'process-type': stage_cfg['cpu']['process-type'],
                        'threads-per-process': stage_cfg['cpu']['threads-per-process'],
                        'thread-type': stage_cfg['cpu']['thread-type'],
                    }

    # Assumes run_icesheets.py is available in the current directory
    t.upload_input_data = [stage_cfg['arguments']['filename']]

    s.add_tasks(t)

    return s

def get_projection_stage(wcfg):
    return None

def get_post_proc_stage(wcfg):
    return None


def create_pipelines(wcfg):

    p = Pipeline()
    p.name = 'kopp'
    
    pre_proc_stage = get_pre_proc_stage(wcfg)
    if pre_proc_stage:
        p.add_stages(pre_proc_stage)

    calibration_stage = get_calibration_stage(wcfg)
    if calibration_stage:
        p.add_stages(calibration_stage)

    projection_stage = get_projection_stage(wcfg)
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

    pipelines = create_pipelines(wcfg)
    if not isinstance(pipelines, list):
        pipelines = [pipelines]

    amgr.workflow = pipelines

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
    else:
        amgr.run()