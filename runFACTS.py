import sys
import os
import argparse
import yaml
import re
import errno
import time
from pprint import pprint
import FACTS as facts
from radical.entk import Pipeline, Stage, Task, AppManager

def run_experiment(exp_dir, debug_mode, no_total_flag):

    # Initialize a list for experiment steps (each step being a set of pipelines)
    experimentsteps = []
    currentstep = 0

    # Define the configuration and resource file names
    rfile = os.path.join(exp_dir, "resource.yml")
    cfile = os.path.join(exp_dir, "config.yml")

    # Does the experiment configuration file exist?
    if not os.path.isfile(cfile):
        print('{} does not exist'.format(cfile))
        sys.exit(1)

    # Does the resource file exist?
    if not os.path.isfile(rfile):
        print('{} does not exist'.format(rfile))
        sys.exit(1)

    # Load the resource and experiment configuration files
    with open(rfile, 'r') as fp:
        rcfg = yaml.safe_load(fp)
    with open(cfile, 'r') as fp:
        ecfg = yaml.safe_load(fp)

    # Does the output directory exist? If not, make it
    try:
        os.makedirs(os.path.join(exp_dir, "output"))
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise

    # Reserved configuration entries
    reserved_econfig_entries = ["global-options", "total-options", "extremesealevel-options"]

    # Are there global options?
    if "global-options" in ecfg.keys():
        global_options = ecfg["global-options"]
    else:
        global_options = {}


    # Initialize a list for pipelines
    pipelines = []

    # Loop through the user-requested modules
    for this_mod in ecfg.keys():

        # Skip this entry if it's not associated with SLR projection workflow
        if this_mod in reserved_econfig_entries:
            continue

        # if this request is a collection of modules
        if "multistep" in ecfg[this_mod].keys():

            # if already have accumulated pipelines, load them as step and reset
            if len(pipelines) > 0:
                experimentsteps.append(pipelines)
                currentstep += 1
                pipelines = []

            for this_mod_sub in ecfg[this_mod].keys():
                parsed = ParsePipelineConfig(this_mod_sub, ecfg[this_mod][this_mod_sub], global_options=global_options)
                pipelines.append(facts.GeneratePipeline(parsed['pcfg'], parsed['modcfg'], parsed['pipe_name'], exp_dir))
               
        else:
            parsed = ParsePipelineConfig(this_mod, ecfg[this_mod], global_options=global_options)
            pipelines.append(facts.GeneratePipeline(parsed['pcfg'], parsed['modcfg'], parsed['pipe_name'], exp_dir))

    experimentsteps.append(pipelines)

    # Print out PST info if in debug mode
    if debug_mode:
        for this_step in experimentsteps:
            print('EXPERIMENT STEP: ', experimentsteps.index(this_step))
            print('-----------------')
            print_pipeline(this_step)
            print('')

        # Exit
        sys.exit(0)

    # Initialize the EnTK App Manager
    amgr = AppManager(hostname=rcfg['rabbitmq']['hostname'], port=rcfg['rabbitmq']['port'], autoterminate=False)

    # Apply the resource configuration provided by the user
    res_desc = {'resource': rcfg['resource-desc']['name'],
        'walltime': rcfg['resource-desc']['walltime'],
        'cpus': rcfg['resource-desc']['cpus'],
        'queue': rcfg['resource-desc']['queue'],
        'project': rcfg['resource-desc']['project']}
    amgr.resource_desc = res_desc

    # Load the localization list
    if(not os.path.isfile(os.path.join(exp_dir, "location.lst"))):
        with open(os.path.join(exp_dir, "location.lst"), 'w') as templocationfile:
            templocationfile.write("New_York\t12\t40.70\t-74.01")
    amgr.shared_data = [os.path.join(exp_dir, "location.lst")]


    for this_step in experimentsteps:
        print('EXPERIMENT STEP: ', experimentsteps.index(this_step))
 
        # Assign the list of pipelines to the workflow
        amgr.workflow = pipelines

        # Run the SLR projection workflow
        amgr.run()


    # the extreme sea level and totaling workflows are hard coded in right now, 
    # predating implementation of functionality to allow nesting of modules
    # they should be converted over at some point

    # If the user want to run extremesealevel module, must perform a total and flag
    do_extremesl_flag = False
    if "extremesealevel-options" in ecfg.keys():
        no_total_flag = False
        do_extremesl_flag = True

    # Setup and run the totaling workflow
    if not no_total_flag:
        ecfg['total-options']["options"].update(ecfg['global-options'])
        total_pipeline = facts.GenerateTotalPipeline(ecfg['total-options'], exp_dir)
        amgr.workflow = [total_pipeline]
        amgr.run()

    # Setup the extreme sea-level workflow if needed
    if do_extremesl_flag:
        this_mod = "extremesealevel-options"
        pcfg_file = os.path.join(os.path.dirname(__file__), "modules", ecfg[this_mod]['module_set'], ecfg[this_mod]['module'], "pipeline.yml")
        if not os.path.isfile(pcfg_file):
            print('{} does not exist'.format(pcfg_file))
            sys.exit(1)
        with open(pcfg_file, 'r') as fp:
            pcfg = yaml.safe_load(fp)

        # Append the global options to this module
        ecfg[this_mod]["options"].update(global_options)

        # Generate a pipeline for this module
        pipe_name = ".".join((ecfg[this_mod]['module_set'], ecfg[this_mod]['module']))
        esl_pipeline = facts.GeneratePipeline(pcfg, ecfg[this_mod], pipe_name, exp_dir)
        amgr.workflow = [esl_pipeline]
        amgr.run()

    # Close the application manager
    amgr.terminate()


    return(None)

def ParsePipelineConfig(this_mod, modcfg, global_options={}):
    # Load the pipeline configuration file for this module
    pcfg_file = os.path.join(os.path.dirname(__file__), "modules", modcfg['module_set'], modcfg['module'], "pipeline.yml")
    if not os.path.isfile(pcfg_file):
        print('{} does not exist'.format(pcfg_file))
        sys.exit(1)
    with open(pcfg_file, 'r') as fp:
        pcfg = yaml.safe_load(fp)

    # Append the global options to this module
    modcfg["options"].update(global_options)

    # Generate a pipeline for this module
    pipe_name = ".".join((this_mod, modcfg['module_set'], modcfg['module']))
    p = {
        "modlabel": this_mod,
        "pcfg": pcfg,
        "modcfg": modcfg,
        "pipe_name": pipe_name
    }
    return p


def print_pipeline(pipelines):
    
    for p in pipelines:
        print("Pipeline {}:".format(p.name))
        print("################################")
        print(p.as_dict())
        for s in p.stages:
            print("Stage {}:".format(s.name))
            print("============================")
            pprint(s.as_dict())
            for t in s.tasks:
                print("Task {}:".format(t.name))
                print("----------------------------")
                pprint(t.as_dict())



if __name__ == "__main__":

    # Initialize the argument parser
    parser = argparse.ArgumentParser(description="The Framework for Assessing Changes To Sea-level (FACTS)")

    # Add arguments for the resource and experiment configuration files
    parser.add_argument('edir', help="Experiment Directory")
    parser.add_argument('--no-total', help="Disable totaling of global and local sea-level projection wokflow", action="store_true")
    parser.add_argument('--debug', help="Enable debug mode", action="store_true")

    # Parse the arguments
    args = parser.parse_args()

    # Does the experiment directory exist?
    if not os.path.isdir(args.edir):
        print('%s does not exist'.format(args.edir))
        sys.exit(1)

    # Go ahead and try to run the experiment
    run_experiment(args.edir, args.debug, args.no_total)


    #sys.exit(0)
