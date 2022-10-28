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
 
    experimentsteps = facts.ParseExperimentConfig(exp_dir)
    
    # Does the output directory exist? If not, make it
    try:
        os.makedirs(os.path.join(exp_dir, "output"))
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise

    # Print out PST info if in debug mode
    if debug_mode:
        print_experimentsteps(experimentsteps)
        # Exit
        sys.exit(0)

    # Initialize the EnTK App Manager
    amgr = AppManager(hostname=rcfg['rabbitmq']['hostname'], port=rcfg['rabbitmq']['port'], autoterminate=False)

    # Apply the resource configuration provided by the user
    amgr.resource_desc = facts.LoadResourceConfig(exp_dir)

    # Load the localization list
    if(not os.path.isfile(os.path.join(exp_dir, "location.lst"))):
        with open(os.path.join(exp_dir, "location.lst"), 'w') as templocationfile:
            templocationfile.write("New_York\t12\t40.70\t-74.01")
    amgr.shared_data = [os.path.join(exp_dir, "location.lst")]

    for pipelines in experimentsteps:
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

def print_experimentsteps(experimentsteps):

        for this_step in experimentsteps:
            print('EXPERIMENT STEP: ', experimentsteps.index(this_step))
            print('-----------------')
            print_pipeline(this_step)
            print('')


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
