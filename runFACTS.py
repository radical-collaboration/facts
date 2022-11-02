import sys
import os
import argparse
import errno
from pprint import pprint
import FACTS as facts
from radical.entk import Pipeline, Stage, Task


def run_experiment(exp_dir, debug_mode):

    expconfig = facts.ParseExperimentConfig(exp_dir)
    experimentsteps = expconfig['experimentsteps']
    # ecfg = expconfig['ecfg']

    # Print out PST info if in debug mode
    if debug_mode:
        print_experimentsteps(experimentsteps)
        # Exit
        sys.exit(0)

    # Does the output directory exist? If not, make it
    try:
        os.makedirs(os.path.join(exp_dir, "output"))
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise

    # Initialize the EnTK App Manager
    amgr = AppManager(hostname=rcfg['rabbitmq']['hostname'], port=rcfg['rabbitmq']['port'], autoterminate=False)

    # Apply the resource configuration provided by the user
    amgr.resource_desc = facts.LoadResourceConfig(exp_dir)

    # Load the localization list
    if not os.path.isfile(os.path.join(exp_dir, "location.lst")):
        with open(os.path.join(exp_dir, "location.lst"), 'w') as templocationfile:
            templocationfile.write("New_York\t12\t40.70\t-74.01")
    amgr.shared_data = [os.path.join(exp_dir, "location.lst")]

    for pipelines in experimentsteps:
        # Assign the list of pipelines to the workflow
        amgr.workflow = pipelines

        # Run the SLR projection workflow
        amgr.run()

    # Close the application manager
    amgr.terminate()


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
    parser.add_argument('--debug', help="Enable debug mode", action="store_true")

    # Parse the arguments
    args = parser.parse_args()

    # Does the experiment directory exist?
    if not os.path.isdir(args.edir):
        print('%s does not exist'.format(args.edir))
        sys.exit(1)

    # Go ahead and try to run the experiment
    run_experiment(args.edir, args.debug)

    #sys.exit(0)
