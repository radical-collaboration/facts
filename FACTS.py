import numpy as np
import os
import argparse
import yaml
import sys
from pprint import pprint
from radical.entk import Pipeline, Stage, Task, AppManager


def GeneratePipeline(wcfg, pcfg):
	
	# Initialize the pipeline object
	p = Pipeline()
	
	# Give the pipeline a name
	p.name = "_".join((pcfg['module_set'], pcfg['module']))
	
	# Loop through the necessary stages for this module
	for this_stage in wcfg.keys():
		print(this_stage)
	p.add_stages(GenerateStage(wcfg, pcfg, p.name))
	
	return(p)


def GenerateStage(wcfg, scfg, pipe_name):
	
	s = Stage()
	s.name="Test Stage"
	s.add_tasks(GenerateTask(wcfg, scfg))
	
	return(s)


def GenerateTask(wcfg, tcfg):
	
	t = Task()
	t.name="Test Task"
	#t = []
	
	return(t)



def match_options(wopts, eopts):
	
	# Initialize return list
	opt_list = []
	
	# Find the experiment's options that match this particular workflow
	matched_option_names = [i for i in wopts if i in eopts.keys()]
	
	# Generate the appropriate flag and append the matched options to the output list
	for matched_opt in matched_option_names:
		opt_list.append("--{0}".format(matched_opt))
		opt_list.append(eopts[matched_opt])
	
	# Return the matched options list
	return(opt_list)



def run_experiment(exp_dir, debug_mode):
	
	# Initialize a list of pipelines
	pipelines = []
	
	# Define the configuration and resource file names
	rfile = os.path.join(exp_dir, "resource.yml")
	cfile = os.path.join(exp_dir, "config.yml")
	
	# Does the experiment configuration file exist?
	if not os.path.isfile(cfile):
		print '%s does not exist' % cfile
		sys.exit(1)

	# Does the resource file exist?
	if not os.path.isfile(rfile):
	    print '%s does not exist' % rfile
	    sys.exit(1)
	
	# Load the resource and experiment configuration files
	with open(rfile, 'r') as fp:
		rcfg = yaml.safe_load(fp)
	with open(cfile, 'r') as fp:
		ecfg = yaml.safe_load(fp)
	
	# Loop through the user-requested modules
	for this_mod in ecfg.keys():
		
		# Load the workflow configuration file for this module
		wcfg_file = os.path.join(os.path.dirname(__file__), "modules", ecfg[this_mod]['module_set'], ecfg[this_mod]['module'], "workflow.yml")
		if not os.path.isfile(wcfg_file):
			print '%s does not exist' % wcfg_file
			sys.exit(1)
		with open(wcfg_file, 'r') as fp:
			wcfg = yaml.safe_load(fp)
		
		# Generate a pipeline for this module
		pipelines.append(GeneratePipeline(wcfg, ecfg[this_mod]))
	
	# Print out PST info if in debug mode
	if debug_mode:
		for p in pipelines:
			print("Pipeline {}:".format(p.name))
			print("################################")
			pprint(p.to_dict())
			for s in p.stages:
				print("Stage {}:".format(s.name))
				print("============================")
				pprint(s.to_dict())
				for t in s.tasks:
					print("Task {}:".format(t.name))
					print("----------------------------")
					pprint(t.to_dict())
				
		
		# Exit
		sys.exit(0)
	'''
    amgr = AppManager(hostname=rcfg['rabbitmq']['hostname'], port=rcfg['rabbitmq']['port'])
	
	res_desc = {'resource': rcfg['resource-desc']['name'],
		'walltime': rcfg['resource-desc']['walltime'],
		'cpus': rcfg['resource-desc']['cpus'],
		'queue': rcfg['resource-desc']['queue'],
		'project': rcfg['resource-desc']['project']}
	
	amgr.resource_desc = res_desc
	amgr.workflow = pipelines
	amgr.run()
	'''
	return(None)


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
		print '%s does not exist' % args.edir
		sys.exit(1)
	
	# Go ahead and try to run the experiment
	run_experiment(args.edir, args.debug)

	'''	
	for this_stage in wcfg.keys():
		print(this_stage)
		stage_options = wcfg[this_stage]['options']
		user_options = ucfg['icesheets']['options']
		matched_options = match_options(stage_options, user_options)
		print("{} {} {}".format(wcfg[this_stage]['executable'], wcfg[this_stage]['script'], " ".join(str(e) for e in matched_options)))
	'''
	
	sys.exit(0)