import numpy as np
import os
import argparse
import yaml
import sys
import re
import errno
from pprint import pprint
from radical.entk import Pipeline, Stage, Task, AppManager

# Magic variable replacement functions ================================
def mvar_replace_var(mvar, sub, string):
	mvar = "%{}%".format(mvar)
	return(re.sub(mvar, str(sub), string))

def mvar_replace_dict(dict, string):
	for var in dict.keys():
		string = mvar_replace_var(var, dict[var], string)
	return(string)

def mvar_replace_list(dict, list):
	return([mvar_replace_dict(dict, x) for x in list])

# =====================================================================

def GeneratePipeline(pcfg, ecfg, pipe_name, exp_dir):
	
	# Append the exp_dir to the ecfg dictionary to simplify things a bit
	ecfg['exp_dir'] = exp_dir
	
	# Append the input file to the list of options (if need be)
	if "input_data_file" in ecfg.keys():
		ecfg['options']['input_data_file'] = ecfg['input_data_file']
	
	# Append the pipeline id to the list of options
	ecfg['options']['pipeline_id'] = pipe_name
	
	# Initialize the pipeline object
	p = Pipeline()
	
	# Give the pipeline a name
	p.name = pipe_name
	
	# Loop through the necessary stages for this module
	stage_names = ["pre-process", "fit", "project", "post-process"]
	for this_stage in stage_names:
		if this_stage in pcfg.keys():
			
			# Populate the pipeline with the stages
			p.add_stages(GenerateStage(pcfg[this_stage], ecfg, p.name, this_stage))
	
	return(p)


def GenerateStage(scfg, ecfg, pipe_name, stage_name):
	
	# Initialize a stage object
	s = Stage()
	
	# Provide a name for this stage
	s.name=stage_name
	
	# Loop through the tasks for this stage
	for this_task in scfg.keys():
		
		# Populate the stage object with the tasks
		s.add_tasks(GenerateTask(scfg[this_task], ecfg, pipe_name, stage_name, this_task))
	
	# Return the stage object
	return(s)


def GenerateTask(tcfg, ecfg, pipe_name, stage_name, task_name):
	
	# Initialize a task object
	t = Task()
	
	# Define magic variable dictionary
	mvar_dict = {"PIPELINE_ID": pipe_name}
	
	# Give this task object a name
	t.name=task_name
	
	# Pre exec let you load modules, set environment before executing the workload
	if tcfg['pre_exec'] != "":
		t.pre_exec = [tcfg['pre_exec']]
	
	# Executable to use for the task
	t.executable = tcfg['executable']
	
	# If there's a user-defined input file (likely for genmod modules), add it to the
	# options list and upload file list if needed
	if "input_data_file" in tcfg['options']:
		tcfg['upload_input_data'].append(os.path.join(ecfg['exp_dir'], "input", ecfg['input_data_file']))

	# List of arguments for the executable
	t.arguments = [tcfg['script']] + match_options(tcfg['options'], ecfg['options'])

	# CPU requirements for this task
	t.cpu_threads = {
						'processes': tcfg['cpu']['processes'],
						'process-type': tcfg['cpu']['process-type'],
						'threads-per-process': tcfg['cpu']['threads-per-process'],
						'thread-type': tcfg['cpu']['thread-type'],
					}

	# Upload data from your local machine to the remote machine
	# Note: Remote machine can be the local machine
	t.upload_input_data = tcfg['upload_input_data']

	# Copy data from other stages/tasks for use in this task
	copy_list = []
	if "copy_input_data" in tcfg.keys():
		for copy_stage in tcfg['copy_input_data'].keys():
			for copy_task in tcfg['copy_input_data'][copy_stage].keys():
				loc = "$Pipeline_{0}_Stage_{1}_Task_{2}".format(pipe_name, copy_stage, copy_task)
				copy_list.extend(['%s/%s'%(loc, mvar_replace_dict(mvar_dict,x)) for x in tcfg['copy_input_data'][copy_stage][copy_task]])
	
	# Append the copy list (if any) to the task object	
	t.copy_input_data = copy_list
	
	# Set the download data for the task
	download_list = []
	outdir = os.path.join(ecfg['exp_dir'], "output")
	if "download_output_data" in tcfg.keys():
		download_list.extend(['%s > %s/%s'%(mvar_replace_dict(mvar_dict,x),outdir,mvar_replace_dict(mvar_dict,x)) for x in tcfg['download_output_data']]) 
	
	# Append the download list to this task
	t.download_output_data = download_list

	# Return the task object
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
		
	# Does the output directory exist? If not, make it
	try:
		os.makedirs(os.path.join(exp_dir, "output"))
	except OSError as e:
		if e.errno != errno.EEXIST:
			raise
	
	# Loop through the user-requested modules
	for this_mod in ecfg.keys():
		
		# Load the pipeline configuration file for this module
		pcfg_file = os.path.join(os.path.dirname(__file__), "modules", ecfg[this_mod]['module_set'], ecfg[this_mod]['module'], "pipeline.yml")
		if not os.path.isfile(pcfg_file):
			print '%s does not exist' % pcfg_file
			sys.exit(1)
		with open(pcfg_file, 'r') as fp:
			pcfg = yaml.safe_load(fp)
		
		# Generate a pipeline for this module
		pipe_name = "-".join((this_mod, ecfg[this_mod]['module_set'], ecfg[this_mod]['module']))
		pipelines.append(GeneratePipeline(pcfg, ecfg[this_mod], pipe_name, exp_dir))
	
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
	
	# Initialize the EnTK App Manager
	amgr = AppManager(hostname=rcfg['rabbitmq']['hostname'], port=rcfg['rabbitmq']['port'])
	
	# Apply the resource configuration provided by the user
	res_desc = {'resource': rcfg['resource-desc']['name'],
		'walltime': rcfg['resource-desc']['walltime'],
		'cpus': rcfg['resource-desc']['cpus'],
		'queue': rcfg['resource-desc']['queue'],
		'project': rcfg['resource-desc']['project']}
	amgr.resource_desc = res_desc
	
	# Assign the list of pipelines to the workflow
	amgr.workflow = pipelines
	
	# Run the workflow
	amgr.run()

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

	
	#sys.exit(0)