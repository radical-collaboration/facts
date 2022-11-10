
import os
import yaml
import sys
import re
from radical.entk import Pipeline, Stage, Task


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
def GeneratePipeline(pcfg, ecfg, pipe_name, exp_dir, stage_names=None, workflow_name="", scale_name=""):

    if not stage_names:
        stage_names = ["preprocess", "fit", "project", "postprocess"]

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
    if 'stages' in ecfg.keys():
        stage_names = ecfg['stages']

    for this_stage in stage_names:
        if this_stage in pcfg.keys():

            # Populate the pipeline with the stages
            p.add_stages(GenerateStage(pcfg[this_stage], ecfg, p.name, this_stage, workflow_name=workflow_name, scale_name=scale_name))

    return(p)


def GenerateStage(scfg, ecfg, pipe_name, stage_name, workflow_name="", scale_name=""):

    # Initialize a stage object
    s = Stage()

    # Provide a name for this stage
    s.name = stage_name

    # Loop through the tasks for this stage
    for this_task in scfg.keys():

        # Populate the stage object with the tasks
        s.add_tasks(GenerateTask(scfg[this_task], ecfg, pipe_name, stage_name, this_task, workflow_name=workflow_name, scale_name=scale_name))

    # Return the stage object
    return(s)


def GenerateTask(tcfg, ecfg, pipe_name, stage_name, task_name, workflow_name="", scale_name=""):

    # Initialize a task object
    t = Task()

    # Define magic variable dictionary
    module_path = os.path.join('.', 'modules', ecfg['module_set'],ecfg['module'] )
    mvar_dict = {"PIPELINE_ID": pipe_name, "WORKFLOW_NAME": workflow_name, "SCALE_NAME": scale_name,
        "MODULE_SET_NAME": ecfg['module_set'], "MODULE_NAME": ecfg['module'],
        "MODULE_PATH": module_path,
        "CLIMATE_DATA_FILE": ecfg['options']['climate_data_file'], "CLIMATE_GSAT_FILE": ecfg['options']['climate_gsat_data_file'],
        "CLIMATE_OHC_FILE": ecfg['options']['climate_ohc_data_file'],
        "EXP_DIR": ecfg['exp_dir'], "EXPERIMENT_NAME": ecfg['options']['experiment_name']}
    # Give this task object a name
    t.name = task_name

    # If there's a user-defined input file (likely for genmod modules), add it to the
    # options list and upload file list if needed
    # if "input_data_file" in tcfg['options']:
    #     tcfg['upload_input_data'].extend(os.path.join(ecfg['exp_dir'], "input", ecfg['input_data_file']))

    if "input_data_file" in ecfg.keys():
        for x in ecfg['input_data_file']:
            fp = os.path.join(ecfg['exp_dir'], "input", x)
            if not os.path.isfile(fp):
                raise(FileNotFoundError("input_data_file: " + fp + " not found!"))
            tcfg['upload_input_data'].append(fp)

    # Upload data from your local machine to the remote machine
    # Note: Remote machine can be the local machine
    t.upload_input_data = []
    for this_file0 in tcfg['upload_input_data']:
        this_file = mvar_replace_dict(mvar_dict,this_file0)
        if not os.path.isfile(this_file):
            # inelegant, but don't raise an exception if we are uploading workflows.yml, which will be
            # created later
            if not os.path.basename(this_file) == "workflows.yml":
                raise(FileNotFoundError(pipe_name + "." + stage_name + ": upload_input_data: " + this_file + " not found!"))
        t.upload_input_data.append(this_file)

    # Copy data from other stages/tasks for use in this task
    copy_list = []
    if "copy_input_data" in tcfg.keys():
        for copy_stage in tcfg['copy_input_data'].keys():
            for copy_task in tcfg['copy_input_data'][copy_stage].keys():
                loc = "$Pipeline_{0}_Stage_{1}_Task_{2}".format(pipe_name, copy_stage, copy_task)
                copy_list.extend(['{0}/{1}'.format(loc, mvar_replace_dict(mvar_dict, x))
                                 for x in tcfg['copy_input_data'][copy_stage][copy_task]])

    # Copy data from the shared directory
    if "copy_shared_data" in tcfg.keys():
        copy_list.extend(['{0}'.format(mvar_replace_dict(mvar_dict, x))
                         for x in tcfg['copy_shared_data']])


    # Append the copy list (if any) to the task object
    t.copy_input_data = copy_list

    # Pre exec let you load modules, set environment before executing the workload
    t.pre_exec = []

    # If there's a data file to upload and extract, add it to upload and
    # add the extraction command to pre-exec
    if "upload_and_extract_input_data" in tcfg.keys():
        for this_file0 in tcfg['upload_and_extract_input_data']:
            this_file = mvar_replace_dict(mvar_dict,this_file0)
            if not os.path.isfile(this_file):
                raise(FileNotFoundError(pipe_name + "." + stage_name + ": upload_and_extract_input_data: " + this_file + " not found!"))
            t.pre_exec.append('tar -xvf ' + os.path.basename(this_file) + ' 2> /dev/null; rm ' + os.path.basename(this_file))
            t.upload_input_data.append(this_file)

   # if their are python dependencies, add pip call to pre_exec
    if "python_dependencies" in tcfg.keys():
        t.pre_exec.append('pip install --upgrade pip; pip install ' + tcfg['python_dependencies'])

    if "pre_exec" in tcfg.keys():
        if len(tcfg['pre_exec']) > 0:
            t.pre_exec.append(tcfg['pre_exec'])

    # Executable to use for the task
    t.executable = tcfg['executable']




    # List of arguments for the executable
    # t.arguments = [tcfg['script']] + match_options(tcfg['options'], ecfg['options'])
    t.arguments = [tcfg['script']]
    if "arguments" in tcfg.keys():
        t.arguments += [mvar_replace_dict(mvar_dict,x)  for x in tcfg['arguments']]
    for x in match_options(tcfg['options'], ecfg['options']):
        if type(x)==str:
            t.arguments.append(mvar_replace_dict(mvar_dict,x))
        else:
            t.arguments.append(x)

    # CPU requirements for this task
    t.cpu_reqs = {
                     'cpu_processes': tcfg['cpu']['processes'],
                     'cpu_process_type': tcfg['cpu']['process-type'],
                     'cpu_threads': tcfg['cpu']['threads-per-process'],
                     'cpu_thread_type': tcfg['cpu']['thread-type'],
                 }


    # Send the global and local files to the shared directory for totaling
    copy_output_list = []
    if "copy_output_data" in tcfg.keys():
        copy_output_list.extend(['{0} > $SHARED/{0}'.format(mvar_replace_dict(mvar_dict, x))
                                for x in tcfg['copy_output_data']])

    if "climate_output_data" in tcfg.keys():
        copy_output_list.extend(['{0} > $SHARED/climate/{0}'.format(mvar_replace_dict(mvar_dict, x))
                                for x in tcfg['climate_output_data']])

    if "global_total_files" in tcfg.keys():
        copy_output_list.extend(['{0} > $SHARED/to_total/global/{0}'.format(mvar_replace_dict(mvar_dict, x))
                                for x in tcfg['global_total_files']])

    if "local_total_files" in tcfg.keys():
        copy_output_list.extend(['{0} > $SHARED/to_total/local/{0}'.format(mvar_replace_dict(mvar_dict, x))
                                for x in tcfg['local_total_files']])

    if "totaled_files" in tcfg.keys():
        copy_output_list.extend(['{0} > $SHARED/totaled/{0}'.format(mvar_replace_dict(mvar_dict, x))
                                for x in tcfg['totaled_files']])

    # Set the download data for the task
    download_list = []
    outdir = os.path.join(ecfg['exp_dir'], "output")
    if "download_output_data" in tcfg.keys():
        download_list.extend(['{0} > {1}/{0}'.format(mvar_replace_dict(mvar_dict, x), outdir)
                             for x in tcfg['download_output_data']])

    # Append the "total" lists to the copy output list
    t.copy_output_data = copy_output_list

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

def ParsePipelineConfig(this_mod, modcfg, global_options={}, relabel_mod=''):

    # Load the pipeline configuration file for this module

    if not "pipeline_file" in modcfg.keys():
        modcfg['pipeline_file'] = "pipeline.yml" 

    if 'module_set' in modcfg.keys():
        pcfg_file = os.path.join(os.path.dirname(__file__), "modules",
                       modcfg['module_set'], modcfg['module'], modcfg['pipeline_file'])
    else:
        pcfg_file = os.path.join(os.path.dirname(__file__), "modules", modcfg['module'], modcfg['pipeline_file'])

    if not os.path.isfile(pcfg_file):
        raise(FileNotFoundError(pcfg_file + " does not exist"))

    with open(pcfg_file, 'r') as fp:
        pcfg = yaml.safe_load(fp)

    if "options" not in modcfg.keys():
        modcfg["options"] = {}
    if "options_allowoverwrite" not in modcfg.keys():
        modcfg["options_allowoverwrite"] = {}
    modcfg["options_allowoverwrite"].update(global_options)

    # Append the global options to this module
    for this_opt in modcfg["options_allowoverwrite"]:
        if not this_opt in modcfg["options"].keys():
            if this_opt in global_options.keys():
                modcfg["options"][this_opt] = global_options[this_opt]

    if len(relabel_mod) == 0:
        relabel_mod = this_mod

    # Generate a pipeline for this module
    if 'module_set' in modcfg.keys():
        pipe_name = ".".join((relabel_mod, modcfg['module_set'], modcfg['module']))
    else:
        pipe_name = ".".join((relabel_mod, modcfg['module']))
    
    if "experiment_name" in global_options.keys():
        pipe_name = ".".join((global_options['experiment_name'],pipe_name))

    p = {
        "modlabel": this_mod,
        "pcfg": pcfg,
        "modcfg": modcfg,
        "pipe_name": pipe_name,
        "module_set": modcfg['module_set'],
        "module": modcfg['module']
    }
    return p


def IdentifyOutputFiles(pcfg,pipe_name):
    p={'global': [], 'local': []}

    # Define magic variable dictionary
    mvar_dict = {"PIPELINE_ID": pipe_name}

    for this_stage in pcfg:
        for this_task in pcfg[this_stage]:
            tcfg = pcfg[this_stage][this_task]
            if "global_total_files" in tcfg.keys():
                p['global'].extend(['global/{0}'.format(mvar_replace_dict(mvar_dict,x))
                                   for x in tcfg['global_total_files']])
            if "local_total_files" in tcfg.keys():
                p['local'].extend(['local/{0}'.format(mvar_replace_dict(mvar_dict,x))
                                  for x in tcfg['local_total_files']])
    return p

def IdentifyClimateOutputFiles(pcfg,pipe_name):
    pd={'climate': [], 'gsat': [], 'ohc': []}

    # Define magic variable dictionary
    mvar_dict = {"PIPELINE_ID": pipe_name}

    for this_stage in pcfg:
        for this_task in pcfg[this_stage]:
            tcfg = pcfg[this_stage][this_task]
            if "climate_output_data" in tcfg.keys():
                for this_file in tcfg['climate_output_data']:
                    if this_file.__contains__('climate.nc'):
                        pd['climate'] = '$SHARED/climate/' + mvar_replace_dict(mvar_dict, this_file)
                    elif this_file.__contains__('gsat.nc'):
                        pd['gsat'] = '$SHARED/climate/' + mvar_replace_dict(mvar_dict, this_file)
                    elif this_file.__contains__('ohc.nc'):
                        pd['ohc'] = '$SHARED/climate/' + mvar_replace_dict(mvar_dict, this_file)
 
    return pd

def ParseExperimentConfig(exp_dir):
    # Initialize a list for experiment steps (each step being a set of pipelines)
    experimentsteps = {}

    # Define the configuration file names
    cfile = os.path.join(exp_dir, "config.yml")

    # Does the experiment configuration file exist?
    if not os.path.isfile(cfile):
        raise(FileNotFoundError(cfile + ' does not exist'))

    # Load the resource and experiment configuration files
    with open(cfile, 'r') as fp:
        ecfg = yaml.safe_load(fp)

    # Reserved configuration entries
    reserved_econfig_entries = ["global-options", "total-options", "extremesealevel-options", "multistep", "include_in_workflow"]

    # Are there global options?
    if "global-options" in ecfg.keys():
        global_options = ecfg["global-options"]
    else:
        global_options = {}
        ecfg["global-options"] = global_options

    # set up global options for climate data files
    climate_data_files = []
    global_options['climate_data_file'] = None
    global_options['climate_gsat_data_file'] = None
    global_options['climate_ohc_data_file'] = None
    
    # add experiment name to global options
    global_options['experiment_name'] = os.path.basename(os.path.dirname(exp_dir))

    # Initialize a list for pipelines
    pipelines = []
    workflows_to_include = {}

    # Loop through the user-requested modules
    for this_mod in ecfg.keys():

        # Skip this entry if it's not associated with SLR projection workflow
        if this_mod in reserved_econfig_entries:
            continue

        for this_mod_sub in ecfg[this_mod].keys():
            if (this_mod_sub in reserved_econfig_entries):
                continue

            parsed = ParsePipelineConfig(this_mod_sub, ecfg[this_mod][this_mod_sub], global_options=global_options)
            #print(parsed['pipe_name'])

            # loop over workflows/scales if requested
            if "loop_over_workflows" in ecfg[this_mod][this_mod_sub].keys():
                for this_workflow in workflows_to_include:
                    if "loop_over_scales" in ecfg[this_mod][this_mod_sub].keys():
                        for this_scale in workflows_to_include[this_workflow]:
                            if len(workflows_to_include[this_workflow][this_scale]) > 0:
                                pipelines.append(GeneratePipeline(parsed['pcfg'], parsed['modcfg'], parsed['pipe_name'] + "." + this_workflow + "." + this_scale, exp_dir, workflow_name=this_workflow, scale_name=this_scale))
                    else:
                        pipelines.append(GeneratePipeline(parsed['pcfg'], parsed['modcfg'], parsed['pipe_name'] + "." + this_workflow , exp_dir, workflow_name=this_workflow))
            # specify workflow/scale if requested
            elif "workflow" in ecfg[this_mod][this_mod_sub].keys():
                this_workflow = ecfg[this_mod][this_mod_sub]['workflow']
                if "scale" in ecfg[this_mod][this_mod_sub].keys():
                    this_scale = ecfg[this_mod][this_mod_sub]['scale']
                    pipelines.append(GeneratePipeline(parsed['pcfg'], parsed['modcfg'], parsed['pipe_name'] + "." + this_workflow + "." + this_scale, exp_dir, workflow_name=this_workflow, scale_name=this_scale))
                else:
                    pipelines.append(GeneratePipeline(parsed['pcfg'], parsed['modcfg'], parsed['pipe_name'] + "." + this_workflow, exp_dir, workflow_name=this_workflow))
            else:
                pipelines.append(GeneratePipeline(parsed['pcfg'], parsed['modcfg'], parsed['pipe_name'], exp_dir))

            if "include_in_workflow" in ecfg[this_mod][this_mod_sub].keys():
                outfiles = IdentifyOutputFiles(parsed['pcfg'], parsed['pipe_name'])
                for this_wf in ecfg[this_mod][this_mod_sub]['include_in_workflow']:
                    if not this_wf in workflows_to_include.keys():
                        workflows_to_include[this_wf] = {'global': [], 'local': []}
                    for this_scale in outfiles:
                        workflows_to_include[this_wf][this_scale].extend(outfiles[this_scale])

            if "generates_climate_output" in ecfg[this_mod][this_mod_sub].keys():
                climate_data_files = IdentifyClimateOutputFiles(parsed['pcfg'], parsed['pipe_name'])
                global_options['climate_data_file'] = climate_data_files['climate']
                global_options['climate_gsat_data_file'] = climate_data_files['gsat']
                global_options['climate_ohc_data_file'] = climate_data_files['ohc']


        experimentsteps[this_mod] = pipelines
        pipelines = []

    return {'experimentsteps': experimentsteps, 'ecfg': ecfg, 'workflows': workflows_to_include, 'climate_data_files': climate_data_files}


def LoadResourceConfig(exp_dir, rcfg_name):

    if rcfg_name:
        rcfg_fname = 'resource_%s.yml' % rcfg_name
    else:
        rcfg_fname = 'resource.yml'

    # Define the configuration and resource file names
    rfile = os.path.join(exp_dir, rcfg_fname)

    # Does the resource file exist?
    if not os.path.isfile(rfile):
        raise(FileNotFoundError(rfile + '{} does not exist'))

    # Load the resource and experiment configuration files
    with open(rfile, 'r') as fp:
        rcfg = yaml.safe_load(fp)

    return rcfg


if __name__ == "__main__":

    print("")
    print("FACTS.py is intended to be called as a library, not from the command line.")
    print("See runFACTS.py for an example.")
    sys.exit(1)

    # Initialize the argument parser
    # parser = argparse.ArgumentParser(description="The Framework for Assessing Changes To Sea-level (FACTS)")

    # # Add arguments for the resource and experiment configuration files
    # parser.add_argument('edir', help="Experiment Directory")
    # parser.add_argument('--no-total', help="Disable totaling of global and local sea-level projection wokflow", action="store_true")
    # parser.add_argument('--debug', help="Enable debug mode", action="store_true")

    # # Parse the arguments
    # args = parser.parse_args()

    # # Does the experiment directory exist?
    # if not os.path.isdir(args.edir):
    #     print('%s does not exist'.format(args.edir))
    #     sys.exit(1)

    # Go ahead and try to run the experiment
    # run_experiment(args.edir, args.debug, args.no_total)

    #sys.exit(0)
