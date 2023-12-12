.. _chapter_configuration:

Experiment configuration
========================

Each experiment directory includes a config.yml file, which defines the configuration of the experiment.


Global Options
--------------

The **global-options** section of config.yml specifies options that should, by default, apply to every module (but can be overwritten by module-level configuration). Common global options include:

* **nsamps**: number of samples to run
* **scenario**: climate scenario to use
* **pyear_start**: year in which to start projection
* **pyear_step**: step size between projection years
* **pyear_end**: year in which to end projections
* **baseyear**: zero point for sea level projections
* **pipeline_file**: default pipeline file name (optional; defaults to 'pipeline.yml'). See *pipeline_file* under experiment step configuration.

In addition, this section can specify:

* **rcfg-name**: the name of the resource configuration file to be used. FACTS will look in the resources/ directory for file with name *resource_(rcfg_name).yml*. (If not specified, *resource.yml* will be used.)

Global options can be overwritten when runFACTS.py is invokved from the command line by passing a dictionary along with the ``--global_options`` parameter; for example::

    python3 runFACTS.py --global_options '{"rcfg-name": "localhost" }' experiments/dummy/


Experiment step and module configuration
----------------------------------------

All other top-level entries in config.yml specify experiment steps. The label of these are arbitrary, though common labels are *climate_step*, *sealevel_step*, *totaling_step* and *esl_step*.

The second-level entries under the experiment steps specify the modules to be run. Module labels are arbitrary, and are used together with the experiment name to label output files.

The following third-level entries are used under the module label:

* **module_set** (required): the name of the module set (first level hierarchy in the modules/ directory)

* **module** (required): the name of the module (second level hierarchy in the climate directory)

* **options**: Entries under this header are passed as options to the modules if they match options defined in the module pipeline.yml. These overwrite global options.

* **options_allowoverwrite**: Entries under this header are passed as options to the modules if they match options defined in the module pipeline.yml. These can be overwriten by global options. Options that make use of magic variables should be specified here. A particularly useful example is the specification ``climate_data_file: "%CLIMATE_DATA_FILE%"``, which allows the identified of the climate data file generate in the climate stage to be passed onto sea level modules that use it as input.

* **input_data_file**: Identifies a data file to be uploaded into the sandbox from the input/ subdirectory of the experiment directory.

* **input_compressed_data_file**: Identifies a .tgz file to be uploaded into the sandbox from the input/ subdirectory of the experiment directory and extracted.

* **stages**: Specifies stages from the pipeline.yml file to be run. Defaults to 'preprocess', 'fit', 'project', 'postprocess'.

* **pipeline_file**: Pipeline file name. Defaults to 'pipeline.yml'. Alternatives can be useful for special cases (e.g., using only Antarctic ice sheet output from a module that produces both Greenland and Antarctic ice sheet output.) Many modules include a 'pipeline.global.yml' file for running the module without producing localized output.

* **include_in_workflow**: A list of all workflows the module output should be included in at the totaling steps.

* **loop_over_workflows**: If defined, replicate the module for all workflows defined to date.

* **loop_over_scales**: If defined, replicate the module for both global and local scale (e.g., for a totaling module).

* **generates_climate_output**: Module generates climate output data. Globally set ```climate_output_data``` based on the module's pipeline configuration.

* **climate_output_data**, **global_total_files**, **local_total_files**, **totaled_files**: See definitions in pipeline configuration file specification. Can be used together with the `facts/dummy` module and `input_data_file` to upload outputs produced by prior Experiment Steps and have them filed in an appropriate directory for subsequent modules to find. 


Example experiment file
^^^^^^^^^^^^^^^^^^^^^^^

::

    global-options:
        nsamps: 2000
        scenario: ssp585
        pyear_start: 2020
        pyear_end: 2150
        pyear_step: 10
        baseyear: 2005

    climate_step:
        temperature:
            module_set: "fair"
            module: "temperature"
            generates_climate_output: true
    #        input_data_file:
    #              - "emissions.csv"
    #        options:
    #              rcmip_file: emissions.csv

    sealevel_step:
        GrIS1f:
            module_set: "FittedISMIP"
            module: "GrIS"
            include_in_workflow:
                - "wf1f"
                - "wf2f"
                - "wf3f"

        deconto21:
            module_set: "deconto21"
            module: "AIS"
            include_in_workflow:
                - "wf3e"   
                - "wf3f"

        bamber19:
            module_set: "bamber19"
            module: "icesheets"
            include_in_workflow:
                - "wf4" 

        emuAIS:
            module_set: "emulandice"
            module: "AIS"
            include_in_workflow:
                - "wf1e"      
            options:
                pyear_end: 2100

        emuGrIS:
            module_set: "emulandice"
            module: "GrIS"
            include_in_workflow:
                - "wf1e"  
                - "wf2e"  
                - "wf3e"  
            options:
                pyear_end: 2100

        emuglaciers:
            module_set: "emulandice"
            module: "glaciers"
            include_in_workflow:
                - "wf1e"   
                - "wf2e"   
                - "wf3e"
            options:
                pyear_end: 2100

        larmip:
            module_set: "larmip"
            module: "AIS"
            include_in_workflow:
                - "wf2e"
                - "wf2f"

        ar5glaciers:
            module_set: "ipccar5"
            module: "glaciers"
            options:
                gmip: 2
            include_in_workflow:
                - "wf1f"
                - "wf2f"
                - "wf3f"
                - "wf4"

        ar5AIS:
            module_set: "ipccar5"
            module: "icesheets"
            pipeline_file: "pipeline.AIS.yml"
            include_in_workflow:
                - "wf1f"

        ocean:
            module_set: "tlm"
            module: "sterodynamics"
            include_in_workflow:
                - "wf1f"
                - "wf1e"
                - "wf2e"
                - "wf2f"
                - "wf3e"
                - "wf3f"
                - "wf4"

        k14vlm:
            module_set: "kopp14"
            module: "verticallandmotion"
            include_in_workflow:
                - "wf1f"
                - "wf1e"
                - "wf2e"
                - "wf2f"
                - "wf3e"
                - "wf3f"
                - "wf4"

        lws:
            module_set: "ssp"
            module: "landwaterstorage"
            options:
                scenario: "ssp5"
                dcrate_lo: -0.4
            include_in_workflow:
                - "wf1f"
                - "wf1e"
                - "wf2e"
                - "wf2f"
                - "wf3e"
                - "wf3f"
                - "wf4"

    totaling_step:
        total:
            module_set: "facts"
            module: "total"
            loop_over_workflows: true
            loop_over_scales: true
            stages:
                - workflow

    esl_step:
        extremesealevel:
            loop_over_workflows: true
            module_set: "extremesealevel"
            module: "pointsoverthreshold"
            options:
                target_years: 2050,2100
                total_localsl_file: "$SHARED/totaled/%EXPERIMENT_NAME%.total.workflow.%WORKFLOW_NAME%.local.nc" 



Pipeline configuration
======================

The top level entries in the pipeline configuration identify the stages of the module. Most commonly, these are 'preprocess', 'fit', 'project' and 'postprocess', but any labels are permitted. If other than these four, the stages must be specifically identified in the config.yml file.

The next level entry should be "task1:"

The third level entry defines the task, including the executable to be run, parameters to be passed to it, and associated files. Entries include:

* **executable** (required): The name of the executable to run (e.g., 'python3').
* **upload_input_data**: Files to be uploaded prior to run. If path is not specified, looks in the module directory.
* **upload_and_extract_input_data**: .tgz files to be uploaded and extracted prior to run. If path is not specified, looks in the modules-data directory.
* **cpu**: The start of an entry defining computational requirements. Fields are:
    * **cpu_processes**
    * **cpu_process_type**
    * **cpu_threads**
    * **cpu_thread_type**
* **python_dependencies**: Python modules to be installed via pip prior to execution, contained as a single string separated by spaces. (Only needs to be included for the first stage in a module, as these will persist to subsequent changes.)
* **script**: The name of the script file to run with executable. Will be uploaded from the module directory.
* **script_noupload**: The name of the script file to run with executable, but do not look in the module directory for the file to upload. (Useful if script is in a different location.)
* **arguments**: Arguments to be passed to the script file.
* **options**: A list of option names to be passed to the script file if their value is specified in config.yml.
* **pre_exec**: Any commands, not otherwise specified, to be run before execution.
* **copy_input_data**: A hierarchy (next level being stage name, subsequent level being 'task1', third level being file names) of files to copy from previous stages.
* **copy_shared_data**: A list of data to be copied from a shared storage area for this experiment.
* **climate_output_data**: A list of climate output data to be coped to a shared (cross-module) directory for use in subsequent stages and downloaded.
* **global_total_files**: A list of global output files to be copied to a shared (cross-module) directory for use in the totaling stage and downloaded.
* **local_total_files**: A list of local output files to be copied to a shared (cross-module) directory for use in the totaling stage and downloaded.
* **totaled_files**: A list of total sea level files to be copied to a shared (cross-module) directory for use in post-totaling stage and downloaded.
* **copy_output_data**: A listing of output files, not otherwise specified, to be copied to a shared (cross-module) directory for subsequent use.
* **download_output_data**: A listing of additional output files to be downloaded.

Example pipeline.yml file
-------------------------

::

    preprocess:
    task1:
        executable: "python3"
        python_dependencies: "numpy scipy netCDF4 pyyaml matplotlib"
        script: "bamber19_preprocess_icesheets.py"
        options:
        - "pipeline_id"
        - "pyear_start"
        - "pyear_end"
        - "pyear_step"
        - "scenario"
        - "baseyear"
        upload_and_extract_input_data:
        - "bamber19_icesheets_preprocess_data.tgz"


    fit:
    task1:
        executable: "python3"
        script: "bamber19_fit_icesheets.py"
        options:
        - "pipeline_id"

    project:
    task1:
        executable: "python3"
        script: "bamber19_project_icesheets.py"
        options:
        - "nsamps"
        - "seed"
        - "replace"
        - "pipeline_id"
        copy_input_data:
        preprocess:
            task1:
            - "%PIPELINE_ID%_data.pkl"
        global_total_files:
        - "%PIPELINE_ID%_GIS_globalsl.nc"
        - "%PIPELINE_ID%_AIS_globalsl.nc"
        download_output_data:
        - "%PIPELINE_ID%_EAIS_globalsl.nc"
        - "%PIPELINE_ID%_WAIS_globalsl.nc"

    postprocess:
    task1:
        executable: "python3"
        script: "bamber19_postprocess_icesheets.py"
        options:
        - "locationfile"
        - "pipeline_id"
        upload_input_data:
        - 'read_locationfile.py'
        - 'AssignFP.py'
        - 'ReadFingerprint.py'
        upload_and_extract_input_data:
        - 'bamber19_icesheets_postprocess_data.tgz'
        copy_shared_data:
        - '$SHARED/location.lst'
        copy_input_data:
        project:
            task1:
            - "%PIPELINE_ID%_projections.pkl"
        local_total_files:
        - "%PIPELINE_ID%_GIS_localsl.nc"
        - "%PIPELINE_ID%_AIS_localsl.nc"
        download_output_data:
        - "%PIPELINE_ID%_WAIS_localsl.nc"
        - "%PIPELINE_ID%_EAIS_localsl.nc".


Magic variables
===============

The following magic variables can be used in options specified in the config.yml and pipeline.yml files.

+-------------------------+---------------------------------------------------------------------------------------------+
| Variable                | Definition                                                                                  |
+=========================|=============================================================================================|
| ``%PIPELINE_ID%``       | Pipeline ID                                                                                 |
+-------------------------+---------------------------------------------------------------------------------------------+
| ``%WORKFLOW_NAME%``     | Workflow Name                                                                               |
+-------------------------+---------------------------------------------------------------------------------------------+
| ``%SCALE_NAME%``        | Name of scale (global, local)                                                               |
+-------------------------+---------------------------------------------------------------------------------------------+
| ``%MODULE_SET_NAME%``   | Module set name                                                                             |
+-------------------------+---------------------------------------------------------------------------------------------+
| ``%MODULE_NAME%``       | Module name                                                                                 |
+-------------------------+---------------------------------------------------------------------------------------------+
| ``%MODULE_PATH%``       | Module path                                                                                 |
+-------------------------+---------------------------------------------------------------------------------------------+
| ``%CLIMATE_DATA_FILE%`` | Climate data file produced by climate step, to be used as input (includes both GSAT and OHC)|
+-------------------------+---------------------------------------------------------------------------------------------+
| ``%CLIMATE_GSAT_FILE%`` | GSAT data file produced by climate step                                                     |
+-------------------------+---------------------------------------------------------------------------------------------+
| ``%CLIMATE_OHC_FILE%``  | OHC data file produced by climate step                                                      |
+-------------------------+---------------------------------------------------------------------------------------------+
| ``%EXP_DIR%``           | Experiment path                                                                             |
+-------------------------+---------------------------------------------------------------------------------------------+
| ``%EXPERIMENT_NAME%``   | Experiment name                                                                             |
+-------------------------+---------------------------------------------------------------------------------------------+

In addition, the variable $SHARED, used in a file path, directs scripts to a shared storage area in the sandbox.
