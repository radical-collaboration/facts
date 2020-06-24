import os
import copy
from radical.entk import Pipeline, Stage, Task, AppManager

def run():
	
	# Initialize the EnTK App Manager
	amgr = AppManager(hostname="localhost", port=5672)
	
	# Apply the resource configuration provided by the user
	# Amarel
	res_desc = {'resource': "rutgers.amarel",
		'walltime': 30,
		'cpus': 16,
		'queue': "kopp_1",
		'project': "",
		'schema': "local"}
	# Local
	#res_desc = {'resource': "local.localhost",
	#	'walltime': 60,
	#	'cpus': 2,
	#	'queue': "",
	#	'project': "",
	#	'schema': ""}
	amgr.resource_desc = res_desc
	
	# Push the input data to the shared directory
	amgr.shared_data = ['CMIP6_CanESM5_Omon_piControl_r1i1p1f1_zos_6000-6199.nc', 'xarray_script.py']
	
	# New pipeline
	p1 = Pipeline()
	p1.name = "Test-pipeline1"
	p2 = Pipeline()
	p2.name = "Test-pipeline2"
	
	# First stage with two tasks
	s1 = Stage()
	s1.name = "Test-stage1"
	
	s2 = Stage()
	s2.name = "Test-stage2"
	
	t1 = Task()
	t1.name = "Test-task1"
	t1.cpu_reqs = {"processes": 1, "process_type": None, "threads_per_process": 1, "thread_type": None}
	t1.pre_exec = ["pip3 install --upgrade; pip3 install pandas zarr cftime toolz \"dask[complete]\" bottleneck xarray"]
	t1.executable = 'python3'
	t1.arguments = ['xarray_script.py']
	t1.copy_input_data = ["$SHARED/CMIP6_CanESM5_Omon_piControl_r1i1p1f1_zos_6000-6199.nc", "$SHARED/xarray_script.py"]
	t1.download_output_data = ["test_netcdf_file.nc > test_netcdf_file1.nc"]
	
	t2 = copy.deepcopy(t1)
	t2.name = "Test-task2"
	t2.download_output_data = ["test_netcdf_file.nc > test_netcdf_file2.nc"]
	
	t3 = copy.deepcopy(t1)
	t3.name = "Test-task3"
	t3.download_output_data = ["test_netcdf_file.nc > test_netcdf_file3.nc"]
	
	t4 = copy.deepcopy(t1)
	t4.name = "Test-task4"
	t4.download_output_data = ["test_netcdf_file.nc > test_netcdf_file4.nc"]
	
	# Assign tasks and stages to pipeline
	s1.add_tasks(t1)
	s1.add_tasks(t2)
	p1.add_stages(s1)
	
	s2.add_tasks(t3)
	s2.add_tasks(t4)
	p2.add_stages(s2)
	
	# Assign the pipeline to the workflow and run
	amgr.workflow = [p1, p2]
	amgr.run()
	
	# Done
	return(None)


if __name__ == "__main__":
	
	run()
	
	exit()