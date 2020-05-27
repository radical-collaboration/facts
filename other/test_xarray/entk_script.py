import os
from radical.entk import Pipeline, Stage, Task, AppManager

def run():
	
	# Initialize the EnTK App Manager
	amgr = AppManager(hostname="localhost", port=5672)
	
	# Apply the resource configuration provided by the user
	res_desc = {'resource': "local.localhost",
		'walltime': 30,
		'cpus': 1,
		'queue': "",
		'project': ""}
	amgr.resource_desc = res_desc
	
	# New pipeline
	p1 = Pipeline()
	p1.name = "Test-pipeline"
	
	# First stage with one task
	s1 = Stage()
	s1.name = "Test-stage"
	t1 = Task()
	t1.name = "Test-task"
	t1.pre_exec = ["pip3 install --upgrade; pip3 install pandas zarr cftime toolz \"dask[complete]\" bottleneck xarray"]
	t1.executable = 'python3'
	t1.arguments = ['xarray_script.py']
	t1.upload_input_data = ["CMIP6_CanESM5_Omon_piControl_r1i1p1f1_zos_6000-6199.nc", "xarray_script.py"]
	t1.download_output_data = ["test_netcdf_file.nc"]
	
	# Assign tasks and stages to pipeline
	s1.add_tasks(t1)
	p1.add_stages(s1)
	
	# Assign the pipeline to the workflow and run
	amgr.workflow = [p1]
	amgr.run()
	
	# Done
	return(None)


if __name__ == "__main__":
	
	run()
	
	exit()