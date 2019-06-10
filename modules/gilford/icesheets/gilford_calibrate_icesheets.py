import argparse

''' gilford_calibrate_icesheets.py

This script runs the icesheet calibration task for the Gilford icesheet emulator. 

Parameters: 
None

Output:
None

Notes: At this point, the emulator is being generated prior to its introduction into
the framework.  This stage will be necessary, but at a later date.  For now, it's
simply a NULL process.

'''

def gilford_calibrate_icesheets():
	
	return(None)


if __name__ == '__main__':
	
	# Initialize the command-line argument parser
	parser = argparse.ArgumentParser(description="Run the ice sheet calibration stage for the Gilford Ice Sheet Emulator workflow",\
	epilog="Note: This is meant to be run as part of the Gilford Ice Sheet module within the Framework for the Assessment of Changes To Sea-level (FACTS)")
	
	# Define the command line arguments to be expected
	
	# Parse the arguments
	args = parser.parse_args()
	
	# Run the calibration stage
	gilford_calibrate_icesheets()
	
	exit()