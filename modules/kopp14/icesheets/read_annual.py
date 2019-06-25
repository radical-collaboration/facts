import numpy as np
import os


''' read_annual.py

Reads in a directory of PSMSL data akin to the "read_annual.m" code provided
by PSMSL.

Parameters:
dataDir = Unzipped directory of PSMSL data
metaonly = Extract only the meta data for all the sites in dataDir

Return:
data = Object containing the PSMSL information


'''

def read_annual(dataDir, metaonly=False):
	
	class Site:
		def __init__(self, id, lat, lon, name, coastline, stationcode, stationflag, year=None, height=None, interpolated=None, dataflag=None, isMtl=None):
			self.name = name.strip()
			self.lat = float(lat)
			self.lon = float(lon)
			self.id = int(id)
			self.coastline = int(coastline)
			self.stationcode = int(stationcode)
			self.stationflag = stationflag.strip()
			if year is None:
				year = []
			self.year = year
			if height is None:
				height = []
			self.height = height
			if interpolated is None:
				interpolated=[]
			self.interpolated = interpolated
			if dataflag is None:
				dataflag = []
			self.dataflag = dataflag
			if isMtl is None:
				isMtl = []
			self.isMtl = isMtl
			
		def __repr__(self):
			out = "Name: {0}\n".format(self.name)
			out = out + "ID: {0}\n".format(self.id)
			out = out + "Latitude: {0}\n".format(self.lat)
			out = out + "Longitude: {0}\n".format(self.lon)
			out = out + "Coastline: {0}\n".format(self.coastline)
			out = out + "Station Code: {0}\n".format(self.stationcode)
			out = out + "Station Flag: {0}".format(self.stationflag)
			return(out)
		
	# Do tests to make sure dataDir exists, is a directory, and contains
	# everything that is expected to be there
	
	# Load the catalogue file
	psmsl_sites = []
	f = open(os.path.join(dataDir, "filelist.txt"), 'r')
	for line in f:
		
		# Get the metadata for this site
		(id,lat,lon,name,coastline,stationcode,stationflag) = line.split(';')

		# Initialize variables to hold data
		years = []
		heights = []
		interpolated = []
		isMt1 = []
		dataflag = []
		
		# Load the site's data file
		if(metaonly):
			# Store everything in an instance of the Sites class
			psmsl_sites.append(Site(id, lat, lon, name, coastline, stationcode, stationflag))
		else:
			df = open(os.path.join(dataDir, "data", id.strip()+".rlrdata"), 'r')
		
			# Extract the data from the file
			for dline in df:
				(tyear, theight, tinterp, tdflag) = dline.split(';')
				years.append(int(tyear))
				heights.append(np.nan) if int(theight) == -99999 else heights.append(int(theight))
				interpolated.append(tinterp == 'Y')
				isMt1.append(tdflag[1] == '1')
				dataflag.append(tdflag[2] == '1')
			df.close()
			
			# Store everything in an instance of the Sites class
			psmsl_sites.append(Site(id, lat, lon, name, coastline, stationcode, stationflag, years, heights, interpolated, dataflag, isMt1))

	f.close()
	
	return(psmsl_sites)