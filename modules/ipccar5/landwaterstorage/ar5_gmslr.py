# GMSLR projection program used for IPCC WG1 AR5
# Translated from IDL to Python 2.7 by Jonathan Gregory 23.10.19

import os,os.path,fnmatch
import cf
import numpy

class ProjectionError(Exception):
  pass

def mSLEoGt():
# Conversion factor for Gt to m SLE
  return 1e12/3.61e14*1e-3
  
def endofhistory():
  return 2006

def vlikely_range(data):
# Compute median and 5-95% range for the first (or only) axis of data.
# Return array (stat[,dim1,...]), where stat=0,1,2 for 50-,5-,95-percentile,
# and dim1,... are any remaining axes of data.
# NB model 5-95% range is judged to be "likely" for the AR5 projections
# data -- array-like
  return numpy.percentile(data,[50,5,95],0)

def actual_range(data):
# Compute mean and actual range for the first (or only) axis of data
# Return array (stat[,dim1,...]), where stat=0,1,2 for mean,minimum,maximum
# of data, and dim1,... are any remaining axes of data.
# data -- array-like
  return numpy.array([numpy.mean(data,0),numpy.amin(data,0),\
    numpy.amax(data,0)])

def dant():
# m SLE from Antarctica during 1996 to 2005 according to AR5 chapter 4
  return (2.37+0.13)*1e-3

def dgreen():
# m SLE from Greenland during 1996 to 2005 according to AR5 chapter 4
  return (3.21-0.30)*1e-3

def fgreendyn():
# Fraction of SLE from Greenland during 1996 to 2005 assumed to result from
# rapid dynamical change, with the remainder assumed to result from SMB change
  return 0.5

def report(quantity,field,output=None,uniform=False):
# Report the likely range of a projected quantity in the last timestep and
# optionally save the timeseries of likely range and median as CF-netCDF files
# quantity -- str, printed name of quantity
# field -- cf.Field, containing the data of the quantity, assumed to have time
#   as its last dimension
# output -- str, optional, name of directory for netCDF output files; if output
#   is not supplied, no files are written
# uniform -- bool, optional, indicates that the quantity has a uniform
#   distribution, for which the likely range is given by the extrema; by
#   default False, in which case the likely range is 5-95%

# Reshape as two-dimensional with time as the second dimension
  nyr=field.shape[-1]
  data=field.array.reshape(field.size/nyr,-1)

# Print likely range
  vformat="%20s %6.3f [%6.3f to %6.3f]"
  if uniform:
    datarange=actual_range(data)
  else:
    datarange=vlikely_range(data)
  print(vformat%tuple([quantity]+list(datarange[:,-1])))
  
# Optionally write output files
  if output:
    statfield=cf.Field(data=cf.Data(numpy.empty(nyr)))
    if quantity=="GMSLR":
      statfield.standard_name="global_average_sea_level_change"
    elif quantity=="expansion":
      statfield.standard_name="global_average_thermosteric_sea_level_change"
    else:
      statfield.long_name="GMSLR contribution from "+\
        dict(glacier="decrease of glacier mass",
        greensmb="decrease of Greenland ice sheet mass due to change in SMB",
        greendyn="decrease of Greenland ice sheet mass "+\
        "due to rapid dynamical change",
        greennet="decrease of Greenland ice sheet mass",
        antsmb="decrease of Antarctic ice sheet mass due to change in SMB",
        antdyn="decrease of Antarctic ice sheet mass "+\
        "due to rapid dynamical change",
        antnet="decrease of Antarctic ice sheet mass",
        landwater="decrease of land water storage",
        sheetdyn="decrease of ice sheet mass due to rapid dynamical change")\
        [quantity]
    statfield.insert_dim(field.dim('T'))
    statfield.unit='m'
    statfield.ncvar=quantity
    stats=dict(mid=0,lower=1,upper=2)
    for stat in stats:
      statfield.data[:]=datarange[stats[stat],:]
      cf.write(statfield,output+quantity+"_"+stat+".nc")

def project(input=None,scenarios=None,output=None,**kwargs):
# input -- str, path to directory containing input files. The directory should
#   contain files named SCENARIO_QUANTITY_STATISTIC.nc, where QUANTITY is
#   temperature or expansion, and STATISTIC is mean or sd. Each file contains
#   one field, having a single dimension of time.
# scenarios -- list of str, scenarios for which projections are to be made,
#   by default all those represented in the input directory
# output, str, optional -- path to directory in which output files are to be
#   written. It is created if it does not exist. No files are written if this
#   argument is omitted.
# seed -- optional, for numpy.random, default zero
# nt -- int, optional, number of realisations of the input timeseries for each
#   scenario, default 450
# nm -- int, optional, number of realisations of components and sum for each
#   realisation of the input timeseries, default 1000

# Check input directory
  if input is None:
    raise ProjectionError('input must be specified')
  input=os.path.expandvars(os.path.expanduser(input))
  if not os.path.isdir(input):
    raise ProjectionError('input must be an existing directory')

  if scenarios is None:
# Obtain list of scenarios from the input filenames
    bname=fnmatch.filter(os.listdir(input),'*_*.nc')
    scenarios=[tname.split('_',1)[0] for tname in bname]
    scenarios=sorted(list(set(scenarios)))
  elif isinstance(scenarios,str): # basestring
    scenarios=[scenarios]

  for scenario in scenarios:
    print(scenario)
    project_scenario(input,scenario,output=output,**kwargs)

def project_scenario(input,scenario,output=None,seed=0,nt=450,nm=1000):
# Make GMSLR projection for the specified single scenario
# Arguments are all the same as project() except for:
# scenario -- str, name of the scenario

  if not isinstance(scenario,str):
    raise ProjectionError('scenario must be a single string')

  numpy.random.seed(seed)

  startyr=endofhistory() # year when the timeseries for integration begin
    
# Read the input fields of ensemble statistics for temperature and expansion
# into txin, whose four elements are temperature mean, temperature sd,
# expansion mean, expansion sd. Check that each field is one-dimensional in
# time, that the mean and sd fields for each quantity have equal time axes,
# that temperature applies to calendar years (indicated by its time bounds),
# that expansion applies at the ends of the calendar years of temperature, and
# that there is no missing data.
  quantin=['temperature','expansion'] # input quantities
  it=0; ix=1 # indices to input quantities
  nqi=len(quantin)
  statin=['mean','sd'] # input statistics
  txin=[]
  for quant in quantin:
    for stat in statin:
      key=quant+'_'+stat
      file=input+'/'+scenario+'_'+key+'.nc'
      file=os.path.expandvars(os.path.expanduser(file))
      if not os.path.isfile(file):
        raise ProjectionError('missing input file: '+file)
      field=cf.read(file)[0]
      if field.ndim!=1:
        raise ProjectionError('field is not one-dimensional in file '+file)
      if not field.dims().values()[0].T:
        raise ProjectionError('field does not have a time axis in file '+file)
      if any(field.data==field._FillValue):
        raise ProjectionError('missing data is not allowed in file '+file)
      field.override_units('1',i=True)
      txin.append(field)
  for ii in [it,ix]:
    if not txin[ii*2].dim('T').equals(txin[ii*2+1].dim('T')):
      raise ProjectionError('time axes of mean and sd fields disagree for '+\
        quantin[ii]+' in scenario '+scenario)
  ttime=txin[0].dim('T')
  tbounds=ttime.bounds
  tupper=ttime.upper_bounds
  if (tbounds.month!=1 or tbounds.day!=1 \
    or tbounds.hour!=0 or tbounds.minute!=0 or tbounds.second!=0 \
    or tupper.year!=(ttime.lower_bounds.year+1)).any():
    raise ProjectionError('temperature values must be for calendar years')
  if (tupper.year[0]-1)!=startyr:
    raise ProjectionError('temperature must begin at '+str(startyr))
  if tupper.year[-1]>2100:
    raise ProjectionError('temperature input must not go beyond 2100')
  time=txin[2].dim('T') # expansion input supplies the output time coords
  nyr=time.size
  if (time!=tupper).any():
    raise ProjectionError('expansion must be for the ends of calendar years')

# Integrate temperature to obtain K yr at ends of calendar years
  itin=[txin[0].copy(),txin[1].copy()]
  for field in itin:
    field.data=cf.Data(numpy.cumsum(field.data))
    field.remove_item('T','d')
    field.insert_dim(txin[2].dim('T'))

# Generate a sample of perfectly correlated timeseries fields of temperature,
# time-integral temperature and expansion, each of them [realisation,time]
  z=numpy.random.standard_normal(nt)
  z=cf.Field(data=cf.Data(z))
  climdim=cf.DimensionCoordinate(data=cf.Data(numpy.arange(nt)),\
    properties=dict(standard_name='climate_realization'))
  z.insert_dim(climdim)
# For each quantity, mean + standard deviation * normal random number
  zt=txin[0]+txin[1]*z
  zx=txin[2]+txin[3]*z
  zit=itin[0]+itin[1]*z

# Create a cf.Field with the shape of the quantities to be calculated
# [component_realization,climate_realization,time]
  template=cf.Field(data=cf.Data(numpy.full([nm,nt,nyr],numpy.nan)))
  template.insert_dim(txin[2].dim('T'),key='dim2')
  template.insert_dim(climdim,key='dim1')
  template.insert_dim(cf.DimensionCoordinate(data=cf.Data(numpy.arange(nm)),\
    properties=dict(standard_name='component_realization')),key='dim0')

# Obtain ensembles of projected components as cf.Field objects and add them up
  expansion=zx
  glacier=project_glacier(itin[0],zit,template)
  greensmb=project_greensmb(zt,template)
  greendyn=project_greendyn(scenario,template)
  greennet=greensmb+greendyn
  fraction=numpy.random.rand(nm*nt) # correlation between antsmb and antdyn
  antsmb=project_antsmb(zit,template,fraction=fraction)
  antdyn=project_antdyn(template,fraction=fraction)
  antnet=antdyn+antsmb
  sheetdyn=greendyn+antdyn
  landwater=project_landwater(template)
  gmslr=expansion+glacier+greensmb+greendyn+antnet+landwater

# Report the range of the final year and write output files if requested
  if output:
    output=os.path.expandvars(os.path.expanduser(output))
    if not os.path.isdir(output): os.mkdir(output)
    output=output+"/"+scenario+"_"
  report("expansion",expansion,output)
  report("glacier",glacier,output)
  report("greensmb",greensmb,output)
  report("antsmb",antsmb,output)
  report("greendyn",greendyn,output,uniform=True)
  report("antdyn",antdyn,output,uniform=True)
  report("landwater",landwater,output,uniform=True)
  report("GMSLR",gmslr,output)
  report("greennet",greennet,output)
  report("antnet",antnet,output)
  report("sheetdyn",sheetdyn,output)

  return

def project_glacier(it,zit,template):
# Return projection of glacier contribution as a cf.Field
# it -- cf.Field, time-integral of median temperature anomaly timeseries
# zit -- cf.Field, ensemble of time-integral temperature anomaly timeseries
# template -- cf.Field with the required shape of the output

  startyr=int(template.dim('T').year.data[0])-1

  dmzdtref=0.95 # mm yr-1 in Marzeion's CMIP5 ensemble mean for AR5 ref period
  dmz=dmzdtref*(startyr-1996)*1e-3 # m from glacier at start wrt AR5 ref period
  cvgl=0.20 # random methodological error
  glmass=412.0-96.3 # initial glacier mass, used to set a limit, from Tab 4.2
  glmass=1e-3*glmass # m SLE

  nr=template.shape[0]
  glparm=[dict(name='Marzeion',factor=4.96,exponent=0.685),\
    dict(name='Radic',factor=5.45,exponent=0.676),\
    dict(name='Slangen',factor=3.44,exponent=0.742),\
    dict(name='Giesen',factor=3.02,exponent=0.733)]
  ngl=len(glparm) # number of glacier methods
  if nr%ngl:
    raise ProjectionError('number of realisations '+\
      'must be a multiple of number of glacier methods')
  nrpergl=nr/ngl # number of realisations per glacier method
  r=cf.Field(data=cf.Data(numpy.random.standard_normal(nr)))
  r.insert_dim(template.dim('dim0'))

# Make an ensemble of projections for each method
  glacier=template.copy()
  for igl in range(ngl):
# glacier projection for this method using the median temperature timeseries
    mgl=project_glacier1(it,glparm[igl]['factor'],glparm[igl]['exponent'])
# glacier projections for this method with the ensemble of timeseries
    zgl=project_glacier1(zit,glparm[igl]['factor'],glparm[igl]['exponent'])
    ifirst=igl*nrpergl
    ilast=ifirst+nrpergl
    glacier[ifirst:ilast,...]=zgl+mgl*r[ifirst:ilast]*cvgl

  glacier+=dmz
  glacier.where(glacier>glmass,glmass,i=True)

  return glacier

def project_glacier1(it,factor,exponent):
# Return projection of glacier contribution by one glacier method
  scale=1e-3 # mm to m
  return scale*factor*(it.where(it<0,0)**exponent)

def project_greensmb(zt,template):
# Return projection of Greenland SMB contribution as a cf.Field
# zt -- cf.Field, ensemble of temperature anomaly timeseries
# template -- cf.Field with the required shape of the output

  dtgreen=-0.146 # Delta_T of Greenland ref period wrt AR5 ref period  
  fnlogsd=0.4 # random methodological error of the log factor
  febound=[1,1.15] # bounds of uniform pdf of SMB elevation feedback factor

  nr=template.shape[0]
# random log-normal factor
  fn=numpy.exp(numpy.random.standard_normal(nr)*fnlogsd)
# elevation feedback factor
  fe=numpy.random.sample(nr)*(febound[1]-febound[0])+febound[0]
  ff=cf.Field(data=cf.Data(fn*fe))
  ff.insert_dim(template.dim('dim0'))
  
  ztgreen=zt-dtgreen
  greensmbrate=fettweis(ztgreen)*ff

  greensmb=template.copy()
  greensmb.data[:]=numpy.cumsum(greensmbrate.data,axis=2)[:]
  greensmb+=(1-fgreendyn())*dgreen()

  return greensmb

def fettweis(ztgreen):
# Greenland SMB in m yr-1 SLE from global mean temperature anomaly
# using Eq 2 of Fettweis et al. (2013)
  return (71.5*ztgreen+20.4*(ztgreen**2)+2.8*(ztgreen**3))*mSLEoGt()

def project_antsmb(zit,template,fraction=None):
# Return projection of Antarctic SMB contribution as a cf.Field
# zit -- cf.Field, ensemble of time-integral temperature anomaly timeseries
# template -- cf.Field with the required shape of the output
# fraction -- array-like, random numbers for the SMB-dynamic feedback

  antsmb=template.copy()
  nr,nt,nyr=antsmb.shape

# The following are [mean,SD]
  pcoK=[5.1,1.5] # % change in Ant SMB per K of warming from G&H06
  KoKg=[1.1,0.2] # ratio of Antarctic warming to global warming from G&H06

# Generate a distribution of products of the above two factors
  pcoKg=(pcoK[0]+numpy.random.standard_normal([nr,nt,1])*pcoK[1])*\
    (KoKg[0]+numpy.random.standard_normal([nr,nt,1])*KoKg[1])
  meansmb=1923 # model-mean time-mean 1979-2010 Gt yr-1 from 13.3.3.2
  moaoKg=-pcoKg*1e-2*meansmb*mSLEoGt() # m yr-1 of SLE per K of global warming

  if fraction is None:
    fraction=numpy.random.rand(nr,nt)
  elif fraction.size!=nr*nt:
    raise ProjectionError('fraction is the wrong size')
  else:
    fraction.shape=(nr,nt,1)

  smax=0.35 # max value of S in 13.SM.1.5
  ainterfactor=1-fraction*smax
  
  antsmb.data[:]=moaoKg*ainterfactor*zit.array.reshape(1,nt,-1)[:]
  return antsmb

def project_greendyn(scenario,template):
# Return projection of Greenland rapid ice-sheet dynamics contribution
# as a cf.Field
# scenario -- str, name of scenario
# template -- cf.Field with the required shape of the output

# For SMB+dyn during 2005-2010 Table 4.6 gives 0.63+-0.17 mm yr-1 (5-95% range)
# For dyn at 2100 Chapter 13 gives [20,85] mm for rcp85, [14,63] mm otherwise

  if scenario in ['rcp85','ssp5_85']:
    finalrange=[0.020,0.085]
  else:
    finalrange=[0.014,0.063]
  return time_projection(0.63*fgreendyn(),\
    0.17*fgreendyn(),finalrange,template)+fgreendyn()*dgreen()

def project_antdyn(template,fraction=None):
# Return projection of Greenland rapid ice-sheet dynamics contribution
# as a cf.Field
# template -- cf.Field with the required shape of the output
# fraction -- array-like, random numbers for the dynamic contribution

# For SMB+dyn during 2005-2010 Table 4.6 gives 0.41+-0.24 mm yr-1 (5-95% range)
# For dyn at 2100 Chapter 13 gives [-20,185] mm for all scenarios

  return time_projection(0.41,0.20,[-0.020,0.185],template,fraction=fraction)+\
    dant()

def project_landwater(template):
# Return projection of land water storage contribution as a cf.Field

# The rate at start is the one for 1993-2010 from the budget table.
# The final amount is the mean for 2081-2100.
  nyr=2100-2081+1 # number of years of the time-mean of the final amount

  return time_projection(0.38,0.49-0.38,[-0.01,0.09],template,nyr)

def time_projection(startratemean,startratepm,finalrange,template,\
  nfinal=1,fraction=None):
# Return projection of a quantity which is a quadratic function of time
# in a cf.Field.
# startratemean, startratepm -- rate of GMSLR at the start in mm yr-1, whose
#   likely range is startratemean +- startratepm
# finalrange -- two-element list giving likely range in m for GMSLR at the end
# template -- cf.Field with the required shape of the output
# nfinal -- int, optional, number of years at the end over which finalrange is
#   a time-mean; by default 1 => finalrange is the value for the last year
# fraction -- array-like, optional, random numbers in the range 0 to 1,
#   by default uniformly distributed

  nr,nt,nyr=template.shape
  if fraction is None:
    fraction=numpy.random.rand(nr,nt)
  elif fraction.size!=nr*nt:
    raise ProjectionError('fraction is the wrong size')
  fraction=cf.Field(data=cf.Data(fraction.reshape(nr,nt)))
# The following will be correct regardless of the ordering of the
# component_realization and climate_realization axes
  fraction.insert_dim(template.dim('dim0'),key='dim0')
  fraction.insert_dim(template.dim('dim1'),key='dim1')

# For terms where the rate increases linearly in time t, we can write GMSLR as
#   S(t) = a*t**2 + b*t
# where a is 0.5*acceleration and b is start rate. Hence
#   a = S/t**2-b/t
  momm=1e-3 # convert mm yr-1 to m yr-1
  startrate=(startratemean+\
    startratepm*numpy.array([-1,1],dtype=numpy.float))*momm
  finalyr=numpy.arange(nfinal)-nfinal+nyr+1 # last element ==nyr
# If nfinal=1, the following is equivalent to
# numpy.array(finalrange,dtype=numpy.float)/nyr**2-startrate/nyr
  acceleration=(numpy.array(finalrange,dtype=numpy.float)-\
    startrate*finalyr.mean())/(finalyr**2).mean()

# Create a field of elapsed time in years
  tdim=template.dim('T')
  time=tdim.year.data
  time=time-time[0]+1 # years since start
  time=cf.Field(data=time)
  time.insert_dim(tdim)

# Calculate two-element list containing fields of the minimum and maximum
# timeseries of projections, then calculate random ensemble within envelope
  range=[float(acceleration[i])*(time**2)+float(startrate[i])*time \
    for i in [0,1]]
  projection=template.copy()
  projection.data[:]=(range[0]*(1-fraction)+range[1]*fraction)[:]

  return projection
