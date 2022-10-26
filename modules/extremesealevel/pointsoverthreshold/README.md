# extremesealevel/pointsoverthreshold

This module fits a General Pareto Distribution on preprocessed GESLA2 data. It also calculates the associated covariance matrix. Using input from the total module, it generates samples of local msl change and GPD parameters. From those samples it
calculates historical and future return curves at user defined return periods. The return curves are used to calculate the amplification factor and allowance 
for a given station at user defined percentiles. The analysis is based on the MATLAB code of Thomas Frederikse used for SROCC and LocalizeSL (Buchanan et al. 2016). A Peak-Over-Threshold is used. 
Above the threshold a Pareto distribution is used to model the extremes. Below, a Gumbel distribution is assumed and cut off at MHHW. MHHW is calculated as the
long-term mean of daily maxima.

Code by Tim Hermans
