# ar5/glaciers

This module implements the glacier projection approach adopted by IPCC AR5 WG1. It was translated from IDL to Python 2.7 by Jonathan Gregory 23.10.19 and adapted for use in FACTS by Gregory Garner 20 November 2019.

See IPCC AR5 WG1 13.SM.1.4 - 13.SM.1.6 for the original implementation. This implementation also includes as options calibration parameters derived from the GlacierMIP and GlacierMIP2 studies as described in AR6 WG1 9.SM.4.5. As described therein:
                          
The glacier contribution is the integral of $fI(t)^ρ$, where $I(t)$ is the time integral of GSAT from 2006 to time t in degrees Celsius year, and the constants f and ρ are calibrated for each glacier model. The spread of the results around this median projection has a coefficient of variation (standard deviation divided by the mean) σ which is determined on a per-model basis. This variation is incorporated by taking for each Monte Carlo sample a normally distributed random number. This number is multiplied by the time-dependent standard deviation and added to the sample. All models are equally weighted.

