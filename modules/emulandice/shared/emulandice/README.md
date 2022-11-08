# emulandice ///////////////////////////////
Emulator of multi-model land ice projections (Edwards et al., 2021, Nature).

Input data files (inst/extdata/):

Calibrated with ISMIP6 and GlacierMIP model projections: 20201106_SLE_SIMULATIONS.csv.

Forced with GSAT projections from:

- FAIR ensemble of McKenna et al. (2020): 20201105_CLIMATE_FORCING.csv, or

- IPCC AR6 Chapter 7 two layer model ensemble: 20210215_CLIMATE_FORCING_IPCC.csv

Prior distributions on ice sheet parameters derived from:

- Greenland retreat: data_for_tamsin.txt                

- Antarctic basal melt: output_gamma0_NonLocal_MeanAnt.dat and output_gamma0_NonLocal_PIGL.dat

Below is the list of commands to run the analyses. OK to have this message on running: 
"There were 50 or more warnings (use warnings() to see the first 50)"


#______________________________________________________
# MAIN PROJECTIONS

# Edwards et al. (2021) ////////////////////////////////////

# Main projections: 2100
Set N_temp because default for timeslice is 1000 (for sensitivity tests)

main(N_temp = 5000L)

# Main projections: full timeseries
main("timeseries")

# Antarctic projections: subset of models in LARMIP and ISMIP6
main(N_temp = 5000L, select_ism == "larmip")

# IPCC AR6 Chapter 9 ///////////////////////////////////////

# Main projections: 2100
Set N_temp as above

main(dataset = "IPCC", N_temp = 10000L)  

# Main projections: full timeseries
Run 2-3 SSPs at once for timeseries because 4x larger ensembles than main FAIR projections

main("timeseries", dataset = "IPCC", fair_ssps = c("SSP126", "SSP585"))   

main("timeseries", dataset = "IPCC", fair_ssps = c("SSP119", "SSP245", "SSP370"))   

# Antarctic projections: subset of models in LARMIP and ISMIP6
main(dataset = "IPCC", N_temp = 10000L, select_ism = "larmip")

main("timeseries", dataset = "IPCC", fair_ssps = c("SSP126", "SSP245", "SSP585"), select_ism = "larmip")

#______________________________________________________
# COVARIANCE FUNCTION TESTS
Each sets same covariance function for all regions

# Default covariance (plus model comparison)
main("SA", do_model_comp = TRUE, do_covar_fn = "matern_5_2")

# Alternative covariances
main("SA", do_covar_fn = "matern_3_2")

main("SA", do_covar_fn = "pow_exp", do_covar_alpha = 1.9)

main("SA", do_covar_fn = "pow_exp", do_covar_alpha = 1.0)

main("SA", do_covar_fn = "pow_exp", do_covar_alpha = 0.1)

#________________________________________________________
# SENSITIVITY ANALYSIS

# Main sensitivity analysis 
main("SA")

# Simulation data plots
main("sim_only")

#________________________________________________________
# SENSITIVITY TESTS: Extended Data Figure 4

# S2 CMIP6 ensemble for GSAT
main("CMIP6")

# S3 Fixed GSAT
main("fixed_T")

# S4 Fixed ice sheet retreat and basal melt
main("fixed_melt")

# S5-8 Basal melt priors
main("gamma0_MeanAnt", ice_sources = "AIS")

main("gamma0_PIGL", ice_sources = "AIS")

main("gamma0_unif", ice_sources = "AIS")

main("gamma0_unif_high", ice_sources = "AIS")

# S9-10 Ice shelf collapse priors
main(collapse_prior = "off")

main(collapse_prior = "on")

# S11 Risk-averse projections: SA, 2100 and full timeseries
main("SA", risk_averse = TRUE)

main(risk_averse = TRUE, N_temp = 5000L)

main("timeseries", risk_averse = TRUE)

#________________________________________________________
# ROBUSTNESS TESTS: Extended Data Figure 5

# SA are to check nothing strange happens using smaller datasets
# (check Tdep, meltdep plots)

# R2 Minimum ice sheet resolution
main("SA", min_res_is = c(8,30))

main("high_res")

# R3 More balanced ice sheet model designs
main("SA", select_ism = "balanced")

main(select_ism = "balanced")

# R4 One ice sheet model
main("SA", select_ism = "single")

main(select_ism = "single")

# R5 Further Antarctic ice sheet model (low sensitivity)
main("SA", select_ism = "single_low")

main(select_ism = "single_low")

# R6 Subset of Antarctice ice sheet models with high sensitivity
main("SA", select_ism = "high_sensitivity")

main(select_ism = "high_sensitivity")

# R7 Subsets of climate models driving Antarctic ice sheet models
main("SA", select_gcm = "high_five")

main(select_gcm = "high_five")

main("SA", select_gcm = "high_two")

main(select_gcm = "high_two")

# R8 History matching with Antarctic observational trend
main("SA", history_match = TRUE)

main(history_match = TRUE)

# R9 Exclude Antarctic models not using parameterisation
main("SA", exclude_open = TRUE)

main(exclude_open = TRUE)

# R10 Use higher value for missing basal melt 
main("SA", impute_high = TRUE)

main(impute_high = TRUE)

