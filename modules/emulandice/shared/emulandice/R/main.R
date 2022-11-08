#' Environment for useful variables.
e <- new.env(parent = emptyenv())

# TO DO
# Check all function arguments written to log file
# Check behaviour if N_temp < number off models in prior

# start --------------------------------------

print('____________________________________________________')
print('Welcome to emulandice')
print(' ', quote = FALSE)
print('Can ignore the warnings from optimize')
print('____________________________________________________')

# No arguments: runs defaultx 2100
# Change expt: will alter arg(s)
# Change args: less well tested...
# N_temp = default for tests etc; use 5000 for projections at 2100
main <- function(expt = "default",
                 ice_sources = c("GrIS", "AIS", "Glaciers"),
                 years = 2100,
                 dataset = "main",
                 N_temp = 1000L,
                 N_FACTS =  2237L,
                 outdir = "results",
                 temp_prior = "FAIR",
                 fair_ssps = NA,
                 mean_temp = FALSE,
                 gamma0_prior = "joint",
                 mean_melt = FALSE,
                 collapse_prior = "both",
                 risk_averse = FALSE,
                 min_res_is = c(NA, NA),
                 select_ism = "all",
                 select_gcm = "all",
                 history_match = FALSE,
                 exclude_open = FALSE,
                 impute_high = FALSE,
                 do_model_comp = FALSE,
                 do_covar_fn = NA,
                 do_covar_alpha = NA,
                 packagename = "emulandice") {

  #' Main analysis steering function
  #' @param expt Analysis to run: e.g. "SA", "timeseries", "decades"
  #' @param ice_sources Ice sources: GrIS, AIS, Glaciers
  #' @param years Year(s) to predict: default is 2100, time series is 2015:2100
  #' @param dataset Forcing dataset: 2019, main, IPCC, FACTS
  #' @param N_temp Number of climate values in prior: 501 code testing, 1000L default for tests, 5000L projections, over-ridden for timeseries or decades
  #' @param N_FACTS Number of  FAIR time series samples in forcing file passed by FACTS
  #' @param outdir output directory
  #' @param temp_prior Climate ensemble for prior: FAIR, CMIP6
  #' @param fair_ssps Restrict FAIR SSPs run: NA; or e.g. c("SSP126", "SSP585"); must be set for IPCC timeseries or decades runs
  #' @param mean_temp Use mean temperature value or ice sheets: T/F
  #' @param gamma0_prior Gamma0 prior distribution for AIS: "joint", "MeanAnt", "PIGL", "unif", "unif_high"
  #' @param mean_melt Use mean kappa/gamma0 value for ice sheets: T/F
  #' @param collapse_prior Ice shelf collapse prior: both, on, off
  #' @param risk_averse Risk-averse settings: T/F
  #' @param min_res_is Set minimum spatial resolution for ice sheets in km (GIS, AIS): default c(8, 32)
  #' @param select_ism Select subset of ice sheet models by defined list: balanced, single, single_low, high_sensitivity, larmip
  #' @param select_gcm Select subset of gcms by defined list: high_five, high_two
  #' @param history_match Select subset of ice sheet models using IMBIE observations: T/F
  #' @param exclude_open Exclude open parameterisation Antarctic ice sheet models: T/F
  #' @param impute_high Use higher value of basal melt for open parameterisation Antarctic models: T/F
  #' @param do_model_comp Run stepwise model comparison: T/F
  #' @param do_covar_fn Set fixed covariance function for all regions: matern_5_2, matern_3_2, pow_exp
  #' @param do_covar_alpha Set fixed pow_exp exponent for all regions: 0.1, 1.0, 1.9
  #' @param packagename Set package name

  # EXPERIMENT OPTIONS: each changes one of the other options
  stopifnot(expt %in% c("default", "timeseries", "decades", # --> projections for timeslice or full time series
                        "sim_only", "SA", # --> plot simulations or param dep, NOT projections
                        "high_res", # ice sheet model selection
                        "CMIP6", "gamma0_MeanAnt", "gamma0_PIGL", "gamma0_unif", "gamma0_unif_high",
                        "fixed_melt", "fixed_T")) # -> priors

  # options --------------------------------------

  e$packagename <- packagename

  # Ice sources: default is to predict for all land ice
  e$ice_source_list <- ice_sources

  # Ice sheet only tests
  if ( expt == "sim_only" || expt == "high_res" || select_ism != "all" ||
       ( ! is.na(min_res_is[1]) || ! is.na(min_res_is[2]) ) ) e$ice_source_list <- c("GrIS", "AIS")

  # Antarctic only tests
  if ( collapse_prior != "both" ||
       select_ism == "single_low" || select_ism == "high_sensitivity" || select_ism == "larmip" ||
       select_gcm != "all" ||
       history_match == TRUE ||
       exclude_open == TRUE ||
       impute_high == TRUE ) e$ice_source_list <- "AIS"

  # Number of IPCC runs for each SSP
  N_IPCC <- 2237L # v.0.1.0 i.e. 20210215_CLIMATE_FORCING_IPCC.csv (was 2000 in v.0.0.0)
 
  # Number of T/melt samples in 2100 projections and SA
  # If equal to number of FAIR projections, uses each one, otherwise samples
  # 501 for testing, 1000 for SA tests, ~2000 for FAIR 2LM IPCC, 5000 for FAIR main projections
  stopifnot(N_temp %in% c(501L, 1000L, N_IPCC, N_FACTS, 5000L, 10000L))

  # Collapse prior
  stopifnot(collapse_prior %in% c("both", "on", "off"))

  # Model comparison: stepwise BIC for ice sheets
  # Covariance comparison: set single covariance function for all regions so quick to test
  if (!is.na(do_covar_fn)) {
    stopifnot(do_covar_fn %in% c("matern_5_2", "matern_3_2", "pow_exp"))
    if ( do_covar_fn == "pow_exp") {
      stopifnot( !is.na(do_covar_alpha) )
      stopifnot( do_covar_alpha > 0 && do_covar_alpha < 2 )
    }
  }

  # Dummy variables
  e$add_dummy <- ifelse(select_ism %in% c("balanced", "single"), "none", "melt") # Only adds for Greenland now
  e$dummy_melt_pred <- 0 # for predicting 'open' type models (1), standard (0) or a random mix (0.5)
  stopifnot(e$add_dummy %in% c("model", "group", "melt", "none"))
  stopifnot(e$dummy_melt_pred %in% c(0,0.5,1))

  # YEARS TO PREDICT i.e. output to CSV files
  e$years_pred <- years # some years added later: e$years_data
  if (expt == "timeseries") e$years_pred <- 2016:2100
  if (expt == "decades") e$years_pred <- seq(from=2020, to=2100, by=10)
  stopifnot( min(e$years_pred) >= 2016 && max(e$years_pred) <= 2100 )

  # PARAMETER SAMPLING
  annual_resample <- FALSE # sample melt annually (T) or once at start of time series (F, default)

  # Basal melt by region
  one_sample_AIS <- risk_averse # If TRUE, use same gamma for all regions

  # Don't include emulator uncertainty in Monte Carlo sample
  no_emulator_uncertainty_mc <- FALSE # default = FALSE!

  # Quantiles to predict
  q_list <- c(0.50, 0.05, 0.95, 0.17, 0.83, 0.25, 0.75 )

  # DATASET SELECTION

  # Risk-averse projections
  if (risk_averse) {
    select_ism <- "high_sensitivity"
    select_gcm <- "high_five"
  }

  # Select GCM forcing
  stopifnot(select_gcm %in% c("all", "high_five", "high_two"))
  e$select_gcm <- ifelse(select_gcm == "all", FALSE, TRUE)
  e$cond_gcm <- list()

  # Greenland: no GCM selection
  e$cond_gcm[["GrIS"]] <- NA

  # Antarctica: 5 or 2 GCMs
  if (select_gcm == "high_five") e$cond_gcm[["AIS"]] <- c("HadGEM2-ES", "UKESM1-0-LL", "MIROC-ESM-CHEM", "NorESM1-M", "CCSM4")
  if (select_gcm == "high_two") e$cond_gcm[["AIS"]] <- c("NorESM1-M", "CCSM4")

  # Select ISM
  stopifnot(select_ism %in% c("all", "balanced", "single", "single_low", "high_sensitivity", "larmip"))
  e$select_ism <- ifelse(select_ism == "all", FALSE, TRUE)
  e$cond_ism <- list()

  # Default if unchanged below
  e$cond_ism[["GrIS"]] <- NA
  e$cond_ism[["AIS"]] <- NA

  # More balanced design
  # 10 GrIS models with 14 or more runs
  if (select_ism == "balanced") e$cond_ism[["GrIS"]] <- c("AWI__ISSM1", "AWI__ISSM2", "AWI__ISSM3",
                                                          "ILTS_PIK__SICOPOLIS1", "ILTS_PIK__SICOPOLIS2", "IMAU__IMAUICE1",
                                                          "JPL__ISSM", "LSCE__GRISLI2", "NCAR__CISM", "VUB__GISMHOMv1")

  # 1 model with most runs
  if (select_ism == "single") e$cond_ism[["GrIS"]] <- "IMAU__IMAUICE1"

  # Checks against lists copied from write_sealevel.R
  stopifnot( is.na(e$cond_ism[["GrIS"]]) ||
               e$cond_ism[["GrIS"]] %in% c("AWI__ISSM1", "AWI__ISSM2", "AWI__ISSM3",
                                           "BGC__BISICLES", "GSFC__ISSM",
                                           "ILTS_PIK__SICOPOLIS1", "ILTS_PIK__SICOPOLIS2",
                                           "IMAU__IMAUICE1", "IMAU__IMAUICE2",
                                           "JPL__ISSM", "JPL__ISSMPALEO",
                                           "LSCE__GRISLI2", "MUN__GSM2601", "MUN__GSM2611",
                                           "NCAR__CISM", "UAF__PISM1", "UAF__PISM2",
                                           "UCIJPL__ISSM1", "UCIJPL__ISSM2", "VUB__GISMHOMv1", "VUW__PISM") )

  # More balanced design
  # 4 AIS models with most runs
  if (select_ism == "balanced") e$cond_ism[["AIS"]] <- c("ILTS_PIK__SICOPOLIS", "JPL1__ISSM", "LSCE__GRISLI", "NCAR__CISM")

  # LARMIP-2 comparison: groups/models common to both (including all variants of model)
  if (select_ism == "larmip") e$cond_ism[["AIS"]] <- c( "AWI__PISM1", "CPOM__BISICLES", "DOE__MALI",
                                                        "ILTS_PIK__SICOPOLIS", "IMAU__IMAUICE1", "IMAU__IMAUICE2",
                                                        "JPL1__ISSM", "LSCE__GRISLI",
                                                        "NCAR__CISM", "PIK__PISM1", "PIK__PISM2",
                                                        "UCIJPL__ISSM", "ULB__fETISh_16km", "ULB__fETISh_32km",
                                                        "VUB__AISMPALEO", "VUW__PISM")

  # Single models - high and low sensitivity, then 4 highest sensitivity
  if (select_ism == "single") e$cond_ism[["AIS"]] <- "ILTS_PIK__SICOPOLIS"
  if (select_ism == "single_low") e$cond_ism[["AIS"]] <- "LSCE__GRISLI"
  if (select_ism == "high_sensitivity") e$cond_ism[["AIS"]] <- c("ILTS_PIK__SICOPOLIS", "ULB__fETISh_16km",
                                                                 "ULB__fETISh_32km", "DOE__MALI")

  stopifnot( is.na(e$cond_ism[["AIS"]]) ||
               e$cond_ism[["AIS"]] %in% c( "AWI__PISM1", "CPOM__BISICLES", "DOE__MALI",
                                           "ILTS_PIK__SICOPOLIS", "IMAU__IMAUICE1", "IMAU__IMAUICE2",
                                           "JPL1__ISSM", "LSCE__GRISLI", "NCAR__CISM",
                                           "PIK__PISM1", "PIK__PISM2", "UCIJPL__ISSM",
                                           "ULB__fETISh_16km", "ULB__fETISh_32km", "UTAS__ElmerIce",
                                           "VUB__AISMPALEO", "VUW__PISM") )
  # Not used in main analysis
  e$cond_ism[["Glaciers"]] <- NA

  # Exclude open melt Antarctic
  e$exclude_open_AIS <- exclude_open

  # Use PIGL_med for imputed melt value instead of ensemble mean
  e$impute_high_melt_AIS <- impute_high

  # IMBIE Antarctic mass trend
  e$history_match <- history_match

  # RESOLUTION OF MODELS TO EXCLUDE
  # Excludes this grid cell size and larger
  # use e.g. 16 to exclude 16km upwards
  if (expt == "high_res") min_res_is <- c(8, 32)

  # Which Antarctic collapse prior?
  # default: 0.5 for random sampling from c(0,1); 0 or 1 for conditional on this
  if (collapse_prior == "both" ) collapse_prior <- 0.5
  if (collapse_prior == "on" || risk_averse == TRUE) collapse_prior <- 1
  if (collapse_prior == "off" ) collapse_prior <- 0
  stopifnot(collapse_prior %in% c(0,0.5,1))

  # Run on 2019 SLE CSV to check impact of corrected and extended dataset
  old_data <- FALSE # xxx improve: remove old_data capacity once code rewrite finished

  # FORCING csv: old, main or new 2LM forcing
  # stopifnot(dataset %in% c("2019", "main", "IPCC", "FACTS"))

  # This is used to override N_temp for timeseries runs, so that each is a trajectory
  N_FAIR <- N_FACTS
  if(dataset == "IPCC") N_FAIR <- N_IPCC

  # End of anything changed by hand
  #__________________________________________________________________________________________________
  #__________________________________________________________________________________________________

  # Random seed
  set.seed(2020)

  # process choices --------------------------------------

  imbie_range <- "2012_2017"
  stopifnot(imbie_range %in% c("1992_2017", "2012_2017"))

  # Final year of period (beginning 2015) to compare rate
  e$y_obs <- 2020 # used later in date selection
  n_yr <- e$y_obs - 2015 # not + 1 because 2015 value is zero

  # Number of sigma to exclude from mean
  e$n_sig <- 5

  # IMBIE rate in Gt/yr
  if (imbie_range == "1992_2017") e$imbie_obs <- c(-109, 56)
  if (imbie_range == "2012_2017") e$imbie_obs <- c(-219, 43)
  names(e$imbie_obs) <- c("mean", "sd")

  # Convert to cm/yr SLE contribution
  # Multiply by number of years to get contribution in y_obs
  e$imbie_obs <- n_yr * e$imbie_obs / (362.5*10)
  e$imbie_obs["mean"] = -e$imbie_obs["mean"]

  # Store resolution thresholds in list
  min_res <- list()
  min_res[["GrIS"]] = min_res_is[1]
  min_res[["AIS"]] = min_res_is[2]

  # ANTARCTIC MELT PRIOR
  if (expt == "gamma0_MeanAnt") gamma0_prior <- "MeanAnt"
  if (expt == "gamma0_PIGL" || risk_averse == TRUE) gamma0_prior <- "PIGL"
  if (expt == "gamma0_unif") gamma0_prior <- "unif"
  if (expt == "gamma0_unif_high") gamma0_prior <- "unif_high"
  stopifnot( gamma0_prior %in% c("joint", "MeanAnt", "PIGL", "unif", "unif_high") )

  # CLIMATE PRIOR
  if (expt == "CMIP6") temp_prior <- "CMIP6"

  # N_temp = N_FAIR for timeseries projections, otherwise use N_temp
  if (temp_prior == "FAIR" && expt == "timeseries") N_temp <- N_FAIR
  if (temp_prior == "FAIR" && expt == "decades") N_temp <- N_FAIR

  # Number of melt samples per temperature
  if (expt == "SA") N_melt_Tdep <- N_temp # T-dep plots
  else N_melt_Tdep <- 1L # Projections etc

  # Density estimate to get more climate values
  climate_prior_kde <- ifelse(temp_prior == "CMIP6", TRUE, FALSE)

  # Use mean melt instead of sampling distribution
  if (expt == "fixed_melt") mean_melt <- TRUE
  # Later: mean_melt_value[[is]] for value, e.g. default (kappa_50 and MeanAnt_50)

  # Use mean temp instead of sampling distribution?
  if (expt == "fixed_T") mean_temp <- TRUE
  stopifnot(mean_melt == FALSE || mean_temp == FALSE) # fixed temp and fixed melt doesn't work

  # SCENARIOS
  # Names match CSV entries
  scenario_list <- list()

  # Full possible lists
  if (dataset == "2019") scenario_list[["FAIR"]] <- c("SSP119", "SSP126", "SSP245", "SSP370", "SSP585")
  else if (dataset == "main") scenario_list[["FAIR"]] <- c("SSP119", "SSP126", "SSP245", "SSPNDC", "SSP370", "SSP585")
  else if (dataset == "IPCC") scenario_list[["FAIR"]] <- c("SSP119", "SSP126", "SSP245", "SSP370", "SSP585")
  else scenario_list[["FAIR"]] <- c("FACTS")

  # Do not allow running all SSPs: too slow
  if (dataset == "IPCC" && expt == "timeseries") {
    stopifnot( ! is.na(fair_ssps[1]) )
  }

  # Subselection of FAIR SSPs, if given
  if ( !is.na(fair_ssps[1]) ) {
    scenario_list[["FAIR"]] <- scenario_list[["FAIR"]][ which(scenario_list[["FAIR"]] %in% fair_ssps ) ]
    print("Selecting scenarios:")
    print(scenario_list[["FAIR"]])
    stopifnot(length(scenario_list[["FAIR"]]) >= 1)
  }

  scenario_list[["CMIP5"]] <- c("RCP26", "RCP85", "RCP45", "RCP60")
  scenario_list[["CMIP6"]] <- c("SSP126", "SSP245", "SSP370", "SSP585") # not enough models for SSP119

  # SSP/RCP names to expect for each ensemble
  # Needs to be same order as scenario_list above
  # XXX Could have done this as a lookup table
  e$scen_name_list <- list()

  # Full lists
  if (dataset == "2019") e$scen_name_list[["FAIR"]] <- c("SSP1-19", "SSP1-26", "SSP2-45", "SSP3-70", "SSP5-85")
  else if (dataset == "main") e$scen_name_list[["FAIR"]] <- c("SSP1-19", "SSP1-26", "SSP2-45", "NDCs", "SSP3-70", "SSP5-85")
  else if (dataset == "IPCC") e$scen_name_list[["FAIR"]] <- c("SSP1-19", "SSP1-26", "SSP2-45", "SSP3-70", "SSP5-85")
  else e$scen_name_list[["FAIR"]] <- c("FACTS")

  # Get subset of names if selecting SSPs
  if ( !is.na(fair_ssps[1]) ) {
    e$scen_name_list[["FAIR"]] <- e$scen_name_list[["FAIR"]][ which(scenario_list[["FAIR"]] == fair_ssps ) ]
  }

  e$scen_name_list[["CMIP5"]] <- c("RCP2.6", "RCP8.5", "RCP4.5", "RCP6.0")
  e$scen_name_list[["CMIP6"]] <- c("SSP1-26", "SSP2-45", "SSP3-70", "SSP5-85")


  # output files --------------------------------------

  # OUTPUT DIR
  e$outdir <- outdir

  # OUTPUT TEXT FILE
  e$log_file <- file( paste0(e$outdir,"/output.txt"), "w" )
  e$sink_file <- file( paste0(e$outdir,"/stats.txt"), "w" )

  sink( file = e$sink_file, append = TRUE )

  cat("LAND ICE ANALYSIS\n\n", file = e$log_file)

  cat("Experiment:", expt, "\n", file = e$log_file)
  if (expt == "SA") cat("Melt samples for ice sheet T-dep plots: N = ", N_melt_Tdep, "\n", file = e$log_file)

  if (no_emulator_uncertainty_mc) cat("No emulator uncertainty sampling in MC projections!\n", file = e$log_file)

  # OUTPUT CSV FILES
  csv_mme <- list()
  csv_mme[["RCP26"]] <- paste0( e$outdir,"/mme_RCP26.csv")
  csv_mme[["RCP45"]] <- paste0( e$outdir,"/mme_RCP45.csv") # This is for glaciers only
  csv_mme[["RCP85"]] <- paste0( e$outdir, "/mme_RCP85.csv")
  for (scenario in names(csv_mme)) cat("ice_source,region,year,mean,sd\n", file = csv_mme[[scenario]])

  if ( expt != "SA" && expt != "sim_only" ) {

    csv_full <- list()
    csv_summary <- list()

    # CSV FILES FOR EACH SSP PROJECTION
    for (scen in scenario_list[[temp_prior]]) {

      csv_summary[[scen]] <- paste0( e$outdir, "/summary_", temp_prior, "_", scen, ".csv")
      cat( paste("ice_source,region,year,", paste( paste0("q",q_list), collapse = ","),
                 ",sample_mean,sample_sd,sample_min, sample_max\n", sep = "" ),
           file = csv_summary[[scen]] )

      csv_full[[scen]] <- paste0( e$outdir, "/projections_", temp_prior, "_", scen, ".csv")
      cat( "ice_source,region,year,sample,GSAT,melt,collapse,SLE\n", file = csv_full[[scen]] )

    }

  }

  # region names --------------------------------------

  r_nums <- 1:19 # Glacier regions

  # Match region names in SLE CSV
  e$region_list <- list()
  e$region_list[["GrIS"]] <- "ALL"
  e$region_list[["AIS"]] <- c("WAIS", "EAIS", "PEN")
  e$region_list[["Glaciers"]] <- paste("region", r_nums, sep = "_")

  # For plots
  e$region_name_list <- list()
  e$region_name_list[["GrIS"]] <- "Greenland"
  e$region_name_list[["AIS"]] <- c("West Antarctica", "East Antarctica", "Antarctic Peninsula")
  regionnames.file <- system.file("extdata", "regionnames.txt", package = e$packagename, mustWork = TRUE)
  e$region_name_list[["Glaciers"]] <- as.character(unlist(read.csv(regionnames.file, header = FALSE)))

  # Add 'peripherals' for ice sheet glaciers to avoid ambiguity
  e$region_name_list[["Glaciers"]][5] <- paste(e$region_name_list[["Glaciers"]][5], "periphery")
  e$region_name_list[["Glaciers"]][19] <- paste(e$region_name_list[["Glaciers"]][19], "periphery")

  #____________________________________________________________

  # MELT VALUES: kappa (GrIS) and gamma0 (AIS)
  # Values for vertical lines in read_melt()
  # and med used for using/finding default in sensitivity_analysis()
  e$melt_values <- list()
  e$melt_values[["GrIS"]] <- c(-0.37, -0.17, -0.06, -0.9705, 0.0079 )
  names(e$melt_values[["GrIS"]]) <- c("high", "med", "low", "pc05", "pc95") # Note names refer to response, not value
  e$sc <- 1000 # Scaling factor for neater plots
  e$melt_values[["AIS"]] <- c(9619, 14477, 21005, 86984, 159188, 471264)/e$sc
  names(e$melt_values[["AIS"]]) <- c("low", "med", "high", "PIGL_low", "PIGL_med", "PIGL_high")

  # output files --------------------------------------

  # OUTPUT DIR
  e$outdir <- outdir

  # OUTPUT TEXT FILE
  e$log_file <- file( paste0(e$outdir,"/output.txt"), "w" )
  e$sink_file <- file( paste0(e$outdir,"/stats.txt"), "w" )

  sink( file = e$sink_file, append = TRUE )

  cat("LAND ICE ANALYSIS\n\n", file = e$log_file)

  cat("Experiment:", expt, "\n", file = e$log_file)
  if (expt == "SA") cat("Melt samples for ice sheet T-dep plots: N = ", N_melt_Tdep, "\n", file = e$log_file)

  # SLE values at which density estimates are made
  # for T-dependent plots and projections
  sle_lim <- seq(from = -50, to = 100, by = 0.001)
  e$smid <- sle_lim[-length(sle_lim)] + diff(sle_lim) / 2

  # Store emulator for each region (only most recent year)
  e$input_centre <- list()
  e$input_scale <- list()
  e$emulator <- list()

  # Max temperature of plots (depends on year): used in read_forcing()
  e$max_temp <- list()
  e$max_temp[["2025"]] <- 6
  e$max_temp[["2050"]] <- 6
  e$max_temp[["2075"]] <- 7.5
  e$max_temp[["2100"]] <- 7.5

  # Melt parameter limits: used in read_melt() and sensitivity_analysis()
  e$range_melt <- list()
  e$range_melt[["GrIS"]] <- c(-2, 0.5)
  e$range_melt[["AIS"]] <- c(-50000, 700000)/e$sc

  # SSP/SCP scenario colours - uses last two characters of scenario name
  e$scen_col <- list()
  e$scen_col[["19"]] <- rgb(30, 150, 132, maxColorValue = 255) # SSP1-19
  e$scen_col[["26"]] <- rgb(29, 51, 84, maxColorValue = 255) # SSP1-26
  e$scen_col[["45"]] <- rgb(234, 221, 61, maxColorValue = 255) # SSP2-45
  e$scen_col[["60"]] <- rgb(232, 136, 49, maxColorValue = 255) # SSP4-60
  e$scen_col[["70"]] <- rgb(242, 17, 17, maxColorValue = 255) # SSP3-70
  e$scen_col[["85"]] <- rgb(132, 11, 34, maxColorValue = 255) # SSP5-85
  e$scen_col[["DC"]] <- grey(0.2) # SSPNDC / NDC

  # Same but 60% transparency so can overlap
  e$scen_col_trans <- list()
  e$scen_col_trans[["19"]] <- rgb(30, 150, 132, maxColorValue = 255, alpha = 153) # SSP1-19
  e$scen_col_trans[["26"]] <- rgb(29, 51, 84, maxColorValue = 255, alpha = 153) # SSP1-26
  e$scen_col_trans[["45"]] <- rgb(234, 221, 61, maxColorValue = 255, alpha = 153) # SSP2-45
  e$scen_col_trans[["60"]] <- rgb(232, 136, 49, maxColorValue = 255, alpha = 153) # SSP4-60
  e$scen_col_trans[["70"]] <- rgb(242, 17, 17, maxColorValue = 255, alpha = 153) # SSP3-70
  e$scen_col_trans[["85"]] <- rgb(132, 11, 34, maxColorValue = 255, alpha = 153) # SSP5-85
  e$scen_col_trans[["DC"]] <- grey(0.2, alpha = 0.2) # SSPNDC / NDC

  # Same but 20% transparency
  e$scen_col_trans2 <- list()
  e$scen_col_trans2[["19"]] <- rgb(30, 150, 132, maxColorValue = 255, alpha = 51) # SSP1-19
  e$scen_col_trans2[["26"]] <- rgb(29, 51, 84, maxColorValue = 255, alpha = 51) # SSP1-26
  e$scen_col_trans2[["45"]] <- rgb(234, 221, 61, maxColorValue = 255, alpha = 51) # SSP2-45
  e$scen_col_trans2[["60"]] <- rgb(232, 136, 49, maxColorValue = 255, alpha = 51) # SSP4-60
  e$scen_col_trans2[["70"]] <- rgb(242, 17, 17, maxColorValue = 255, alpha = 51) # SSP3-70
  e$scen_col_trans2[["85"]] <- rgb(132, 11, 34, maxColorValue = 255, alpha = 51) # SSP5-85
  e$scen_col_trans2[["DC"]] <- grey(0.8, alpha = 0.2) # SSPNDC / NDC

  # read data --------------------------------------

  # Years to extract from forcing and sea level datasets
  # 2015 is for calculating anomalies
  # e$y_obs is for history matching / plotting with IMBIE
  # 2099 is to fill some 2100 gaps; 2100 is to print details when reading in
  # Global environment because used in multiple functions
  e$years_data <- sort(unique(c(2015, e$y_obs, 2099, 2100, e$years_pred)))
  stopifnot(min(e$years_data) >= 2015)

  # Change from numeric to dataset column format
  e$years_data <- paste0("y", e$years_data)
  e$years_pred <- paste0("y", e$years_pred)

  # READ CLIMATE PRIOR
  stopifnot(temp_prior %in% c( "FAIR", "CMIP6"))
  read_forcing(scenario_list, temp_prior, N_temp, climate_prior_kde, mean_temp, dataset)

  # READ SEA LEVEL PROJECTIONS
  # Read file and return SL dataset
  sle.filename <- ifelse(old_data == TRUE, "20191216_SLE_SIMULATIONS.csv", "20201106_SLE_SIMULATIONS.csv")
  e$sle_data <- read_sealevel( min_res, old_data, sle.filename )

  # Fixed_melt value for sensitivity analysis
  mean_melt_value <- list()
  for (is in c("GrIS", "AIS")) {
    # mean_melt_value[[is]] <- e$open_melt[[is]] # mean of ensemble - not that meaningful really
    mean_melt_value[[is]] <- e$melt_values[[is]][["med"]] # default value
  }

  # Grab col_names for parameters
  sle_col_names <- names(e$sle_data)[1:e$ncol_param]

  # READ MELT PRIORS
  read_melt(gamma0_prior)

  # emulators --------------------------------------

  cat("\n\nemulators ---------------------------\n", file = e$log_file)

  # EMULATOR PACKAGE: DiceKriging (old) or RobustGaSP (default)
  e$emul_type <- "RG"
  stopifnot(e$emul_type %in% c("DK", "RG"))

  # Default in RG is to bound; set to FALSE because warning when searching
  # "NA/Inf replaced by maximum positive value"
  # (stil searches even if not used!)
  # and to check for inert inputs
  if (e$emul_type == "RG") bound_corr_lengths = FALSE

  # Formula terms by ice source or AIS region CHANGE select() in loop too!
  # N.B. Not set up to do interactions
  mean_fn <- list()
  mean_fn[["ALL"]] <- "~ temp + melt"
  mean_fn[["WAIS"]] <- "~ temp + melt + collapse" # v1 data and DK: "~ temp + melt"
  mean_fn[["EAIS"]] <- "~ temp + melt + collapse"
  mean_fn[["PEN"]] <- "~ temp + melt + collapse" # v1 data and DK: "~ temp + collapse"

  # Separate for each region in case adding dummy model names (see below)
  for (reg in e$region_list[["Glaciers"]]) mean_fn[[reg]] <- "~ temp"

  # Option: add dummy variable for open vs standard ice sheet melt
  # xxx improve: only do this if "is" is in e$ice_source_list
  # otherwise confusing to add to output.txt
  if (e$add_dummy == "melt") {

    cat("\nAdd dummy variable for standard/open melt\n", file = e$log_file)

    for (is in "GrIS") { # Can add others in e$ice_source_list
      cat("\n", is, ":\n", file = e$log_file)
      for (reg in e$region_list[[is]] ) {
        cat("\n", reg, ":\n", file = e$log_file)
        mean_fn[[reg]] <- sprintf( "%s + %s",  mean_fn[[reg]], "melt0" )
        cat(mean_fn[[reg]], "\n", file = e$log_file)
      }
    }

  }

  # Option: add dummy variable for model or group (overfits though)
  if (e$add_dummy %in% c("model", "group")) {

    # These come from read_sealevel()
    cat("\nAdd unique model names for each ice sheet to formula\n", file = e$log_file)
    cat("(drop first model name to avoid dummy variable trap)\n", file = e$log_file)

    for (is in e$ice_source_list) {

      if (is == "Glaciers" && e$add_dummy == "group") next

      for (reg in e$region_list[[is]] ) {

        # If region included (redundant):
        if (length(e$ice_models[[reg]]) > 0) {

          cat("\n", reg, ":\n", file = e$log_file)

          # Drop first [-1] when adding to avoid dummy variable trap
          # If change to choosing one to exclude, need to change excluded column in priors too
          mean_fn[[reg]] <- sprintf("%s + %s",  mean_fn[[reg]],
                                    paste(paste(e$add_dummy, unlist(e$ice_models[[reg]])[-1], sep = "_"), collapse = " + " ))

          cat(mean_fn[[reg]], "\n", file = e$log_file)

        }
      }
    }

  }

  # Covariance functions by region - note over-ride option after for testing
  e$covariance_fn <- list()
  e$alpha_powexp <- list()

  # Greenland
  e$covariance_fn[["ALL"]] <- ifelse(e$emul_type == "RG", "pow_exp", "powexp")
  e$alpha_powexp[["ALL"]] <- 0.1 # power for pow_exp covariance in RobustGaSP

  # Antarctic regions
  e$covariance_fn[["WAIS"]] <- ifelse(e$emul_type == "RG", "pow_exp" , "powexp")
  e$alpha_powexp[["WAIS"]] <- 0.1

  e$covariance_fn[["EAIS"]] <- ifelse(e$emul_type == "RG", "pow_exp", "powexp")
  e$alpha_powexp[["EAIS"]] <- 0.1

  e$covariance_fn[["PEN"]] <- ifelse(e$emul_type == "RG", "pow_exp", "exp")
  e$alpha_powexp[["PEN"]] <- 0.1

  # Glacier regions
  e$covariance_fn[["region_1"]] <- ifelse(e$emul_type == "RG", "pow_exp", "powexp")
  e$alpha_powexp[["region_1"]] <- 1.0
  e$covariance_fn[["region_2"]] <- ifelse(e$emul_type == "RG", "pow_exp", "powexp")
  e$alpha_powexp[["region_2"]] <- 1.9
  e$covariance_fn[["region_3"]] <- ifelse(e$emul_type == "RG", "pow_exp", "matern3_2")
  e$alpha_powexp[["region_3"]] <- 1.9
  e$covariance_fn[["region_4"]] <- ifelse(e$emul_type == "RG", "pow_exp", "exp")
  e$alpha_powexp[["region_4"]] <- 0.1
  e$covariance_fn[["region_5"]] <- ifelse(e$emul_type == "RG", "pow_exp", "powexp")
  e$alpha_powexp[["region_5"]] <- 0.1
  e$covariance_fn[["region_6"]] <- ifelse(e$emul_type == "RG", "pow_exp", "matern5_2")
  e$alpha_powexp[["region_6"]] <- 1.0
  e$covariance_fn[["region_7"]] <- ifelse(e$emul_type == "RG", "pow_exp", "powexp")
  e$alpha_powexp[["region_7"]] <- 0.1
  e$covariance_fn[["region_8"]] <- ifelse(e$emul_type == "RG", "pow_exp", "matern5_2")
  e$alpha_powexp[["region_8"]] <- 1.0
  e$covariance_fn[["region_9"]] <- ifelse(e$emul_type == "RG", "pow_exp", "matern5_2")
  e$alpha_powexp[["region_9"]] <- 1.0
  e$covariance_fn[["region_10"]] <- ifelse(e$emul_type == "RG", "pow_exp", "powexp")
  e$alpha_powexp[["region_10"]] <- 1.0
  e$covariance_fn[["region_11"]] <- ifelse(e$emul_type == "RG", "pow_exp", "matern3_2")
  e$alpha_powexp[["region_11"]] <- 1.0
  e$covariance_fn[["region_12"]] <- ifelse(e$emul_type == "RG", "matern_3_2", "matern3_2")
  e$covariance_fn[["region_13"]] <- ifelse(e$emul_type == "RG", "pow_exp", "exp")
  e$alpha_powexp[["region_13"]] <- 0.1
  e$covariance_fn[["region_14"]] <- ifelse(e$emul_type == "RG", "pow_exp", "exp")
  e$alpha_powexp[["region_14"]] <- 1.9
  e$covariance_fn[["region_15"]] <- ifelse(e$emul_type == "RG", "pow_exp", "matern5_2")
  e$alpha_powexp[["region_15"]] <- 0.1
  e$covariance_fn[["region_16"]] <- ifelse(e$emul_type == "RG", "pow_exp", "exp")
  e$alpha_powexp[["region_16"]] <- 0.1
  e$covariance_fn[["region_17"]] <- ifelse(e$emul_type == "RG", "matern_5_2", "matern3_2")
  e$covariance_fn[["region_18"]] <- ifelse(e$emul_type == "RG", "pow_exp", "matern5_2")
  e$alpha_powexp[["region_18"]] <- 0.1
  e$covariance_fn[["region_19"]] <- ifelse(e$emul_type == "RG", "matern_5_2", "matern5_2")

  # Over-ride for testing all regions at once
  if ( ! is.na(do_covar_fn) ) {
    for ( cov_name in names(e$covariance_fn)) {
      e$covariance_fn[[cov_name]] <- do_covar_fn
      e$alpha_powexp[[cov_name]] <- do_covar_alpha
    }
  }

  # glacier caps --------------------------------------
  # Table 3 Farinott et al.: first element is total volume, second SLE (mm)
  # Calculate SLE from total volume if same method as Marzeion et al.
  # For now use SLE i.e. volume above flotation
  e$max_glaciers <- list()
  e$max_glaciers[["region_1"]] <- c(18.98, 43.3)
  e$max_glaciers[["region_2"]] <- c(1.06, 2.6)
  e$max_glaciers[["region_3"]] <- c(28.33, 64.8)
  e$max_glaciers[["region_4"]] <- c(8.61, 20.5)
  e$max_glaciers[["region_5"]] <- c(15.69, 33.6)
  e$max_glaciers[["region_6"]] <- c(3.77, 9.1)
  e$max_glaciers[["region_7"]] <- c(7.47, 17.3)
  e$max_glaciers[["region_8"]] <- c(0.30, 0.7)
  e$max_glaciers[["region_9"]] <- c(14.64, 32.0)
  e$max_glaciers[["region_10"]] <- c(0.14, 0.3)
  e$max_glaciers[["region_11"]] <- c(0.13, 0.3)
  e$max_glaciers[["region_12"]] <- c(0.06, 0.2)
  e$max_glaciers[["region_13"]] <- c(3.27, 7.9)
  e$max_glaciers[["region_14"]] <- c(2.87, 6.9)
  e$max_glaciers[["region_15"]] <- c(0.88, 2.1)
  e$max_glaciers[["region_16"]] <- c(0.10, 0.2)
  e$max_glaciers[["region_17"]] <- c(5.34, 12.8)
  e$max_glaciers[["region_18"]] <- c(0.07, 0.2)
  e$max_glaciers[["region_19"]] <- c(46.47, 69.4)

  # Convert SLE from mm to cm
  # Use this loop if using volume instead
  for (nn in names(e$max_glaciers)) {
    e$max_glaciers[[nn]][2] <- e$max_glaciers[[nn]][2]/10.0
  }

  # plotting details --------------------------------------

  # Max SLE of plots
  e$range_sle <- list()
  e$range_sle[["ALL"]] <- c(-12, 37) # 1.1-34.6; emulator lower
  e$range_sle[["WAIS"]] <- c(-10, 30) # -4.2 to 26.1; emulator
  if (old_data) e$range_sle[["WAIS"]] <- c(-15, 30)
  e$range_sle[["EAIS"]] <- c(-10, 20) # -6.9 to 15.3
  if (old_data) e$range_sle[["EAIS"]] <- c(-20, 20)
  e$range_sle[["PEN"]] <- c(-3, 4) # -1.4 to 3.4
  if (old_data) e$range_sle[["PEN"]] <- c(-3, 7)
  e$range_sle[["region_1"]] <- c(-1, 6) # Alaska **
  e$range_sle[["region_2"]] <- c(0, 0.4) # Western Canada and U.S.
  e$range_sle[["region_3"]] <- c(-1.5, 6.5) # Arctic Canada (North) **
  e$range_sle[["region_4"]] <- c(-0.7, 3.5) # Arctic Canada (South)
  e$range_sle[["region_5"]] <- c(-1.5, 4.5) # Greenland peripherals *
  e$range_sle[["region_6"]] <- c(-0.7, 1.2) # Iceland
  e$range_sle[["region_7"]] <- c(-0.7, 3) # Svalbard *
  e$range_sle[["region_8"]] <- c(-0.02, 0.1) # Scandinavia
  e$range_sle[["region_9"]] <- c(-0.7, 4.5) # Russian Arctic *
  e$range_sle[["region_10"]] <- c(-0.02, 0.08) # North Asia
  e$range_sle[["region_11"]] <- c(0.0, 0.06) # Central Europe x
  e$range_sle[["region_12"]] <- c(0.0, 0.026) # Caucasus x
  e$range_sle[["region_13"]] <- c(-0.2, 1.5) # Central Asia
  e$range_sle[["region_14"]] <- c(-0.2, 1.2) # South Asia (West)
  e$range_sle[["region_15"]] <- c(-0.1, 0.5) # South Asia (East)
  e$range_sle[["region_16"]] <- c(-0.02, 0.05) # Low Latitudes
  e$range_sle[["region_17"]] <- c(-0.5, 1.5) # Southern Andes
  e$range_sle[["region_18"]] <- c(-0.01, 0.03) # New Zealand
  e$range_sle[["region_19"]] <- c(-2, 7) # Antarctic peripherals **

  # RCP name/colour lookup
  e$rcp_name_list <- list()
  e$rcp_name_list[["26"]] <- "RCP2.6"
  e$rcp_name_list[["45"]] <- "RCP4.5"
  e$rcp_name_list[["60"]] <- "RCP6.0"
  e$rcp_name_list[["85"]] <- "RCP8.5"

  # Subfigure labels for LOO glaciers
  subfig_list <- paste(letters[5:23])

  # Store prior and predictions for year
  temp_sample <- list()
  e$pred_mean <- list()
  prior_df <- list()

  # mme plots --------------------------------------

  cat("\n\nmme plots ---------------------------\n", file = e$log_file)

  if (expt == "sim_only") {
    plot_mme(is = "GrIS", reg = "ALL")
    plot_mme(is = "AIS", reg = "WAIS")
    plot_mme(is = "AIS", reg = "EAIS")
    plot_mme(is = "AIS", reg= "PEN")
    # plot_mme(is = "AIS", reg = "ALL") # obsolete and may not work
  }

  # year loop --------------------------------------

  # One time or annually resampled
  melt_sample <- list()
  collapse_sample <- list()

  for (yy in e$years_pred) {

    yy_num <- substr(yy, 2, nchar(yy))

    cat("\n\n_______________________________________________________________\n", file = e$log_file)
    cat("Year:", yy_num, "\n", file = e$log_file)
    cat("_______________________________________________________________\n", file = e$log_file)

    # Get this year and add column for forcing
    sim_yy <- e$sle_data %>%
      select( c(sle_col_names, yy) ) %>%
      cbind(NA)
    colnames(sim_yy)[ dim(sim_yy)[2] ] <- "temp"

    # Get temperature anomaly for this year - by row for each experiment (GCM/scenario pair)
    # xxx improve: better row-wise code?
    for ( ss in 1:dim(sim_yy)[1] ) {
      print(ss)
      print(sim_yy[ss, "temp"] )
      print( select(filter(e$forcing_calib, GCM == sim_yy[ss, "GCM"] & scenario == sim_yy[ss, "scenario"]), yy ))
      sim_yy[ss, "temp"] <- e$forcing_calib %>%
        filter(GCM == sim_yy[ss, "GCM"] & scenario == sim_yy[ss, "scenario"]) %>%
        select(yy) %>%
        unlist()
      stopifnot( !is.na(sim_yy[ss,"temp"]) )
    }

    if ( expt != "sim_only" && expt != "SA") {

      # Get prior forcing
      for (scen in names(e$forcing_prior)) {
        cat("\nGet temperature prior for", scen, "\n", file = e$log_file)

        # Get for this year only
        if (length(e$years_pred) == 1) temp_sample[[scen]] <- e$forcing_prior[[scen]] # timeslice run
        else {
          cat("for year", yy_num, "\n", file = e$log_file)
          temp_sample[[scen]] <- unlist(select(e$forcing_prior[[scen]], yy)) # timeseries run
        }

        # Make sure don't have multiple years
        stopifnot( is.null(dim(temp_sample[[scen]])) )
        cat("N temperatures:", length(temp_sample[[scen]]), "\n", file = e$log_file)

      }
    }

    # FOR EACH ICE SOURCE
    for (is in e$ice_source_list) {

      cat("\n_____________________________________________________\n", file = e$log_file)
      cat("Ice source:", is, yy_num, "\n", file = e$log_file)
      cat("_____________________________________________________\n", file = e$log_file)

      # region loop --------------------------------------

      # FOR EACH REGION
      for (rr in 1:length(e$region_list[[is]])) {

        reg <- e$region_list[[is]][rr]
        reg_name <- e$region_name_list[[is]][rr]

        # Print to sink and log files:
        cat("\n______________________________________\n")
        cat("Region:", reg_name, yy_num, "\n\n")

        cat("\n______________________________________\n", file = e$log_file)
        cat("Region:", reg_name, yy_num, "\n", file = e$log_file)
        cat("______________________________________\n", file = e$log_file)

        # ALL SIMULATION DATA FOR REGION
        sim <- sim_yy %>%
          filter(ice_source == is & region == reg) %>%
          as.data.frame()

        # Output summaries of simulation data
        for (scen in names(csv_mme)) {

          # Merge equivalent SSPs with RCPs
          scen2 <- NA
          if (scen == "RCP26") scen2 <- "SSP126"
          if (scen == "RCP45") scen2 <- "SSP245" # glaciers only
          if (scen == "RCP85") scen2 <- "SSP585"

          sim_scen <- filter(sim, scenario %in% c(scen, scen2)) %>%
            select(yy) %>%
            unlist()

          if (length(sim_scen) > 0) {
            cat(sprintf("%s,%s,%s,%.3f,%.3f\n",
                        is, reg, yy_num, mean(sim_scen), sd(sim_scen)),
                file = csv_mme[[scen]], append = TRUE)
          }
        }

        # If not emulating, skip the rest
        if ( expt == "sim_only" ) next

        # Set melt and collapse to NA for writing to csv_full when not used
        # if one_sample_AIS, don't want to reset
        if ( yy == e$years_pred[1] || annual_resample ) {

          # Extra complications for AIS regional or whole ice sheet sampling
          # Assumes WAIS was run before EAIS, but it always should be
          if ( one_sample_AIS && ( reg == "EAIS" || reg == "PEN") ) {
            print(paste("Using WAIS parameter sample values for year", yy_num))
            melt_sample[[reg]] <- melt_sample[["WAIS"]]
            collapse_sample[[reg]] <- collapse_sample[["WAIS"]]
          } else {
            print(paste("Resetting parameter sample values for year", yy_num, "to NA"))
            melt_sample[[reg]] <- rep(NA, N_temp)
            collapse_sample[[reg]] <- rep(NA, N_temp)
          }
        }

        # SEA LEVEL RESPONSE FOR YEAR
        e$output <- select(sim, sle = yy)

        # Get emulator mean and covariance function for ice source / region
        # Not really necessary to rename (legacy code)
        #trend <- ifelse( is == "AIS", mean_fn[[reg]], mean_fn[[is]] )
        trend <- mean_fn[[reg]]
        kernel <- e$covariance_fn[[reg]]
        alpha_reg <- NA
        cat("\nMean function:", trend, "\n", file = e$log_file)
        cat("Covariance function:", kernel, "\n", file = e$log_file)

        if (e$emul_type == "RG" && kernel == "pow_exp") {
          alpha_reg <- e$alpha_powexp[[reg]]
          cat("alpha =", alpha_reg, "\n", file = e$log_file)
        }

        # PARAMETERS
        # Get variables from formula string (convoluted...)
        var_list <- strsplit(sub("~","", gsub(" ", "", trend)), split = "\\+")[[1]]
        cat("Variables:", paste(var_list), "\n")

        e$input <- select(sim, var_list)
        cat("Check number of columns selected:", dim(e$input)[2], "\n\n")

        # model selection --------------------------------------

        if (is %in% c("GrIS", "AIS") && do_model_comp ) {

          cat("Compare models with step-wise BIC\n")
          nsims <- length(unlist(e$output))
          cat("Number of simulations = ", nsims, "\n")
          fulldata <- as.data.frame( cbind(e$input, e$output))

          cat("Range of models:\n")
          min_formula <- sprintf("sle ~ 1")
          lm_formula <- sprintf("sle ~(%s)^2", paste(names(e$input), collapse ="+"))
          cat(min_formula, "\n")
          cat(lm_formula, "\n\n")

          # Initial linear model
          mylm <- lm( as.formula(lm_formula), fulldata )

          # k = log(n) is for BIC, which is more parsimonious than AIC (k = 2)
          mystep <- MASS::stepAIC( mylm,
                                   scope = list(lower = as.formula(min_formula),
                                                upper = as.formula(lm_formula)),
                                   direction = "both", k = log(nsims) )
          print(mystep$anova)
          cat( sprintf("\nStepwise BIC suggests %s model: %s\n\n", reg_name, deparse(mystep$call$formula)))
        }

        cat("\nInputs:", names(e$input),"\n\n", file = e$log_file)
        summarise(tibble(e$input), mean = mean(.data[["temp"]]), sd = sd(.data[["temp"]])) %>%
          unlist %>%
          cat("Unscaled: Mean and s.d. of temp:",.,"\n", file = e$log_file)

        # Scale inputs
        cat("\nCentre and scale inputs (so mean = 0, s.d. = 1)\n", file = e$log_file)
        e$input <- scale(e$input)

        # Store scaling to undo later for plots
        e$input_centre[[reg]] <- attr(e$input,"scaled:center")
        e$input_scale[[reg]] <- attr(e$input,"scaled:scale")
        cat("Centre factors for each input:", e$input_centre[[reg]], "\n", file = e$log_file)
        cat("Scale factors for each input:", e$input_scale[[reg]],"\n\n", file = e$log_file)

        # Scale outputs a matrix, so convert back to data frame
        e$input <- as.data.frame(e$input)

        summarise(tibble(e$input), mean = mean(.data[["temp"]]), sd = sd(.data[["temp"]])) %>%
          unlist %>%
          cat("Scaled: Mean and s.d. of temp:",.,"\n", file = e$log_file)

        # build emulator --------------------------------------

        # Build emulator
        if (e$emul_type == "DK") e$emulator[[reg]] <- DiceKriging::km( formula = as.formula(trend),
                                                                       design = e$input, response = e$output,
                                                                       covtype = kernel, nugget.estim = TRUE,
                                                                       nugget = var(e$output) )
        if (e$emul_type == "RG") {
          # RobustGasp with linear trends and nugget estimation


          input_mat <- as.matrix(e$input)
          output_mat <- as.matrix(e$output)
          trend.rgasp <- cbind(rep(1,dim(input_mat)[1]), input_mat)
          e$emulator[[reg]] <- RobustGaSP::rgasp(design = input_mat, response = output_mat,
                                                 alpha = rep(alpha_reg, dim(as.matrix(input_mat))[2]),
                                                 lower_bound = bound_corr_lengths,
                                                 trend = trend.rgasp, kernel_type = kernel, nugget.est = TRUE)

          show(e$emulator[[reg]])

          if ( ! bound_corr_lengths ) {
            cat("\nCorrelation lengths unbounded. Checking for weakly influential inputs using threshold = 0.1...\n")
            check_inert <- RobustGaSP::findInertInputs(e$emulator[[reg]])
          }

          if (yy == "y2100" && expt == "SA") { # To avoid over-plotting

            # loo validation --------------------------------------
            cat("\nLOO validation\n\n", file = e$log_file)

            # Default RobustGaSP plots
            pdf( file = paste0( e$outdir, "/RobustGaSP_LOO_", reg,".pdf"), width = 9, height = 8)
            RobustGaSP::plot(e$emulator[[reg]])
            dev.off()

            loo_valid <- RobustGaSP::leave_one_out_rgasp(e$emulator[[reg]])
            wrong <- output_mat > (loo_valid$mean + 2*loo_valid$sd) | output_mat < (loo_valid$mean - 2*loo_valid$sd)
            loo_err <- loo_valid$mean - output_mat
            loo_std_err <- loo_err / loo_valid$sd

            # Distance: for choosing covariance matrix
            euclid_dist <- sqrt( sum( ((loo_valid$mean - output_mat)/loo_valid$sd)^2 ) )
            cat("Euclidean distance", euclid_dist, "\n", file = e$log_file)

            col_dots <- rep("deepskyblue4", length(output_mat))
            col_wrong <- rgb(243, 122, 107, maxColorValue = 255)
            col_dots[wrong] <- col_wrong

            frac_right <- 1 - (length(col_dots[wrong])/length(col_dots))

            cat(sprintf("Number within emulator 95%% intervals: %.1f%%\n",
                        frac_right*100.0), file = e$log_file)
            cat("Mean of absolute errors (cm):", mean(abs(loo_err)), "\n",
                file = e$log_file)
            cat(sprintf("Range of absolute errors (cm): [%.1f, %.1f]",
                        min(loo_err), max(loo_err) ), "\n",
                file = e$log_file)
            cat("Mean of standardised errors:", mean(loo_std_err), "\n",
                file = e$log_file)
            cat(sprintf("Range of standardised errors: [%.1f, %.1f]",
                        min(loo_std_err), max(loo_std_err) ), "\n",
                file = e$log_file)

            # ED: LOO plots --------------------------------------

            # Emulated vs simulated
            pdf( file = paste0( e$outdir, "/ED_LOO_", reg,".pdf"),
                 width = 9, height = 4)
            par(mfrow=c(1,2),mar = c(4, 5, 1.5, 0.5))
            plot(output_mat, loo_valid$mean, pch = 19, xaxs = "i", yaxs = "i",
                 xlab = "Simulated (cm SLE)", ylab = "Emulated (cm SLE)",
                 xlim = e$range_sle[[reg]], ylim = e$range_sle[[reg]], col = col_dots,
                 cex = 1.1, cex.main = 1.5, cex.axis = 1.5, cex.lab = 1.5 )
            abline(0, 1, lwd = 0.5)

            arrows( output_mat, loo_valid$mean - 2*loo_valid$sd,
                    output_mat, loo_valid$mean + 2*loo_valid$sd,
                    angle = 90, lwd = 0.5, col = col_dots,
                    length = 0.025, code = 3)

            reg_name_short <- reg_name
            if (rr == 19) reg_name_short <- "Antarctic and Subantarctic"

            text(e$range_sle[[reg]][1],
                 e$range_sle[[reg]][1] + 0.94*(e$range_sle[[reg]][2] - e$range_sle[[reg]][1]), pos = 4, reg_name_short, cex = 1.6 )
            text(e$range_sle[[reg]][1],
                 e$range_sle[[reg]][1] + 0.84*(e$range_sle[[reg]][2] - e$range_sle[[reg]][1]), pos = 4, cex = 1.4,
                 sprintf("%.1f%%", frac_right*100.0, file = e$log_file))
            text(e$range_sle[[reg]][1],
                 e$range_sle[[reg]][1] + 0.76*(e$range_sle[[reg]][2] - e$range_sle[[reg]][1]), pos = 4, cex = 1.4,
                 sprintf("%.2g cm", mean(abs(loo_err))))

            # Standardised errors
            hist(loo_std_err, xlim = c(-6,6), breaks = seq(-14.4, 14.4, by = 0.2),
                 main = NULL, xlab = "Standardised errors",
                 col = "deepskyblue4", border = "deepskyblue4",
                 cex = 1.5, cex.main = 1.5, cex.axis = 1.5, cex.lab = 1.7)
            hist(loo_std_err[wrong], xlim = c(-6,6), breaks = seq(-14.4, 14.4, by = 0.2),
                 col = col_wrong, border = col_wrong, add = TRUE)
            abline( v=0 )
            abline( h=0 )

            myhist <- hist(loo_std_err, breaks = seq(-14.4, 14.4, by = 0.2), plot = FALSE )
            subfig_y <- max(myhist$counts)

            if (reg == "ALL") subfig <- "a"
            if (reg == "WAIS") subfig <- "b"
            if (reg == "EAIS") subfig <- "c"
            if (reg == "PEN") subfig <- "d"
            if (is == "Glaciers") subfig <- subfig_list[rr]

            text( 5, 0.95*subfig_y, pos = 4, font = 2, subfig, cex = 1.7)

            dev.off()
            par(mfrow=c(1,1))

          } # end of LOO validation

        } # end of RobustGaSP

        # sensitivity analysis --------------------------------------

        # SA: T-dep, melt-dep plots for certain years
        if ( expt == "SA" && yy_num %in% names(e$max_temp) ) {

          cat("\n\nSensitivity analysis:",yy_num,"\n\n", file = e$log_file)

          sensitivity_analysis(year = yy, ice_source = is, region = reg,
                               sim = sim, emulator = e$emulator[[reg]], N_melt = N_melt_Tdep)
        }


        # projections --------------------------------------

        if ( expt != "SA" ) {

          # SCENARIO LOOP
          for (scen in scenario_list[[temp_prior]]) {

            cat( "\n", scen, "\n", file = e$log_file)

            # CONSTRUCT DUMMY PRIORS IF USING
            # Open vs standard melt - can rule out glaciers, but could be either ice sheet
            if ( is %in% c("GrIS", "AIS") && e$add_dummy == "melt" ) {
              if ( abs(e$dummy_melt_pred - 0.5) < 0.01 ) dummy_melt <- sample(c(0,1), N_temp, replace = TRUE)
              else dummy_melt <- rep(e$dummy_melt_pred, N_temp)
              dummy_prior <- data.frame(melt0 = dummy_melt)
            }

            # List of model column names (-1 to drop first name)
            # This seemed easier than extracting models from var_list
            # but should change if selecting which model to drop instead of just the first
            if (e$add_dummy %in% c("model", "group")) {

              ice_models <- paste(e$add_dummy, unlist(e$ice_models[[reg]])[-1], sep = "_")

              # Randomly put a 1 in a model column, or none for the dropped model,
              # by sampling from vector of length (nd + 1) made up of a 1 and the rest zeros
              nd <- length(ice_models) # Number of dummy variables (n_model - 1)
              dummy_prior <- as.data.frame( t(replicate(N_temp, sample( c(1,rep(0,nd)), nd, replace = FALSE ))) )
              names(dummy_prior) <- ice_models
            }

            if (is == "Glaciers") {

              # Variables: temperature + glacier model names (no melt or collapse)
              prior_df[[scen]] <- data.frame(temp = temp_sample[[scen]])

              # Add dummy prior columns
              if (e$add_dummy %in% c("group", "model")) prior_df[[scen]] <- cbind( prior_df[[scen]], dummy_prior )

            } else {

              # If first time around, or resampling every year
              if ( (yy == e$years_pred[1] && scen == scenario_list[[temp_prior]][1]) || annual_resample) {

                # If one_sammple_AIS, don't resample but reuse
                if (! one_sample_AIS || (one_sample_AIS && (reg != "EAIS" && reg != "PEN" ) ) )  {

                  print(paste("Sampling parameters for year",yy_num))

                  # Melt sample: one per climate value for book-keeping
                  if (mean_melt) melt_sample[[reg]] <- rep(mean_melt_value[[is]], N_temp) # fixed value
                  else {
                    melt_sample[[reg]] <- sample( unlist(e$melt_prior[[is]]), N_temp, replace = TRUE ) # prior
                  }

                  # Collapse sample: mix of on and off, or else only on / only off
                  if (is == "AIS") {
                    if ( abs(collapse_prior - 0.5) < 0.01 ) collapse_sample[[reg]] <- sample(c(0,1), N_temp, replace = TRUE)
                    else collapse_sample[[reg]] <- rep(collapse_prior, N_temp)
                  }

                }
              }

              prior_df[[scen]] <- data.frame(temp = unlist(temp_sample[[scen]]),
                                             melt = melt_sample[[reg]],
                                             collapse = collapse_sample[[reg]])

              # Add dummy prior columns
              if (e$add_dummy %in% c("group", "model", "melt")) {
                prior_df[[scen]] <- cbind( prior_df[[scen]], dummy_prior )
              }

              # Drop any terms not in trend formula
              prior_df[[scen]] <- select( prior_df[[scen]], names(e$input) )

            } # ice sheets

            # Rescale priors using same factors as when building emulator
            prior_df_scaled <- as.data.frame( scale(prior_df[[scen]],
                                                    center = e$input_centre[[reg]],
                                                    scale = e$input_scale[[reg]] ))

            # PROJECTIONS
            if (e$emul_type == "DK") e$pred_mean[[scen]] <- DiceKriging::predict( e$emulator[[reg]],
                                                                                  newdata = prior_df_scaled,
                                                                                  type = "UK", checkNames = TRUE )

            if (e$emul_type == "RG") {
              prior_df_scaled <- as.matrix(prior_df_scaled)
              trend.test.rgasp <- cbind(rep(1,dim(prior_df_scaled)[1]), prior_df_scaled)
              e$pred_mean[[scen]] <- RobustGaSP::predict(e$emulator[[reg]], prior_df_scaled, testing_trend = trend.test.rgasp)
            }

            # Cap glaciers mean prediction (don't alter uncertainty)
            if (is == "Glaciers" &&
                max(e$pred_mean[[scen]]$mean) > e$max_glaciers[[reg]][2]) {
              cat("Capping glacier", reg, "emulator prediction at", e$max_glaciers[[reg]][2], " cm\n", file = e$log_file)
              cat("Initial range:", range(e$pred_mean[[scen]]$mean), "\n", file = e$log_file)
              e$pred_mean[[scen]]$mean[ e$pred_mean[[scen]]$mean > e$max_glaciers[[reg]][2] ] <- e$max_glaciers[[reg]][2]
              cat("Final range:", range(e$pred_mean[[scen]]$mean), "\n", file = e$log_file)
            }

            # Plot Gaussian/student-t residuals as a function of the upper and lower bounds
            if (scen == "SSP585" && yy == "y2100") {

              resid_upper <- (e$pred_mean[[scen]]$mean + 2 * e$pred_mean[[scen]]$sd) - e$pred_mean[[scen]]$upper95
              resid_lower <- (e$pred_mean[[scen]]$mean - 2 * e$pred_mean[[scen]]$sd) - e$pred_mean[[scen]]$lower95

              plot( e$pred_mean[[scen]]$upper95, resid_upper,
                    ylab = "Tail residual (cm SLE)", xlab = "Upper/lower 95% interval (cm SLE)", col = "darkred",
                    ylim = c(-1,1), xlim = range( c(e$pred_mean[[scen]]$upper95, e$pred_mean[[scen]]$lower95) ),
                    type = "p", pch = 20, main = paste("Gaussian - student-t residuals:", reg, scen ))
              points( e$pred_mean[[scen]]$lower95, resid_lower, pch = 16, col = "darkblue" )
              abline( h = 0 )

              print("", quote = FALSE)
              print("Range upper and lower residuals (Gaussian minus student-t):")
              print( range( c(resid_upper, resid_lower) ) )
            }

            # estimate pdfs --------------------------------------

            # Add emulator uncertainty (e$pred_mean[[ scen ]]$sd)
            # to emulator mean predictions (e$pred_mean[[ scen ]]$mean) in two ways

            # 1. PDF FROM MONTE CARLO SAMPLE

            # Will save for each region & scenario so can sum after
            proj_tag <- paste(reg, scen, sep = "_")

            # Sample once from emulator uncertainty (Gaussian) for each prediction
            # (unless don't)
            if (no_emulator_uncertainty_mc) {
              e$pred_mc[[proj_tag]] <- e$pred_mean[[ scen ]]$mean
            }
            else {
              e$pred_mc[[proj_tag]] <- rnorm( length(e$pred_mean[[ scen ]]$mean),
                                              mean = e$pred_mean[[ scen ]]$mean ,
                                              sd = e$pred_mean[[ scen ]]$sd )
            }

            # Cap glaciers again
            if (is == "Glaciers" &&
                max(e$pred_mc[[proj_tag]]) > e$max_glaciers[[reg]][2]) {
              cat("Capping glacier", reg, "sample prediction at", e$max_glaciers[[reg]][2], " cm\n", file = e$log_file)
              cat("Initial range:", range(e$pred_mc[[proj_tag]]), "\n", file = e$log_file)
              e$pred_mc[[proj_tag]][e$pred_mc[[proj_tag]] > e$max_glaciers[[reg]][2] ] <- e$max_glaciers[[reg]][2]
              cat("Final range:", range(e$pred_mc[[proj_tag]]), "\n", file = e$log_file)
            }

            # 2. PDF BY DETERMINISTIC NUMERICAL INTEGRATION
            # Used for estimating regional quantiles and neat pdf plots

            # Holder for density estimate at each sea level value
            pred_int <- rep(0, length(e$smid))

            # For each emulator prediction, add density of Gaussian dist with emulator mean and sd
            for (pp in 1:length(e$pred_mean[[ scen ]]$mean)) {

              # Add density for this prediction (normal dist with emulator mean and s.d.)
              pred_int <- pred_int + dnorm(e$smid, mean = e$pred_mean[[ scen ]]$mean[pp],
                                           sd = e$pred_mean[[ scen ]]$sd[pp])
            }

            # Truncate glacier pdfs at maximum mass loss before normalising
            if (is == "Glaciers") {
              cat("Capping glacier", reg, "pdf estimate at", e$max_glaciers[[reg]][2], " cm\n", file = e$log_file)
              pred_int[ e$smid > e$max_glaciers[[reg]][2]] <- 0.0
            }

            # Density estimate for this region -> for plots, quantiles
            # Normalise area
            e$pred_dens[[proj_tag]] <- pred_int / sum(pred_int)

            if (yy == "y2100") {

              # Short name of scenario for selecting colour
              sc <- substr(scen, nchar(scen)-1, nchar(scen))

              # Histogram bin width for MME data
              if (is %in% c("GrIS", "AIS")) {
                hist_inc <- 0.5
                if (reg == "PEN") hist_inc <- 0.2
              } else {
                if (reg %in% paste("region", c(2, 8,10:12,16,18), sep = "_")) hist_inc <- 0.002
                else hist_inc <- 0.1
              }

              # hists & pdfs --------------------------------------

              # If MME data exist for this scenario, plot histogram and emulator projections together
              scen2 <- NA
              if (scen == "SSP126") scen2 <- "RCP26"
              if (scen == "SSP245") scen2 <- "RCP45"
              if (scen == "SSP585") scen2 <- "RCP85"

              # Same selection as used for mme_csv region at start of loop
              sim_scen <- filter(sim, scenario %in% c(scen, scen2)) %>%
                select(yy) %>%
                unlist()

              # If data found, plot it
              if (length(sim_scen) > 0) {

                this_file <- paste0("ED_hist_", scen)
                pdf( file = paste0( e$outdir, "/", this_file, "_", reg,".pdf"), width = 9, height = 8)

                cat("Plot pdfs with simulation histograms\n", file = e$log_file)
                cat("Range of data:", paste(range(sim_scen), collapse = " - "), "\n", file = e$log_file)

                h <- hist( sim_scen, breaks = seq(from = e$range_sle[[reg]][1]-5*hist_inc,
                                                  to = e$range_sle[[reg]][2]+5*hist_inc, by = hist_inc), plot = FALSE )
                ymax <- max( c(density(e$pred_mc[[ proj_tag ]])$y, h$density ) )

                hist( sim_scen, freq = FALSE,
                      breaks = seq(from = e$range_sle[[reg]][1]-5*hist_inc,
                                   to = e$range_sle[[reg]][2]+5*hist_inc, by = hist_inc),
                      col = e$scen_col_trans2[[sc]], border = e$scen_col_trans2[[sc]],
                      cex = 1.5, cex.main = 1.5, cex.axis = 1.5, cex.lab = 1.5,
                      xlab = paste("Sea level contribution at", yy_num, "(cm SLE)"),
                      xlim = e$range_sle[[reg]], ylim = c(0, 1.1*ymax),
                      main = NULL, xaxs = "i", yaxs = "i" )

                abline( h = 0 )
                abline( v = mean(sim_scen), col = e$scen_col_trans2[[sc]], lwd = 2.2 )

                # Plot mean predictions i.e. without emulator uncertainty
                lines(density(e$pred_mean[[ scen ]]$mean), col = e$scen_col[[sc]], lwd = 2.2, lty = 3)

                # Plot MC sample
                lines( density(e$pred_mc[[ proj_tag ]]), col = e$scen_col[[sc]], lwd = 2.2 )

                # Region
                text(e$range_sle[[reg]][1], ymax, reg_name, pos = 4, cex = 1.5 )
                sct <- ifelse(is == "Glaciers", e$rcp_name_list[[sc]],
                              paste0(e$rcp_name_list[[sc]],"/",
                                     e$scen_name_list[[temp_prior]][ scenario_list[[temp_prior]] == scen ] ))

                text( e$range_sle[[reg]][2] - 0.4*(e$range_sle[[reg]][2]-e$range_sle[[reg]][1]), ymax,
                      paste("Simulations:", sct), pos = 4, cex = 1.5 )
                rect(e$range_sle[[reg]][2] - 0.45*(e$range_sle[[reg]][2]-e$range_sle[[reg]][1]), 1.02*ymax,
                     e$range_sle[[reg]][2] - 0.4*(e$range_sle[[reg]][2]-e$range_sle[[reg]][1]), 0.99*ymax,
                     col = e$scen_col_trans2[[sc]], border = e$scen_col_trans2[[sc]])

                text( e$range_sle[[reg]][2] - 0.4*(e$range_sle[[reg]][2]-e$range_sle[[reg]][1]),
                      0.93*ymax, "(mean)",
                      pos = 4, cex = 1.5 )
                lines(e$range_sle[[reg]][2] - c(0.45,0.4)*(e$range_sle[[reg]][2]-e$range_sle[[reg]][1]),
                      rep(0.93*ymax, 2), col = e$scen_col_trans2[[sc]], lwd = 2.2)

                text( e$range_sle[[reg]][2] - 0.4*(e$range_sle[[reg]][2]-e$range_sle[[reg]][1]),
                      0.86*ymax, paste("Emulator:", e$scen_name_list[[temp_prior]][ scenario_list[[temp_prior]] == scen ]),
                      pos = 4, cex = 1.5 )
                lines(e$range_sle[[reg]][2] - c(0.45,0.4)*(e$range_sle[[reg]][2]-e$range_sle[[reg]][1]),
                      rep(0.86*ymax, 2), col = e$scen_col_trans[[sc]], lwd = 2.2)

                text( e$range_sle[[reg]][2] - 0.4*(e$range_sle[[reg]][2]-e$range_sle[[reg]][1]),
                      0.79*ymax, "(no emulator uncertainty)",
                      pos = 4, cex = 1.5 )
                lines(e$range_sle[[reg]][2] - c(0.45,0.4)*(e$range_sle[[reg]][2]-e$range_sle[[reg]][1]),
                      rep(0.79*ymax,2), col = e$scen_col[[sc]], lwd = 2.2, lty = 3)

                dev.off()
              }

            } # Only year 2100 for hist pdf plots

            # write projections --------------------------------------

            for (tt in 1:N_temp) {
              cat( sprintf("%s,%s,%s,%i,%.4f,%.4f,%i,%.4f\n", is, reg, yy_num, tt,
                           unlist(temp_sample[[scen]][tt]),
                           unlist(melt_sample[[reg]][tt]), unlist(collapse_sample[[reg]][tt]),
                           unlist(e$pred_mc[[ proj_tag ]][tt])),
                   file = csv_full[[scen]], append = TRUE )
            }

            # SUMMARY FILE ROWS FOR REGION
            # Quantiles are from 'molehill' integration density estimates
            # Mean, sd, min, max from Monte Carlo sample

            # Interpolation of 95th percentile for some small capped glaciers was wrong
            # - discard repeated 0s and 1s at bottom/top of cdf to fix
            pred_cdf <- cumsum(e$pred_dens[[ proj_tag ]])
            pred_cdf_ind <- which(pred_cdf > 1e-9 & pred_cdf < (1 - 1e-9), arr.ind = TRUE)
            pred_cdf_ind <- c(min(pred_cdf_ind)-1, pred_cdf_ind, max(pred_cdf_ind) + 1)

            qq <- c(pred_cdf[pred_cdf_ind]) %>%
              approx(e$smid[pred_cdf_ind], xout = q_list)

            cat( sprintf("%s,%s,%s,%s,%.4f,%.4f,%.4f,%.4f\n",
                         is, reg, yy_num, paste(qq$y, collapse = ","),
                         mean(e$pred_mc[[ proj_tag ]]), sd(e$pred_mc[[ proj_tag ]]),
                         min(e$pred_mc[[ proj_tag]]), max(e$pred_mc[[ proj_tag]])),
                 file = csv_summary[[scen]], append = TRUE )


          } # SCENARIO LOOP

          if (yy %in% c("y2050","y2100")) {

            # MAIN: SSP pdfs --------------------------------------

            pdf( file = paste0( e$outdir, "/MAIN_SSP_pdfs_", reg, "_", yy_num, ".pdf"), width = 9, height = 8 )
            par(mar = c(5, 5, 1.5, 1))

            # Find max density across scenarios to set axis height
            ymax <- 0
            for (scen in scenario_list[[temp_prior]]) {
              proj_tag <- paste(reg, scen, sep = "_")
              ymax <- max(c(ymax, max(e$pred_dens[[ proj_tag ]])))
            }

            plot(1:3, 1:3, type = "n", main = NULL, xaxs = "i", yaxs = "i",
                 xlim = e$range_sle[[reg]], ylim = c(0, 1.2*ymax),
                 cex = 1.5, cex.main = 1.5, cex.axis = 1.8, cex.lab = 2,
                 xlab = paste("Sea level contribution at", yy_num, "(cm SLE)"), ylab = "Density")
            abline( v=0, col = "grey", lwd = 0.5 )

            ytext <- 1.1*ymax
            text( e$range_sle[[reg]][1] + 0.02*(e$range_sle[[reg]][2]-e$range_sle[[reg]][1]),
                  ytext, reg_name, pos = 4, cex = 2.5 )

            subfig <- ""
            if (reg == "ALL") subfig <- "e"
            if (reg == "region_3") subfig <- "f"
            if (reg == "WAIS") subfig <- "h"
            if (reg == "EAIS") subfig <- "i"

            text( e$range_sle[[reg]][2] - 0.1*(e$range_sle[[reg]][2]-e$range_sle[[reg]][1]),
                  ytext, pos = 4, font = 2, subfig, cex = 2.7)

            ytext <- 0.55*ymax
            for (scen in scenario_list[[temp_prior]]) {

              proj_tag <- paste(reg, scen, sep = "_")
              sc <- substr(scen, nchar(scen)-1, nchar(scen))
              scen_name <- e$scen_name_list[[temp_prior]][ scenario_list[[temp_prior]] == scen ]

              # Plot density vs values at which it was calculated, and add to legend
              lines( e$smid, e$pred_dens[[ proj_tag ]], col = e$scen_col[[sc]], lwd = 3 )
              text( e$range_sle[[reg]][2] - 0.25*(e$range_sle[[reg]][2]-e$range_sle[[reg]][1]),
                    ytext, scen_name, pos = 4, col = e$scen_col[[sc]], cex = 2.2 )
              ytext <- ytext - 0.07*ymax

            }

            # Add max glacier
            if (is == "Glaciers") abline(v = e$max_glaciers[[reg]][2], lty = 5, col = "darkgrey")



            dev.off()

          } # 2050 or 2100

        } # expt != SA, i.e. do projections

      } # END REGION LOOP

      # regional sums --------------------------------------

      # Just store Greenland projections MC sample in same format
      if (is == "GrIS") {

        if ( expt != "SA" && expt != "sim_only" ) {
          for (scen in scenario_list[[temp_prior]]) {
            ice_tag <- paste(is, scen, sep = "_")
            proj_tag <- paste("ALL", scen, sep = "_")
            e$pred_mc_is[[ice_tag]] <- e$pred_mc[[ proj_tag ]]
          }
        }

      } else {

        # SSPs that have MME data (RCP/SSP)
        if (is == "AIS") scen_list <- c("SSP126", "SSP585")
        if (is == "Glaciers") scen_list <- c("SSP126", "SSP245", "SSP585")

        # Add Antarctic and glacier regional projections together
        for (scen in scenario_list[[temp_prior]]) {

          ice_tag <- paste(is, scen, sep = "_")

          # Get MME data if it exists
          sim_reg <- NA
          if (scen %in% scen_list) {

            # Merge equivalent RCPs and SSPs
            if (scen == "SSP126") scen2 <- "RCP26"
            if (scen == "SSP245") scen2 <- "RCP45"
            if (scen == "SSP585") scen2 <- "RCP85"

            # Get simulation data for all regions
            sim_reg <- filter(sim_yy, ice_source == is & region %in% e$region_list[[is]]
                              & scenario %in% c(scen, scen2))

          }

          # Add all regions for ice source (do for all SSPs)
          mme_is <- NA
          mme_is_var <- NA
          for ( rr in 1:length(e$region_list[[is]]) ) {

            reg <- e$region_list[[is]][rr]
            proj_tag <- paste(reg, scen, sep = "_")

            # Sum MME regions for existing scenarios
            if (scen %in% scen_list) {

              # Sum over regions for each model
              if (is == "AIS") {
                # Get region simulation and add to total
                sr <- sim_reg %>%
                  filter(region == reg) %>%
                  select(yy) %>%
                  unlist()

                if (rr == 1) mme_is <- sr
                else mme_is <- mme_is + sr

              }

              # Sum means for each region (diff num of models per region)
              if (is == "Glaciers") {
                region_mean <- sim_reg %>%
                  filter(region == reg) %>%
                  select(yy) %>%
                  unlist() %>%
                  mean()
                region_var <- sim_reg %>%
                  filter(region == reg) %>%
                  select(yy) %>%
                  unlist() %>%
                  var()

                # xxx check: sum of regional variances definitely correct? (not used)
                mme_is <- ifelse(rr == 1, region_mean, mme_is + region_mean)
                mme_is_var <- ifelse(rr == 1, region_var, mme_is_var + region_var)

              }
            }

            # Sum emulator projections (all SSPs)
            if ( expt != "SA" && expt != "sim_only" ) {

              # Add predictions for each region, conditional on temperature
              # (Now only one per temperature anyway)
              if (rr == 1) e$pred_mc_is[[ice_tag]] <- e$pred_mc[[ proj_tag ]]
              else e$pred_mc_is[[ice_tag]] <- e$pred_mc_is[[ice_tag]] + e$pred_mc[[ proj_tag ]]
            }


          } # regions

          # Write MME data sums
          if (scen %in% scen_list) {

            # Get / calculate mean and s.d. and write to MME csv
            rm <- ifelse(is == "AIS", mean(mme_is), mme_is)
            rs <- ifelse(is == "AIS", sd(mme_is), sqrt(mme_is_var))

            cat( sprintf("%s,%s,%s,%.3f,%.3f\n", is,"ALL", yy_num, rm, rs),
                 file = csv_mme[[scen2]], append = TRUE)

          }

          # Write emulator projections for all ice sources
          if ( expt != "SA" && expt != "sim_only" ) {

            # Have to estimate quantiles from sample for summed regions (too hard to estimate densities)
            qq <- quantile(e$pred_mc_is[[ice_tag]], q_list)

            # Write to summary CSV: all not region name
            row_data <- sprintf("%s,%s,%s,%s,%.3f,%.3f,%.3f,%.3f\n", is,"ALL", yy_num,
                                paste(qq, collapse = ","),
                                mean(e$pred_mc_is[[ice_tag]]), sd(e$pred_mc_is[[ice_tag]]),
                                min(e$pred_mc_is[[ice_tag]]), max(e$pred_mc_is[[ice_tag]]))
            cat( row_data, file = csv_summary[[scen]], append = TRUE)
          }

        } # scenario loop for regional sum
      } # AIS and Glaciers

    } # END ICE SOURCE LOOP

    # land ice sum --------------------------------------
    if ( expt != "SA" && expt != "sim_only" ) {

      # Add all land ice emulator predictions for year
      for (scen in scenario_list[[temp_prior]]) {

        # Add regions
        for (is in e$ice_source_list) {

          ice_tag <- paste(is, scen, sep = "_")

          # Sum ice sheets and glaciers for each temperature
          if (is == e$ice_source_list[1]) e$land_ice <- e$pred_mc_is[[ice_tag]]
          else e$land_ice <- e$land_ice + e$pred_mc_is[[ice_tag]]

        }

        # Quantiles from summed samples
        qq <- quantile(e$land_ice, q_list)

        # Write to summary CSV
        cat( sprintf("%s,%s,%s,%s,%.3f,%.3f,%.3f,%.3f\n","LAND_ICE","ALL",
                     yy_num, paste(qq, collapse = ","),
                     mean(e$land_ice), sd(e$land_ice),
                     min(e$land_ice), max(e$land_ice)),
             file = csv_summary[[scen]], append = TRUE)

      }

      # table 1 --------------------------------------

      if (yy_num %in% c(2050, 2100) ) {

        cat("\n______________________________________\n", file = e$log_file)
        cat("\nTABLE 1:", yy_num,"\n", file = e$log_file)
        cat("\n______________________________________\n", file = e$log_file)

        # PRINT QUANTILES FOR TABLE 1
        for (is in c("LAND_ICE", e$ice_source_list)) {

          cat("\n", is, "\n", file = e$log_file)

          for (scen in scenario_list[[temp_prior]]) {

            # Emulated projections: percentiles
            pred_is <- read.csv(csv_summary[[scen]], header = TRUE) %>%
              filter(ice_source == is & region == "ALL" & year == yy_num)

            cat(sprintf("%s\t%.0f [%.0f, %.0f]\t[%.0f, %.0f]\n", scen,
                        select(pred_is, "q0.5"), select(pred_is, "q0.05"), select(pred_is, "q0.95"),
                        select(pred_is, "q0.17"), select(pred_is, "q0.83")), file = e$log_file)
          }
        }

      } # 2050 or 2100

    } # projections

  } # YEAR LOOP

  close(e$log_file)

  sink()

}
