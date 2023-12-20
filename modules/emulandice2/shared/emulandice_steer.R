#' ---
#' title: "emulandice2: predict"
#' output:
#'    html_notebook:
#'      toc: true
#'      number_sections: true
#' ---
#' MULTIVARIATE EMULATION OF LAND ICE CONTRIBUTIONS TO SEA LEVEL
# Tamsin Edwards
# based on emulandice (Edwards et al., 2021)
# and multivariate code by Jonathan Rougier
#
# Emulators are built by emulator_build.R using multi-model simulations of
# Greenland, Antarctica, and glacier regions (19) from EU H2020 PROTECT project.
# These simulations start between 1950 and 2005 and end between 2100-2300.
#
# Emulator projections here (main.R) usually begin in 1995, 2000 or 2005 ('cal_start')
# and end between 2100 and 2300 ('final_year').
#
#_______________________________________________________________________________

#' # Get FACTS args
# Get FACTS args ------------------------------------------------------------

print("Hello! Welcome to emulandice2: predict")

# Get arguments from RScript command in emulandice_steer.sh
args <- commandArgs(TRUE)
if (length(args) == 0) {
  # Defaults if no args set (used for testing and Markdown)
  warning("No arguments set - using defaults")
  i_s <- "GLA"
  reg <- "RGI03"
  emu_file <- "emu_file/GLA_RGI03_GloGEM_OGGM_pow_exp_20_EMULATOR.RData"
  climate_data_file <- "emulandice2.ssp585.temperature.fair.temperature_climate.nc"
  facts_ssp <- "ssp585"
} else {
 i_s <- args[1] # ice_source
 reg <- args[2] # region
 emu_file <- args[3] # emulator file
 climate_data_file <- args[4] # climate netcdf
 facts_ssp <- args[5] # ssp
}

print(paste("Ice source:", i_s))
stopifnot(i_s %in% c("GIS", "AIS", "GLA"))

# Region of ice source
print(paste("Region:", reg))
stopifnot(reg %in% c("ALL", paste0("RGI", sprintf("%02i",1:19))))

# Emulator: this is made from model_list and emulator_settings
print(paste("Emulator build:", emu_file))

# Netcdf name
print(paste("Climate file:", climate_data_file))

# SSP (could extract from filename)
scen <- paste0("SSP",substring(facts_ssp,4)) # emulandice expects upper case
print(paste("Scenario:", scen))

#emu_dir <- "./data-raw/"
#print(paste("Emulator build file directory:", emu_dir))

# LOAD EMULATOR AND OTHER STUFF
cat("Loading emulator build file\n")
load( file = emu_file)
inputs_ext <- "."
inputs_preprocess <- "."
outdir <- "./RESULTS/"

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# OTHER SETTINGS
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Write mean emulator projections as well as full projections with noise
# (used as arg for write_outputs(), and in plot_bayesian.R for some reason)
write_mean <- FALSE
print(paste("Write mean projections CSV:", write_mean))



#' # Setup
# Setup ------------------------------------------------------------

# EMULATOR BUILD FILE: constructed from the above settings
# Directory has to match rdatadir in the build file it is loading
#emu_file <- paste0(paste(i_s, reg, emu_name, sep = "_"), "_EMULATOR.RData")

#cat(sprintf("Looking for emulator build file in data/: %s\n", emu_file))
stopifnot(file.exists(emu_file))



#___________________________________________
# RESET SOME SETTINGS FOR FACTS

# Various code uses scenario loops
# so design_pred, myem, projections and projections_quant are set up as lists
# with label [[scen]]
# FACTS uses only one scenario at a time but will keep these loops for now:
# partly for plot scripts also used by emulator_build.R (multiple scenarios),
# partly in case I can use for plotting multiple scenarios later
scenario_list <- scen

# Plots: 0 = none, 1 = main, 2 = all
plot_level <- 0

# Number of 2LM projections of GSAT expected per SSP
# (and therefore total number of samples for book-keeping by GSAT value)
# Checks when reading in netcdf - could get rid of this
N_2LM <- 50L # 2237L for AR6 files


cat("Running...\n")

# Log file from emulator_build.R is same name but _build.txt
logfile_results <- paste0(outdir, out_name,"_results.txt")
cat(sprintf("\nemulandice2: %s %s\n\n", i_s, reg), file = logfile_results)

cat(sprintf("\nLoaded emulator file: %s\n", emu_file), file = logfile_results, append = TRUE)

#' # Create design
# Create design -----------------------------------------------------------------------
# FACTS: read in GSAT projections

# Future projections
# Design = "AR6_2LM" uses 2-layer model (FaIR) GSAT projections for SSPs
# i.e. climate_data_file is used in this function
# This overwrites uniform design_pred from RData build file (as intended, for reusing plot scripts)
design_pred <- load_design_to_pred( "AR6_2LM" )

# Store number of samples per scenario
# Used for outputs
N_temp <- length( design_pred[[1]][ , 1] )

#' # Predict: emulator mean projections
# Predict ----------------------------------------------------------------------
# emulator_predict() calls emu_mv with type = "var"
# FACTS: uses emulator object saved to RData workspace file

# myem <- list() # This is initalised in RData file with ME and uniform projections

# Rescale priors using same scaling factors as for simulator inputs
# Not the most elegant

cat("\nPredict:\n", file = logfile_results, append = TRUE)

for (scen in scenario_list) {

  cat(paste("Scenario:",scen,"\n"), file = logfile_results, append = TRUE)

  design_pred_scaled_cont <- scale(design_pred[[scen]][ , input_cont_list],
                                   center = inputs_centre,
                                   scale = inputs_scale )
  design_pred_scaled <- as.data.frame( design_pred[[scen]]  )
  design_pred_scaled[ , input_cont_list] <- design_pred_scaled_cont
  myem[[scen]] <- emulator_predict( design_pred_scaled )

}

#' ## Cap glacier mean projections
# GLACIER CAP: MEAN PROJECTIONS
for (scen in scenario_list) {
  if (i_s == "GLA" &&
      max( myem[[scen]]$mean ) > max_glaciers[[reg]] ) {
    cat( sprintf("\nCapping %s mean %s projections at %.3f cm SLE\n",
                 reg, scen, max_glaciers[[reg]]), file = logfile_results, append = TRUE)
    cat( sprintf("Initial range: %.3f - %.3f cm SLE\n", min( myem[[scen]]$mean ),
                 max( myem[[scen]]$mean )), file = logfile_results, append = TRUE)
    myem[[scen]]$mean[ myem[[scen]]$mean > max_glaciers[[reg]] ] <- max_glaciers[[reg]]
    cat( sprintf("Final range: %.3f - %.3f cm SLE\n", min( myem[[scen]]$mean ),
                 max( myem[[scen]]$mean)), file = logfile_results, append = TRUE)
  }
}

#' # Predict: emulator final projections (with uncertainty)
# Add emulator uncertainty to projections

# Projections list is initialised in build (unif_temps predictions)
# this overwrites the scenarios
for (scen in scenario_list) {

  # Generate final projections by sampling mean + uncertainty
  projections[[scen]] <- emulator_uncertainty(myem[[scen]])

  #' ## Cap glacier final projections

  # GLACIER CAP: FINAL PROJECTIONS
  if (i_s == "GLA" &&
      max( projections[[ scen ]] ) > max_glaciers[[reg]] ) {
    cat(sprintf("\nCapping %s final %s projections at %.3f cm SLE\n", reg, scen,
                max_glaciers[[reg]]), file = logfile_results, append = TRUE)

    cat( sprintf("Initial range: %.3f - %.3f cm SLE\n",
                 min( projections[[ scen ]]  ),
                 max( projections[[ scen ]] )), file = logfile_results, append = TRUE)

    projections[[ scen ]][ projections[[ scen ]] > max_glaciers[[reg]] ] <- max_glaciers[[reg]]
    cat( sprintf("Final range: %.3f - %.3f cm SLE\n", min( projections[[ scen ]]  ),
                 max( projections[[ scen ]]  )), file = logfile_results, append = TRUE)
  }

}

#' ## Calculate prior (uncalibrated) quantiles

projections_quant <- list()

for (scen in scenario_list) {

  # Projections
  projections_quant[[scen]] <- matrix( nrow = length(q_list), ncol = N_ts)
  colnames(projections_quant[[scen]]) <- paste0("y", years_em)
  rownames(projections_quant[[scen]]) <- paste0("q", q_list)
  for (yy in years_em) {
    projections_quant[[scen]][ , paste0("y", yy) ] <- quantile(projections[[scen]][, paste0("y", yy)], q_list)
  }

}

#' ## Output results

# PRINT SUMMARY
cat("\n_________________________________________\n", file = logfile_results, append = TRUE)
cat("PROJECTIONS: uncalibrated\n\n", file = logfile_results, append = TRUE)

for (yy in c(2100, 2300)) {

  if (yy %in% years_em) {

    cat(paste0("\n",yy,"\n"), file = logfile_results, append = TRUE)

    for (scen in scenario_list) {

      # Mean
      cat ("Mean:\n",  file = logfile_results, append = TRUE)
      cat(sprintf("%s: %.1f [%.1f,%.1f] cm SLE\n",
                  scen,
                  quantile(myem[[scen]]$mean[ , paste0("y", yy) ], probs = 0.5),
                  quantile(myem[[scen]]$mean[ , paste0("y", yy) ], probs = 0.05),
                  quantile(myem[[scen]]$mean[ , paste0("y", yy) ], probs = 0.95)),
          file = logfile_results, append = TRUE)

      # Final
      cat ("Final:\n",  file = logfile_results, append = TRUE)

      cat(sprintf("%s: %.1f [%.1f,%.1f] cm SLE\n",
                  scen,
                  projections_quant[[scen]][ q_list == 0.5, paste0("y", yy) ],
                  projections_quant[[scen]][ q_list == 0.05, paste0("y", yy) ],
                  projections_quant[[scen]][ q_list == 0.95, paste0("y", yy) ]),
          file = logfile_results, append = TRUE)
    }
  }
}

# PRINT UNCALIBRATED PROJECTIONS TO CSV AND NETCDF FILES
write_outputs(write_mean)

#' # Calibration
# Calibrate --------------------------------------------------------------------

#' ## Calculate model-obs differences

# xxx Make this multivariate! and rename because confusing
obs_change <- obs_data[obs_data$Year == cal_end,"SLE"] - obs_data[obs_data$Year == cal_start, "SLE"]
obs_change_err <- total_err[obs_data$Year == cal_end]

cat(paste0("\nObserved change (", cal_start,"-", cal_end, "): "), file = logfile_results, append = TRUE)
cat(sprintf("%.3f +/- %.3f cm SLE (1 sigma obs error)\n\n", obs_change, obs_data[obs_data$Year == cal_end,"SLE_sd"]), file = logfile_results, append = TRUE)
cat(sprintf("     +/- %.3f cm SLE (3 sigma total error)\n\n", 3*obs_change_err), file = logfile_results, append = TRUE)

# Calculate difference between each ensemble member (mean and final) and observations
dist_mean <- list()
dist_proj <- list()

cat("_________________________\n", file = logfile_results, append = TRUE)
cat("CALIBRATION\n", file = logfile_results, append = TRUE)

# xxx Note by not subtracting cal_start we are assuming already baselined to this year!
cat("\nCalculating difference between ensemble members and observations\n",
    file = logfile_results, append = TRUE)
for (scen in scenario_list) {
  dist_mean[[scen]] <- myem[[scen]]$mean[, paste0("y",cal_end) ] - obs_change
  dist_proj[[scen]] <- projections[[scen]][, paste0("y",cal_end) ] - obs_change
}

#' ## History matching

# History matching calibration for mean and final projections
cat("\nHistory matching calibration\n", file = logfile_results, append = TRUE)

# Save NROY projections
myem_nroy <- list()
proj_nroy <- list()

for (scen in scenario_list) {
  cat(paste0("\n", scen,"\n"), file = logfile_results, append = TRUE)
  cat("Calibrating mean projections:\n", file = logfile_results, append = TRUE)
  myem_nroy[[scen]] <- do_calibration(dist_mean[[scen]], "history_matching")
  cat("Calibrating full projections:\n", file = logfile_results, append = TRUE)
  proj_nroy[[scen]] <- do_calibration(dist_proj[[scen]], "history_matching")
}


#' ### Calculate posterior (calibrated) quantiles
projections_nroy_quant <- list()

for (scen in scenario_list) {

  # Calibrated projections
  projections_nroy_quant[[scen]] <- matrix( nrow = length(q_list), ncol = N_ts)
  colnames(projections_nroy_quant[[scen]]) <- paste0("y", years_em)
  rownames(projections_nroy_quant[[scen]]) <- paste0("q", q_list)
  for (yy in years_em) {
    projections_nroy_quant[[scen]][ , paste0("y", yy) ] <-
      quantile(projections[[scen]][proj_nroy[[scen]], paste0("y", yy)], q_list)
  }

}

#' ### Output results

cat("_______________________________________\n", file = logfile_results, append = TRUE)
cat("PROJECTIONS: calibrated (history matching)\n", file = logfile_results, append = TRUE)

for (yy in c(2100, 2300)) {

  if (yy %in% years_em) {

    cat(paste0("\n", yy,"\n"), file = logfile_results, append = TRUE)

    for (scen in scenario_list) {

      # Mean
      # cat(sprintf("NROY: %.1f [%.1f,%.1f] cm SLE",
      #                quantile(myem[[scen]]$mean[ myem_nroy[[scen]] , "y2100"], probs = 0.5),
      #                quantile(myem[[scen]]$mean[ myem_nroy[[scen]] , "y2100"], probs = 0.05),
      #                quantile(myem[[scen]]$mean[ myem_nroy[[scen]] , "y2100"], probs = 0.95)))

      # Final
      cat(sprintf("%s: %.1f [%.1f,%.1f] cm SLE\n",
                  scen,
                  projections_nroy_quant[[scen]][ q_list == 0.5, paste0("y",yy) ],
                  projections_nroy_quant[[scen]][ q_list == 0.05, paste0("y",yy) ],
                  projections_nroy_quant[[scen]][ q_list == 0.95, paste0("y",yy) ]),
          file = logfile_results, append = TRUE)
    }
  }
}

#' ## Bayesian calibration

cat("_______________________________________\n", file = logfile_results, append = TRUE)
cat("PROJECTIONS: calibrated (Bayesian)\n", file = logfile_results, append = TRUE)

cat("\nBayesian calibration\n", file = logfile_results, append = TRUE)

# Save normalised weights
myem_weights <- list()
proj_weights <- list()
for (scen in scenario_list) {
  myem_weights[[scen]] <- do_calibration(dist_mean[[scen]], "Bayesian")
  proj_weights[[scen]] <- do_calibration(dist_proj[[scen]], "Bayesian")
}

# PRINT CALIBRATED PROJECTIONS TO CSV FILES
# xxx add calibrated option to this function!
# write_outputs(write_mean)


# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# END OF FACTS ANALYSIS
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


if (plot_level > 0) {

  #' ### Plot results if requested

  cat("\nPlot uncalibrated projections:\n",file = logfile_results, append = TRUE)

  pdf( file = paste0( outdir, out_name, "_UNCALIBRATED.pdf"),
       width = 9, height = 5)
  plot_designs("prior", plot_level)
  plot_timeseries("prior", plot_level)
  plot_scatter("prior", "AR6_2LM", plot_level)
  plot_distributions("prior", plot_level)
  dev.off()

  cat("\nPlot calibrated projections:\n", file = logfile_results, append = TRUE)
  pdf( file = paste0( outdir, out_name, "_CALIBRATED.pdf"),
       width = 9, height = 5)
  plot_designs("posterior", plot_level)
  # Note no posterior time series plots
  plot_scatter("posterior", "AR6_2LM", plot_level)
  plot_distributions("posterior", plot_level)
  plot_bayesian()
  dev.off()
}

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Save workspace
save.image( paste0(outdir, out_name, "_RESULTS.RData") )

cat("...done.\n")

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


# Misc notes
# xxx Use obs_change and obs_change_err in plot_figures; implaus_thresh = 3 - done?

# w <- seq(1,1000)
# v <- sort(runif(1000))
# AR6_rgb_med[["SSP585"]] <- rgb(132, 11, 34, maxColorValue = 255, alpha = 153) # SSP5-85








