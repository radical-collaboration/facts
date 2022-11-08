# READ DATA FUNCTIONS
#
# Note:
# - Can still read old (2019) datasets to compare with submitted version
#   but will expect SICOPOLIS old names in model selection robustness tests
# - Cannot read intermediate datasets - pre-6th Nov for SLE_SIMULATIONS.csv
#   due to SICOPOLIS name changes
#_____________________________________________________________________
# Read in FORCING
#_____________________________________________________________________

read_forcing <- function(scenario_list, temp_prior, N_temp, climate_prior_kde, mean_temp, dataset) {

  #' Read forcing data from CSV
  #' @param scenario_list Scenarios to read
  #' @param temp_prior Name of GSAT prior ensemble
  #' @param N_temp Number of GSAT samples wanted
  #' @param climate_prior_kde Whether to smooth prior before sampling
  #' @param mean_temp Whether to use mean of GSAT prior instead of sampling
  #' @param dataset Which forcing CSV file to read: 2019, main, IPCC, FACTS
  #' @param temp_prior Which ensemble

  cat("\nread_forcing --------------------------------------\n", file = e$log_file)

  # read data --------------------------------------

  # File to read from package inst/extdata/ folder
  external_file <- FALSE
  if (dataset == "2019") forcing.filename <- "20191217_CLIMATE_FORCING.csv"
  else if (dataset == "main") forcing.filename <- "20201105_CLIMATE_FORCING.csv"
  else if (dataset == "IPCC") forcing.filename <- "20210215_CLIMATE_FORCING_IPCC.csv"
  else if (dataset == "FACTS") forcing.filename <- "FACTS_CLIMATE_FORCING.csv"
  else {
    forcing.filename <- dataset
    dataset <- "FACTS"
    external_file <- TRUE
  }

  # Number of initial columns before numeric data columns
  ncol_param <- 3

  cat("\nread_forcing: READ", forcing.filename, "\n\n", file = e$log_file)

  # READ CSV
  if (external_file) {
    if (file.exists(forcing.filename)) forcing.file <- forcing.filename
    else {
      stop("Forcing file not found")
    }
  }
  else forcing.file <- system.file( "extdata", forcing.filename, package = e$packagename, mustWork = TRUE )

  # tidyverse readr package: better defaults than read.csv; creates a tibble
  fd <- suppressMessages(read_csv( forcing.file ))

  # Add y to start of colname for tidyverse functions
  # Only needed for pre- 4th Oct 2020 datasets only
  if (dataset == "2019") colnames(fd)[ (ncol_param + 1):dim(fd)[2] ] <- paste0( "y", colnames(fd)[ (ncol_param + 1):dim(fd)[2] ] )

  cat("Read in", dim(fd)[1], "climate simulations\n", file = e$log_file)

  # write initial --------------------------------------

  cat("\nread_forcing: INITIAL DATA\n\n", file = e$log_file)

  cat("Global mean temperature (GSAT) change since 1850 (CMIP) or 1850-1900 mean (FAIR):\n\n", file = e$log_file)

  fd %>%
    select(ensemble, scenario, y2100) %>%
    group_by(ensemble, scenario) %>%
    summarise(n = n(),
              n_miss2100 = length( which(is.na(y2100)) ),
              q05 = quantile(y2100, prob = 0.05, na.rm = TRUE),
              median = median(y2100, na.rm = TRUE),
              q95 = quantile(y2100, prob = 0.95, na.rm = TRUE),
              mean = mean(y2100, na.rm = TRUE),
              sd = sd(y2100, na.rm = TRUE),
    ) %>%
    {.} -> climate_summ
  for (ii in 5:9) climate_summ[,ii] <- format(unlist(climate_summ[,ii]), digits = 3)
  suppressWarnings(write.table(climate_summ, file = e$log_file, append = TRUE, quote = FALSE, row.names = FALSE, sep = "\t"))


  # imputation --------------------------------------

  cat("\nread_forcing: DATA IMPUTATION\n\n", file = e$log_file)

  # REPLACE MISSING DATA AT 2100 WITH 2099
  fd %>%
    filter(is.na(y2100)) %>%
    summarise(
      n()
    ) %>%
    unlist %>%
    cat("Replacing",., "missing forcing values at 2100 with value at 2099\n", file = e$log_file)

  fd <- fd %>%
    mutate(y2100 = if_else(is.na(y2100), y2099, y2100))

  # data selection --------------------------------------

  cat("\n\nread_forcing: DATA SELECTION\n\n", file = e$log_file)

  # Drop any rows
  cat("Dropping rows that still have missing data\n", file = e$log_file)
  fd <- drop_na(fd)

  cat("Remaining:", dim(fd)[1], "climate simulations\n\n", file = e$log_file)

  cat("Selecting CMIP5+6 ensembles for calibration and", temp_prior, "ensemble for projections\n",
      file = e$log_file)

  # Keep CMIP5 and 6 (for calibrating emulator) and prior ensemble (for projections)
  if (temp_prior == "FAIR") {
    fd %>%
      filter( ensemble %in% c("CMIP5", "CMIP6", "FAIR"))
  }
  if (temp_prior == "CMIP6") {
    fd %>%
      filter( ensemble %in% c("CMIP5", "CMIP6"))
  }

  cat("Remaining:", dim(fd)[1], "climate simulations\n\n", file = e$log_file)

  # Select relevant scenarios too
  cat("Selecting scenarios (RCPs and SSPs)\n", file = e$log_file)
  fd <- fd[ fd$scenario %in% c(scenario_list[["CMIP5"]], scenario_list[["CMIP6"]],
                               scenario_list[[temp_prior]]), ]

  cat("Remaining:", dim(fd)[1], "climate simulations\n\n", file = e$log_file)

  # write final --------------------------------------

  cat("\nread_forcing: FINAL DATA\n\n", file = e$log_file)

  cat("Global mean temperature (GSAT) change since 1850 (CMIP) or 1850-1900 mean (FAIR):\n\n", file = e$log_file)

  fd %>%
    select(ensemble, scenario, y2100) %>%
    group_by(ensemble, scenario) %>%
    summarise(n = n(),
              n_miss2100 = length( which(is.na(y2100)) ),
              q05 = quantile(y2100, prob = 0.05, na.rm = TRUE),
              median = median(y2100, na.rm = TRUE),
              q95 = quantile(y2100, prob = 0.95, na.rm = TRUE),
              mean = mean(y2100, na.rm = TRUE),
              sd = sd(y2100, na.rm = TRUE)
    ) %>%
    {.} -> climate_summ
  for (ii in 5:9) climate_summ[,ii] <- format(unlist(climate_summ[,ii]), digits = 3)
  suppressWarnings(write.table(climate_summ, file = e$log_file, append = TRUE, row.names = FALSE, quote = FALSE, sep = "\t"))

  cat("\nReturning:", dim(fd)[1], "climate simulations\n\n", file = e$log_file)

  # calculate anomalies --------------------------------------

  # CALCULATE TEMPERATURE ANOMALIES WITH RESPECT TO 2015
  # Zero baseline at start of projections for temperature anomalies

  cat("\nread_forcing: ANOMALIES W.R.T. 2015\n\n", file = e$log_file)

  # Ditch the years we don't want
  fd <- select(fd, ensemble, GCM, scenario, e$years_data)

  # Calculate anomalies w.r.t first year
  fd %>%
    select(e$years_data) %>% # operate only on the numeric columns
    apply(1, function(x) as.numeric(x) - as.numeric(x[1])) %>% # otherwise char
    t() %>%
    {.} -> e$anom

  # Write calculated anomalies over original tibble (data frame) columns
  for (ii in 1:length(e$years_data))  fd[ , ncol_param + ii ] <- e$anom[,ii]

  cat("Global mean temperature (GSAT) change since 2015:\n\n", file = e$log_file)

  fd %>%
    select(ensemble, scenario, y2100) %>%
    group_by(ensemble, scenario) %>%
    summarise(n = n(),
              n_miss2100 = length( which(is.na(y2100)) ),
              q05 = quantile(y2100, prob = 0.05, na.rm = TRUE),
              median = median(y2100, na.rm = TRUE),
              q95 = quantile(y2100, prob = 0.95, na.rm = TRUE),
              mean = mean(y2100, na.rm = TRUE),
              sd = sd(y2100, na.rm = TRUE)
    ) %>%
    {.} -> climate_summ
  for (ii in 5:9) climate_summ[,ii] <- format(unlist(climate_summ[,ii]), digits = 3)
  suppressWarnings(write.table(climate_summ, file = e$log_file, append = TRUE, row.names = FALSE, quote = FALSE, sep = "\t"))

  # Construct prior here
  # i.e. select ensemble and scenario; resampling if needed

  # Time series to 2099: use FaIR as is
  # 2100: use resampling of FaIR or CMIP6 kde [check] -> plot these

  # FIGURE: 1b 2100 TEMPERATURE PRIOR

  # Return calibration forcing dataset
  e$forcing_calib <- filter(fd, ensemble %in% c("CMIP5", "CMIP6"))


  # calculate prior --------------------------------------

  cat("\nread_forcing: PRIOR RESAMPLING\n\n", file = e$log_file)

  cat("Number of temperature values per year requested:", N_temp, "\n", file = e$log_file)

  e$forcing_prior <- list()

  # For each scenario in our prior ensemble
  for (scen in scenario_list[[temp_prior]]) {

    # Get predictive years of prior ensemble for this scenario
    forcing_sc <- fd %>%
      filter(ensemble == temp_prior & scenario == scen) %>%
      select(e$years_pred)

    cat("Found", dim(forcing_sc)[1], "temp values for", scen, "\n", file = e$log_file)

    # If requested number exactly equals existing
    # (N_temp set to N_FAIR when predicting timeseries)
    if (N_temp == dim(forcing_sc)[1]) {

      if (length(e$years_pred) == 1) forcing_sc <- unlist(forcing_sc, use.names = FALSE)

      # Not resampling
      e$forcing_prior[[scen]] <- forcing_sc
      cat("Using prior data without resampling\n", file = e$log_file)

      # Do this check because assume N_temp values when predicting
      stopifnot(dim(e$forcing_prior[[scen]])[1] == N_temp)

      cat( "Returning temp prior of:", N_temp, "values for",
           dim(e$forcing_prior[[scen]])[2],"years\n", file = e$log_file )

    } else {

      # Check not trying to resample for more than one year
      stopifnot(length(e$years_pred) == 1)

      # Resample to get more values (i.e. for timeslice prediction)
      cat("Resampling",dim(forcing_sc)[1],"temp values to give",
          N_temp, "prior values for", scen, "\n", file = e$log_file)

      forcing_sc <- unlist(forcing_sc, use.names = FALSE)

      # Default for CMIP6 prior
      if (climate_prior_kde) {

        # Sample from density estimate: CMIP6 for 2100
        forcing_dens <- density(forcing_sc, n = 10000) # kde
        e$forcing_prior[[scen]] <- sample(forcing_dens$x, N_temp, replace = TRUE,
                                          prob = forcing_dens$y)
      } else {

        # Resample from original (i.e bootstrap): DEFAULT FOR 2100 (FAIR prior)
        stopifnot(dim(forcing_sc)[1] > 1) # Unlikely but check as sample does something unexpected
        e$forcing_prior[[scen]] <- sample(forcing_sc, N_temp, replace = TRUE) # direct
      }

      # Now take mean for fixed_T: repeat N_T times for multiple melt samples and book-keeping
      if (mean_temp == TRUE) {
        T_mean <- mean( e$forcing_prior[[scen]] )
        print(paste("Using mean temperature prior:", scen, T_mean))
        e$forcing_prior[[scen]] <- rep( T_mean, length = length( e$forcing_prior[[scen]] ) )
      }

      # Do this check because assume N_temp values when predicting
      stopifnot(length(e$forcing_prior[[scen]]) == N_temp)

      cat( "Returning temp prior of", N_temp, "values for", e$years_pred, "\n", file = e$log_file )

    }

  } # scenarios

  # MAIN: GSAT priors --------------------------------------

  # Plot climate priors if one year
  if (length(e$years_pred) == 1) {

    yy_num <- substr(e$years_pred, 2, nchar(e$years_pred))
    breaks <- seq(-1, 9, by = 0.2)

    pdf( file = paste0( e$outdir, "/MAIN_GSAT_", temp_prior, ".pdf" ), width = 9, height = 8)
    par(mar = c(6, 5, 1.5, 2.5))

    if (temp_prior == "FAIR") {
      ytext <- length(e$forcing_prior[[1]])*0.42
      yinc <- length(e$forcing_prior[[1]])*0.06
      subfig <- "a"
    }
    if (temp_prior == "CMIP6") {
      ytext <- 150
      yinc <- 23
      subfig <- "b"
    }

    for (scen in scenario_list[[temp_prior]]) {

      # Get colour tag
      sc <- substr(scen, nchar(scen)-1, nchar(scen))
      scen_name <- e$scen_name_list[[temp_prior]][ scenario_list[[temp_prior]] == scen ]

      if (scen == scenario_list[[temp_prior]][1]) {
        hist(e$forcing_prior[[scen]], freq = TRUE, breaks = breaks, xlim = c(-1,8),
             xlab = bquote("Global mean temperature change 2015-"*.(yy_num)~"("*degree*"C)"),
             main = NULL,
             col = e$scen_col_trans[[sc]], border = e$scen_col[[sc]],
             cex.axis = 1.6, cex.lab = 2 )
        text( 7.5, ytext, pos = 4, font = 2, subfig, cex = 2.5)
      } else {
        hist(e$forcing_prior[[scen]], freq = TRUE, breaks = breaks,
             col = e$scen_col_trans[[sc]], border = e$scen_col[[sc]], add = TRUE )
      }

      # Add to legend
      text( 6, ytext*0.5, scen_name, pos = 4, col = e$scen_col[[sc]], cex = 2.2 )
      ytext <- ytext - yinc

    }


    dev.off()
  } # years

} # read_forcing end of function


#_____________________________________________________________________
# Read in SEA LEVEL PROJECTIONS
#_____________________________________________________________________

read_sealevel <- function(min_res, old_data, sle.filename) {

  #' Read in sea level CSV
  #' @param min_res Arg: vector of minimum grid resolution in km for GrIS, AIS models
  #' @param old_data Whether data file is the old one (2019) for backward comparions (N.B. the data is wrong!)
  #' @param sle.filename Name of CSV to read in

  cat("\nread_sealevel --------------------------------------\n", file = e$log_file)

  # read data --------------------------------------

  # File to read
  cat("\nread_sealevel: READ", sle.filename, "\n\n", file = e$log_file)

  # READ CSV from package inst/extdata/ folder
  sle.file <- system.file( "extdata", sle.filename, package = e$packagename, mustWork = TRUE )
  fd <- suppressMessages(read_csv( sle.file ))

  # Add y to start of column names if old dataset
  if ( old_data == TRUE ) {

    # Initial number of columns before numeric data columns in this file
    ncol_param <- 10
    colnames(fd)[ (ncol_param + 1):dim(fd)[2] ] <- paste0( "y", colnames(fd)[ (ncol_param + 1):dim(fd)[2] ] )
  }

  # Get initial columns + projection years
  fd <- select(fd, ice_source, region, group, model, resolution, exp_id, GCM, scenario, melt, collapse, e$years_data)

  cat("Read in", dim(fd)[1], "sea level simulations\n", file = e$log_file)

  # write initial --------------------------------------

  cat("\nread_sealevel: INITIAL DATA\n\n", file = e$log_file)

  cat("Sea level equivalent (SLE) change since 2015:\n\n", file = e$log_file)
  cat("Note: Glacier summary doesn't account for model number diffs across regions\n", file = e$log_file)

  fd %>%
    select(ice_source, region, group, model, exp_id, scenario, y2100) %>%
    group_by(ice_source, scenario) %>%
    summarise(n = n(),
              n_miss2100 = length( which(is.na(y2100)) ),
              q05 = quantile(y2100, prob = 0.05, na.rm = TRUE),
              median = median(y2100, na.rm = TRUE),
              q95 = quantile(y2100, prob = 0.95, na.rm = TRUE),
              mean = mean(y2100, na.rm = TRUE),
              sd = sd(y2100, na.rm = TRUE)
    ) %>%
    {.} -> sealevel_summ
  for (ii in 5:9) sealevel_summ[,ii] <- format(unlist(sealevel_summ[,ii]), digits = 3)
  suppressWarnings(write.table(sealevel_summ, file = e$log_file, append = TRUE, quote = FALSE, row.names = FALSE, sep = "\t"))

  cat("\nNow", dim(fd)[1], "sea level simulations\n\n", file = e$log_file)

  # imputation 2100 --------------------------------------

  # REPLACING MISSING DATA AT 2100 WITH 2099
  # Only needed in old dataset due to date errors, or for BISICLES AIS
  # Moved up before plot
  fd %>%
    filter(is.na(y2100)) %>%
    summarise(
      n()
    ) %>%
    unlist %>%
    cat("Replacing",., "missing sea level values at 2100 with value at 2099\n", file = e$log_file)

  fd <- fd %>%
    mutate(y2100 = if_else(is.na(y2100), y2099, y2100))

  # data selection --------------------------------------
  # note this must come before imputation which uses mean over simulations to impute melt

  cat("\nread_sealevel: DATA SELECTION (IF ANY)\n\n", file = e$log_file)

  # OPTION: USE JUST ONE GCM PER ICE SHEET
  if (e$select_gcm) {
    for ( is in c("GrIS", "AIS") ) {
      if ( !is.na(e$cond_gcm[[is]][1]) ) {
        cat(is, ": use only models forced by", e$cond_gcm[[is]], "\n", file = e$log_file)
        fd <- filter(fd, ! (ice_source == is & ! GCM %in% e$cond_gcm[[is]]) )
        cat("Now", dim(fd)[1], "sea level simulations\n\n", file = e$log_file)
      }
    }
  }

  # OPTION: SELECT ICE SHEET MODELS
  if (e$select_ism) {

    cat("\nSelecting ice sheet / glacier models\n", file = e$log_file)

    for ( is in e$ice_source_list) {
      if ( !is.na(e$cond_ism[[is]][1]) ) {

        if (is == "Glaciers") {
          cat(is, ": selecting",length(cond_ism[[is]]), "models\n", file = e$log_file)

          # Remove any runs where ice_source is Glacier and model is not in glacier model list
          fd <- filter(fd, ! (ice_source == is & ! model %in% e$cond_ism[[is]]) )

        } else {

          # Split by "__" and check none failed
          group_model_split <- strsplit(e$cond_ism[[is]], split = "__")
          stopifnot( length(group_model_split) == length(e$cond_ism[[is]]) )

          cat(is, ": selecting",length(group_model_split), "group__model pairs\n", file = e$log_file)
          cat("List:", paste(e$cond_ism[[is]], sep = ", "), "\n", file = e$log_file)

          fd <- filter(fd, ! (ice_source == is &
                                ! ( group %in% unlist(lapply(group_model_split, `[[`, 1)) &
                                      model %in% unlist(lapply(group_model_split, `[[`, 2)) ) ) )
        }
      }
    }
    cat("Now", dim(fd)[1], "sea level simulations\n\n", file = e$log_file)
  }

  # Exclude open melt AIS runs
  if (e$exclude_open_AIS) {
    cat("Excluding open melt AIS experiments:\n\n", file = e$log_file)
    fd <- filter(fd, ! (ice_source == "AIS" & is.na(melt) ) )
    cat("Now", dim(fd)[1], "sea level simulations\n\n", file = e$log_file)
  }


  #_______________________________________________________________________________
  # OPTION: HISTORY MATCHING

  if (e$history_match) {
    cat("History matching: exclude AIS models with", e$y_obs,
        "contribution not within\n", #imbie_range,
        "IMBIE mean +/-", e$n_sig, "sigma:",
        e$imbie_obs["mean"], "+/-", e$n_sig * e$imbie_obs["sd"], "cm\n",
        "Range:", e$imbie_obs["mean"] + c(-1,1) * e$n_sig * e$imbie_obs["sd"], "cm\n\n",
        file = e$log_file)
  }

  # Column to select
  y_obs_char <- paste0("y", e$y_obs)

  plot(1:3, 1:3, type = "n",
       xlim = c(-0.5,0.7), ylim = c(-20, 50), xaxs = "i", yaxs="i", main = "RCP8.5 and SSP5-85", # "RCP26/SSP126"
       xlab = paste("Sea level contribution at", e$y_obs, "(cm SLE)"),
       ylab = paste("Sea level contribution at 2100 (cm SLE)"))
  abline(h = 0)
  abline(v = 0)
  rect( e$imbie_obs["mean"] - e$n_sig*e$imbie_obs["sd"], -20,
        e$imbie_obs["mean"] + e$n_sig*e$imbie_obs["sd"], 50, col = grey(0.8, alpha = 0.2), border = NA)
  abline(v = e$imbie_obs["mean"], col = "grey")
  gcm_list <- c("CESM2", "CNRM-CM6-1", "CNRM-ESM2-1", "NorESM1-M", "MIROC-ESM-CHEM", "CCSM4", "IPSL-CM5A-MR", "CSIRO-Mk3-6-0", "HadGEM2-ES")

  # Copy dataset
  fd_flag <- fd

  # For saving AIS sum
  e$sle_data_AIS <- fd

  for (ff in 1:dim(fd)[1]) {

    if (fd[ff,"region"] == "WAIS") {

      # Get unique run identifiers
      fd_g <- unlist(fd[ ff, "group"])
      fd_m <- unlist(fd[ff, "model"])
      fd_e <- unlist(fd[ff, "exp_id"])

      # Get all three parts off AIS
      fd_runs <- fd_flag %>%
        filter(ice_source == "AIS" &
                 group == fd_g & model == fd_m & exp_id == fd_e)

      # Calibration year SLE
      fd2 <- fd_runs %>%
        select( y_obs_char ) %>%
        unlist()
      stopifnot(length(fd2) == 3)

      # 2100 SLE for plot
      fd3 <- fd_runs %>%
        select("y2100") %>%
        unlist()
      stopifnot(length(fd3) == 3)

      e$sle_data_AIS[ ff, y_obs_char ] <- sum(fd2)
      e$sle_data_AIS[ ff, "y2100" ] <- sum(fd3)

      # HERE IS THE HISTORY MATCHING
      if (e$history_match) {
        # Add and set to NA if outside range
        if ( sum(fd2) < e$imbie_obs["mean"] - e$n_sig * e$imbie_obs["sd"]
             | sum(fd2) > e$imbie_obs["mean"] + e$n_sig * e$imbie_obs["sd"] ) {
          fd_flag[ fd_flag$ice_source == "AIS" & fd_flag$group == fd_g & fd_flag$model == fd_m & fd_flag$exp_id == fd_e, y_obs_char ] <- NA
        }
      }

      # Plot for high or low scenario
      fd_s <- unlist(fd[ff, "scenario"])
      if (fd_s %in% c("RCP85", "SSP585")) { # {c("RCP26", "SSP126")) {
        # c("rgb(180,221,212)", "rgb(28,94,57)", "rgb(194,223,125)", "rgb(96,64,155)", "rgb(198,192,254)", "rgb(51,84,122)", "rgb(63,211,74)", "rgb(197,8,158)", "rgb(28,224,178)"]
        fd_gcm <- unlist(fd[ff, "GCM"])
        points(sum(fd2), sum(fd3), col = which(gcm_list == fd_gcm), pch = which(gcm_list == fd_gcm) )
      }


    } # WAIS
  } # ff

  # Legend
  yleg <- 47
  for (gg in 1:length(gcm_list)) {
    text(-0.45, yleg, pos = 4, gcm_list[gg], col = gg )
    points(-0.45, yleg, pch = gg, col = gg )
    yleg <- yleg - 2
  }

  # region selection is because sum is stored in WAIS rows
  if (e$history_match) cat(e$y_obs, "range before:", range( select(filter(e$sle_data_AIS, region == "WAIS"), y_obs_char)), "\n", file = e$log_file)

  # Apply history matching to total - just because useful to check and store
  # Note this does nothing if history_match <- FALSE
  e$sle_data_AIS[ which(is.na(fd_flag[ , y_obs_char])), y_obs_char] <- NA
  e$sle_data_AIS <- drop_na(e$sle_data_AIS, y_obs_char)

  # Do history match i.e. drop NA set in loop above
  fd <- fd_flag %>% drop_na( y_obs_char )

  if (e$history_match) {
    cat(e$y_obs, "range after:", range( select(filter(e$sle_data_AIS, region == "WAIS"), y_obs_char)), "\n\n", file = e$log_file)
    cat("Now", dim(fd)[1], "sea level simulations\n\n", file = e$log_file)
  }

  # Get AIS totals for plot_mme.R / general use
  # Note different data.frame to e$sle_data because not post-processed e.g. melt not scaled, collapse = T/F
  e$sle_data_AIS <- filter(e$sle_data_AIS, ice_source == "AIS", region == "WAIS") %>%
    select(group, model, GCM, scenario, melt, collapse, y_obs_char, y2100)

  # End of scatter plot and history matching
  #_______________________________________________________________________________

  # OPTION: REMOVE MODELS WITH LOWER RESOLUTION
  for ( is in c("GrIS", "AIS") ) {

    if ( !is.na(min_res[[is]]) ) {
      cat(is, ": exclude models with resolution:", min_res[[is]], "km or greater\n", file = e$log_file)
      fd <- filter(fd, ! (ice_source == is & resolution >= min_res[[is]]) )
      cat("Now", dim(fd)[1], "sea level simulations\n\n", file = e$log_file)
    }
  }

  # GrIS 16km excludes 6 models (IMAUICE1)
  # AIS 32km excludes 132 models (fetish and two IMAUICE)

  cat("\nread_sealevel: DATA TRANSFORMATION\n\n", file = e$log_file)

  # AIS COLLAPSE: change from TRUE/FALSE to 0/1
  cat("Replacing AIS collapse T/F values with 0/1\n", file = e$log_file)

  # hist experiments have collapse = NA
  fd <- fd %>%
    mutate(collapse = case_when(ice_source == "AIS" & collapse == FALSE ~ 0,
                                ice_source == "AIS" & collapse == TRUE ~ 1 ) )

  # imputation --------------------------------------

  # IMPUTATION
  cat("\nread_sealevel: DATA IMPUTATION\n\n", file = e$log_file)

  # AIS: 'melt' is gamma0 parameter
  # Will set 'open' melt run gamma0 to mean of standard melt runs
  # Print mean of all non-missing gamma0 values (v1 open = 0; new open = NA)
  fd %>%
    filter(ice_source == "AIS" ) %>%
    select(melt) %>%
    summarise (
      mean(melt, na.rm = TRUE)
    ) %>%
    unlist() %>%
    cat("Initial mean of all AIS melt values:", ., "\n", file = e$log_file)

  # Standard melt runs, i.e.:
  # v1 data: open melt = 0
  # New data: open melt = NA
  if ( old_data == TRUE ) {
    fd2 <- fd %>%
      filter(ice_source == "AIS" & melt > 0.1 )
  } else {
    fd2 <- fd %>%
      filter(ice_source == "AIS" & !is.na(melt) )
  }

  # Save the number and mean of these std run gamma0 values
  fd2 %>%
    select(melt) %>%
    summarise (
      n = n(),
      mean = mean(melt, na.rm=FALSE)
    ) %>%
    unlist() %>%
    {.} -> std_melt_AIS

  cat(std_melt_AIS["n"], "AIS experiments with std melt\n", file = e$log_file)
  cat("Mean (imputed value):", std_melt_AIS["mean"], "\n", file = e$log_file)

  # Use alternative imputation: PIGL_med
  if (e$impute_high_melt_AIS) {
    cat("NOTE: replacing imputed mean value with alternative:\n", file = e$log_file)
    #std_melt_AIS["mean"] <- e$melt_values[["AIS"]][["PIGL_med"]]*e$sc # value is scaled
    std_melt_AIS["mean"] <- 150000 # value is scaled
    cat("New imputed value:", std_melt_AIS["mean"], "\n", file = e$log_file)
  }

  # Open melt experiments
  # i.e. the reverse of above
  if ( old_data == TRUE ) {
    fd2 <- fd %>%
      filter(ice_source == "AIS" & melt < 0.1 )
  } else {
    fd2 <- fd %>%
      filter(ice_source == "AIS" & is.na(melt) )
  }

  # Save the number (and, relevant for old_data, mean) of these open gamma0 values
  fd2 %>%
    select(melt) %>%
    summarise (
      n = n(),
      mean = mean(melt, na.rm=FALSE)
    ) %>%
    unlist() %>%
    {.} -> open_melt_AIS

  cat(open_melt_AIS["n"], "experiments with open melt\n", file = e$log_file)
  cat("Mean (values to be replaced):", open_melt_AIS["mean"], "\n", file = e$log_file)

  # Replace open melt values (not very elegant)
  if ( old_data == TRUE ) {
    fd <- fd %>%
      mutate(melt = replace(melt, which( ice_source == "AIS" & melt < 0.1), std_melt_AIS["mean"]))
  } else {
    fd <- fd %>%
      mutate(melt = replace(melt, which( ice_source == "AIS" & is.na(melt) ), std_melt_AIS["mean"]))
  }

  # Check
  fd %>%
    filter(ice_source == "AIS" ) %>%
    select(melt) %>%
    summarise (
      mean(melt, na.rm = TRUE)
    ) %>%
    unlist() %>%
    cat("Final mean of all AIS melt values (should match imputed value):", ., "\n", file = e$log_file)

  # Greenland: 'melt' is kappa parameter
  fd %>%
    filter(ice_source == "GrIS" ) %>%
    select(melt) %>%
    summarise (
      mean(melt, na.rm = TRUE)
    ) %>%
    unlist() %>%
    cat("\nInitial mean of all GrIS melt values:", . , "\n", file = e$log_file)

  # Get number and mean melt of standard melt experiments
  # BGC also did high and low which are treated as std, so need to specify medium melt value for open
  # old_data: list to exclude should match open_melt_models list in write_sealevel.R
  # i.e. GrIS experiments with melt=NA
  if ( old_data == TRUE ) {
    fd2 <- fd %>%
      filter(ice_source == "GrIS") %>%
      filter( ! (group == "UAF" & model == "PISM2"),
              ! (group == "UCIJPL" & model == "ISSM2"),
              ! (group == "VUW"),
              ! (group == "BGC" & abs(melt - -0.17) < 0.01) )
  } else {
    fd2 <- fd %>%
      filter(ice_source == "GrIS" & !is.na(melt) )
  }

  fd2 %>%
    summarise(
      n = n(),
      mean = mean(melt, na.rm=FALSE)
    ) %>%
    unlist() %>%
    {.} -> std_melt_GrIS

  cat(std_melt_GrIS["n"], "GrIS experiments with std melt\n", file = e$log_file)
  cat("Mean (imputed value):",std_melt_GrIS["mean"], "\n", file = e$log_file)

  # Get open melt experiments as above
  if ( old_data == TRUE ) {
    fd2 <- fd %>%
      filter(ice_source == "GrIS" & !is.na(melt)) %>%
      filter( (group == "BGC" & abs(melt - -0.17) < 0.01)
              | (group == "UAF" & model == "PISM2")
              | (group == "UCIJPL" & model == "ISSM2")
              | (group == "VUW") )
  } else {
    fd2 <- fd %>%
      filter(ice_source == "GrIS" & is.na(melt))
  }
  # Get number and mean of open melt
  fd2 %>%
    summarise(
      n = n(),
      mean = mean(melt, na.rm=FALSE)
    ) %>%
    unlist() %>%
    {.} -> open_melt_GrIS

  cat(open_melt_GrIS["n"], "experiments with open melt\n", file = e$log_file)
  cat("Mean (values to be replaced):",open_melt_GrIS["mean"], "\n", file = e$log_file)

  # Set values
  if ( old_data == TRUE ) {
    fd <- fd %>%
      mutate(melt = case_when(ice_source == "GrIS" & !is.na(melt)
                              & ( (group == "UAF" & model == "PISM2")
                                  | (group == "UCIJPL" & model == "ISSM2")
                                  | (group == "VUW") ) ~ std_melt_GrIS["mean"],
                              ice_source == "GrIS" & !is.na(melt)
                              & group == "BGC" & abs(melt - -0.17) < 0.01 ~ std_melt_GrIS["mean"],
                              TRUE ~ melt ) )
  } else {

    fd <- fd %>%
      mutate(melt = case_when(ice_source == "GrIS" & is.na(melt) ~ std_melt_GrIS["mean"],
                              TRUE ~ melt ) )
  }

  # Check ensemble mean matches imputed
  fd %>%
    filter(ice_source == "GrIS" ) %>%
    select(melt) %>%
    summarise (
      mean(melt, na.rm = TRUE)
    ) %>%
    unlist() %>%
    cat("\nFinal mean of all GrIS melt values:", . , "\n", file = e$log_file)

  # Save for later
  e$open_melt <- list()
  e$open_melt[["GrIS"]] <- std_melt_GrIS["mean"]
  cat("\nMELT VAL",e$open_melt[["GrIS"]] , file = e$log_file)
  e$open_melt[["AIS"]] <- std_melt_AIS["mean"]

  cat("\nRemaining:", dim(fd)[1], "sea level simulations\n\n", file = e$log_file)

  # Old dataset also had pre-2015 'historical' simulations
  if ( old_data == TRUE ) {
    cat("Exclude hist scenarios if present\n", file = e$log_file)
    fd <- filter(fd, ! scenario == "hist" )
  }



  # write final --------------------------------------

  cat("\nread_sealevel: FINAL DATA\n\n", file = e$log_file)

  cat("Sea level equivalent (SLE) change since 2015:\n\n", file = e$log_file)
  cat("Note: Glacier summary doesn't account for model number diffs across regions\n", file = e$log_file)

  fd %>%
    select(ice_source, region, group, model, exp_id, scenario, y2100) %>%
    group_by(ice_source, scenario) %>%
    summarise(n = n(),
              n_miss2100 = length( which(is.na(y2100)) ),
              q05 = quantile(y2100, prob = 0.05, na.rm = TRUE),
              median = median(y2100, na.rm = TRUE),
              q95 = quantile(y2100, prob = 0.95, na.rm = TRUE),
              mean = mean(y2100, na.rm = TRUE),
              sd = sd(y2100, na.rm = TRUE)
    ) %>%
    {.} -> sealevel_summ
  for (ii in 5:9) sealevel_summ[,ii] <- format(unlist(sealevel_summ[,ii]), digits = 3)
  suppressWarnings(write.table(sealevel_summ, file = e$log_file, append = TRUE, quote = FALSE, row.names = FALSE, sep = "\t"))

  cat("\nReturning:", dim(fd)[1], "sea level simulations\n\n", file = e$log_file)

  cat("\nRescaling AIS gamma: dividing by scale factor", e$sc, "\n\n", file = e$log_file)

  # Divide AIS melt by scale factor - purely for convenience (readability, plot axes etc)
  fd <- fd %>%
    mutate(melt = case_when(ice_source == "AIS" ~ melt/e$sc,
                            TRUE ~ melt ) )
  e$open_melt[["AIS"]] <- e$open_melt[["AIS"]]/e$sc

  # Save list of ice sheet / glacier model names for each region for emulator
  if (e$add_dummy %in% c("group", "model")) {
    e$ice_models <- list()
    for (is in e$ice_source_list) {

      for (reg in e$region_list[[is]] ) {
        e$ice_models[[reg]] <- unique(select(filter(fd, ice_source == is, region == reg), e$add_dummy))
        cat("\n\n", e$add_dummy, "names for", reg, ":\n", file = e$log_file)
        cat(paste(unlist(e$ice_models[[reg]]), collapse = ", "), file = e$log_file)

      }
    }
  }

  # Add dummy variables for model or group name (note merges same model name across different groups)
  if (e$add_dummy %in% c("group", "model", "melt")) {

    print(paste("Creating dummy variable columns for", e$add_dummy, "names:"))

    if (e$add_dummy %in% c("group", "model")) {
      fd2 <- dummies::dummy.data.frame(as.data.frame(fd), names = e$add_dummy, sep = "_", verbose = TRUE)

      # Model column was deleted above so put it back (at front to keep run info together)
      fd <- cbind(select(fd, e$add_dummy), fd2)
    }

    # Create dummy colummn for open vs standard melt
    if (e$add_dummy == "melt") {

      # Set to 1 when open melt
      melt0 <- rep( 0, dim(fd)[1] )
      melt0[ (fd$ice_source == "GrIS" & abs(fd$melt - e$open_melt[["GrIS"]]) < 1e-6)
             | (fd$ice_source == "AIS" & abs(fd$melt - e$open_melt[["AIS"]]) < 1e-6)  ] <- 1

      # Add to front
      fd <- cbind( melt0, fd )

    }

  }

  # Number of initial columns before numeric data columns
  e$ncol_param <- dim(fd)[2] - length(e$years_data)

  cat("\n\nNumber of columns before data:", e$ncol_param, "\n", file = e$log_file)

  # Return SLE dataset
  fd

}

#_______________________________________________________
# READ IN MELT PRIORS FOR ICE SHEET REGIONS
#_______________________________________________________

read_melt <- function(gamma0_prior) {

  #' Read in parameters
  #' @param gamma0_prior Prior distribution for AIS melt parameter

  cat("\nread_melt --------------------------------------\n", file = e$log_file)

  cat("\nAntarctic prior:", gamma0_prior, file = e$log_file)

  # Keep distributions for sampling and density estimates for plots
  e$melt_prior <- list()
  e$melt_prior_dens <- list()

  # GREENLAND
  # K-distribution from Donald Slater (N = 191)
  # Sample from kernel density estimate of this
  k_dist.file <- system.file("extdata", "data_for_tamsin.txt", package = e$packagename, mustWork = TRUE )
  k_dist <- read.table(k_dist.file)
  e$melt_prior_dens[["GrIS"]] <- density(k_dist[,1], n = 10000, bw = 0.0703652) # Auto is 0.07262901
  e$melt_prior[["GrIS"]] <- sample(e$melt_prior_dens[["GrIS"]]$x, 10000, replace = TRUE,
                                   prob = e$melt_prior_dens[["GrIS"]]$y)

  # ANTARCTICA
  # gamma0 distributions from Nico Jourdain (N = 10000 each)
  # Sample from smoothed kernel density estimate of combined distribution (or individual distributions for SA)
  # Truncate at zero i.e. remove these samples
  gamma0_MeanAnt.file <- system.file("extdata", "output_gamma0_NonLocal_MeanAnt.dat", package = e$packagename, mustWork = TRUE )
  gamma0_PIGL.file <- system.file("extdata", "output_gamma0_NonLocal_PIGL.dat", package = e$packagename, mustWork = TRUE )

  gamma0_MeanAnt <- read.table(gamma0_MeanAnt.file)
  gamma0_PIGL <- read.table(gamma0_PIGL.file)

  gamma0_all <- c(gamma0_MeanAnt[,1], gamma0_PIGL[,1])

  if (gamma0_prior == "joint") {

    # Expert judgement: merge with 3x bandwidth kde to fill the gap, truncate at zero
    e$melt_prior_dens[["AIS"]] <- density(gamma0_all, n = 10000, adjust = 3)
    e$melt_prior[["AIS"]] <- sample(e$melt_prior_dens[["AIS"]]$x, 10000, replace = TRUE,
                                    prob = e$melt_prior_dens[["AIS"]]$y)
    e$melt_prior[["AIS"]] <- e$melt_prior[["AIS"]][ e$melt_prior[["AIS"]] >= 0 ]

  }

  # Or use individual distributions for sensitivity analysis
  if (gamma0_prior %in% c("MeanAnt", "PIGL")) {
    if (gamma0_prior == "MeanAnt") e$melt_prior[["AIS"]] <- gamma0_MeanAnt[,1]
    if (gamma0_prior == "PIGL") e$melt_prior[["AIS"]] <- gamma0_PIGL[,1]
    #    e$melt_prior[["AIS"]] <- sample(gamma0_all, 10000, replace = TRUE) # why did I resample?
  }

  if (gamma0_prior == "unif") {
    e$melt_prior[["AIS"]] <- runif(10000, min = e$melt_values[["AIS"]][["PIGL_med"]]*e$sc, max = e$melt_values[["AIS"]]["PIGL_high"]*e$sc)
  }

  # This kind of range necessary to reproduce LARMIP median:
  if (gamma0_prior == "unif_high") {
    e$melt_prior[["AIS"]] <- runif(10000, min = e$melt_values[["AIS"]][["PIGL_med"]]*e$sc, max = 700.0*e$sc)
  }

  cat("\nRange of AIS gamma prior values before and after rescaling:\n", file = e$log_file)
  cat(range(e$melt_prior[["AIS"]]), "\n", file = e$log_file)

  # Rescale, as for inputs
  e$melt_prior[["AIS"]] <- e$melt_prior[["AIS"]]/e$sc
  cat(range(e$melt_prior[["AIS"]]), "\n", file = e$log_file)

  # Plot empirical estimate: Greenland (data range -5.3 to 3.0)
  sbin <- 0.05
  hist( k_dist[,1], xlim = c(-6, 3), freq = FALSE, breaks = seq(-10, 10, by = sbin),
        xlab = expression(paste("Greenland glacier retreat parameter, ", kappa)),
        main = "Empirical distribution",
        col = "cornflowerblue", border = "darkblue")
  abline( v = e$melt_values[["GrIS"]], lty = c(2,1,2,3,3), lwd = 2)
  lines( e$melt_prior_dens[["GrIS"]], lwd = 2)

  # MAIN: kappa --------------------------------------

  # FIGURE: 1c GREENLAND MELT PRIOR
  pdf( file = paste0( e$outdir, "/MAIN_kappa.pdf" ), width = 9, height = 8)
  par(mar = c(6, 5, 1.5, 2.5))

  xlab <- expression("Greenland glacier retreat parameter," ~ kappa ~
                       "(km ("*m^3 ~ s^-1 *")"^-0.4 ~ degree*C^-1 * ")")
  hist( unlist(e$melt_prior[["GrIS"]]), xlab = xlab,
        xlim = e$range_melt[["GrIS"]], xaxs = "i", breaks = seq(from = -6, to = 6, by = sbin),
        main = NULL, col = "cornflowerblue", border = "cornflowerblue",
        cex.axis = 1.6, cex.lab = 2 )
  abline( v = e$melt_values[["GrIS"]], lty = c(2,1,2, 3, 3), lwd = 2)
  text( 0.3, 1060, pos = 4, font = 2, "b", cex = 2.5)
  dev.off()

  # Plot empirical estimates: Antarctica (data range 5200 to 1570000)
  sbin <- 20000
  hist( gamma0_all, xlim = c(0, 1600000), freq = FALSE,
        breaks = seq(0, 2000000, by = sbin), xlab = expression(paste("Antarctic basal melt parameter, ",gamma[0])),
        main = "Empirical distribution",
        col = "cornflowerblue", border = "darkblue")
  abline(v = e$melt_values[["AIS"]][c("low", "med", "high")], lty = c(2,1,2), lwd = 2) # MeanAnt
  abline(v = e$melt_values[["AIS"]][c("PIGL_low", "PIGL_med", "PIGL_high")], lty = c(2,1,2), col = "grey", lwd = 2)
  lines( e$melt_prior_dens[["AIS"]], lwd = 2)

  # MAIN: gamma --------------------------------------

  # FIGURE: 1d ANTARCTICA MELT PRIOR
  sbin <- sbin/e$sc
  pdf( file = paste0( e$outdir, "/MAIN_gamma.pdf" ), width = 9, height = 8)
  par(mar = c(6, 5, 1.5, 2.5))

  xlab <- expression(paste("Antarctic basal melt parameter, ", gamma, " (x",10^3," m ",a^-1,")"))

  hist( unlist(e$melt_prior[["AIS"]]),
        xlab = xlab,
        xlim = e$range_melt[["AIS"]], xaxs = "i",
        breaks = seq(from = 0, to = 2000000/e$sc, by = sbin),
        main = NULL, col = "cornflowerblue", border = "cornflowerblue", cex.axis = 1.6, cex.lab = 2)
  abline(v = e$melt_values[["AIS"]][c("low", "med", "high")], lty = c(2,1,2), lwd = 2) # MeanAnt
  abline(v = e$melt_values[["AIS"]][c("PIGL_low", "PIGL_med", "PIGL_high")], lty = c(2,1,2), col = "grey", lwd = 2) # PIGL
  text( 650, 950, pos = 4, font = 2, "c", cex = 2.5)
  dev.off()

}
