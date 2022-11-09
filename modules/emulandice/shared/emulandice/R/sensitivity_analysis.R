sensitivity_analysis <- function(year, ice_source, region, sim, emulator, N_melt = 1) {

  #' Sensitivity analysis: dependence of results on temperature and ice sheet parameters
  #' @param year Year
  #' @param ice_source Ice source
  #' @param region Ice source region
  #' @param sim Simulation data for plotting
  #' @param emulator Emulator object for predicting
  #' @param N_melt Number of samples in plots sampling melt (not used in MAIN plots)

  yy <- year
  yy_num <- substr(yy, 2, nchar(yy))

  # Text labels for all GSAT and SLE axes
  tlab <- bquote("Global mean temperature change 2015-"*.(yy_num)~"("*degree*"C)")
  slab <- bquote("Sea level contribution at"~.(yy_num)~"(cm SLE)")
  is <- ice_source
  reg <- region
  reg_name <- e$region_name_list[[is]][which(e$region_list[[is]] == reg, arr.ind = TRUE)]
  print(paste("\nSensitivity Analysis:", reg_name))

  # Get ice_model list
  if (e$add_dummy %in% c("model", "group")) {
    ice_models <- paste(e$add_dummy, unlist(e$ice_models[[reg]])[-1], sep = "_")
    nd <- length(ice_models) # Number of dummy variables (n_model - 1)
  }

  # T dep --------------------------------------

  # Plot ranges (range FAIR at 2100 is -0.23 to 7.9)
  temp_vec <- seq(-0.5, 8, by = 0.1)

  # MAIN: T dep --------------------------------------
  # T-dependence of SLE at 2100

  # Plot filename
  pdf( file = paste0( e$outdir, "/MAIN_Tdep_", reg,".pdf"), width = 9, height = 8)
  par(mar = c(5, 6, 2, 1))

  # START PLOT
  plot( 1:3, 1:3, type = "n", xlim = range(temp_vec), ylim = e$range_sle[[reg]],
        xaxs = "i", yaxs = "i",
        xlab = tlab, ylab = slab,
        main = NULL, cex = 1.5, cex.main = 1.5, cex.axis = 1.6, cex.lab = 1.8 )
  abline( h = 0)

  # Prior samples over all models
  if ( is %in% c("GrIS", "AIS") && e$add_dummy == "melt" ) {
    if ( abs(e$dummy_melt_pred - 0.5) < 0.01 ) dummy_melt <- sample(c(0,1), length(temp_vec), replace = TRUE)
    else dummy_melt <- rep(e$dummy_melt_pred, length(temp_vec))
    dummy_prior <- data.frame(melt0 = dummy_melt)
  }
  if (e$add_dummy %in% c("model", "group")) {
    dummy_prior <- as.data.frame( t(replicate( length(temp_vec), sample( c(1,rep(0,nd)), nd, replace = FALSE ))) )
    names(dummy_prior) <- ice_models
  }

  # Glacier prediction (temp is only input)
  if (is == "Glaciers") {

    temp_vec_df <- data.frame(temp = temp_vec)
    if (e$add_dummy %in% c("model", "group")) temp_vec_df <- cbind( temp_vec_df, dummy_prior)

    temp_vec_df_scaled <- as.data.frame(scale(temp_vec_df, center = e$input_centre[[reg]], scale = e$input_scale[[reg]]))

    if (e$emul_type == "DK") pred_tdep <- DiceKriging::predict( emulator, newdata = temp_vec_df_scaled,
                                                                type = "UK", checkNames = TRUE )
    if (e$emul_type == "RG") {
      temp_vec_mat <- as.matrix(temp_vec_df_scaled) # xxx improve: convoluted
      trend.test.rgasp <- cbind(rep(1,dim(temp_vec_mat)[1]), temp_vec_mat)
      pred_tdep <- RobustGaSP::predict( emulator, temp_vec_mat, testing_trend = trend.test.rgasp)
    }

    lty = 1
    col1 <- grey(0.05, alpha = 0.05)
    lines( temp_vec, pred_tdep$mean, lty = lty, lwd = 0.5)
    lines( temp_vec, pred_tdep$mean + 2*pred_tdep$sd, lty = lty, lwd = 0.5)
    lines( temp_vec, pred_tdep$mean - 2*pred_tdep$sd, lty = lty, lwd = 0.5)
    polygon(c(temp_vec, rev(temp_vec)),
            c(pred_tdep$mean + 2*pred_tdep$sd,
              rev(pred_tdep$mean - 2*pred_tdep$sd)), border = col1,
            col = col1)

    # Add maximum contributionfrom Farinott et al. (2 methods)
    abline(h = e$max_glaciers[[reg]][2], lty = 5, col = "darkgrey")

  }

  # Ice sheet prediction: fixed melt and FIXED collapse = 0
  if (is %in% c("GrIS", "AIS")) {

    if (is == "GrIS") melt_val_list <- c( "med", "high")
    if (is == "AIS") melt_val_list <- c("med", "PIGL_med")
    stopifnot(length(melt_val_list) <= 2)

    col1 <- grey(0.05, alpha = 0.05)
    col2 <- grey(0.05, alpha = 0.03)

    # Loop through multiple fixed melt values
    for (melt_val in melt_val_list) {

      cat("Predicting for melt", melt_val,"\n")
      cat(e$melt_values[[is]][[ melt_val ]], "\n")

      # Collapse random: sample(c(0,1), length(temp_vec), replace = TRUE)
      fixed_melt <- data.frame( temp = temp_vec,
                                melt = rep( e$melt_values[[is]][[ melt_val ]], length(temp_vec) ),
                                collapse = rep( 0, length(temp_vec)) )
      if (e$add_dummy %in% c("model", "group", "melt")) fixed_melt <- cbind( fixed_melt, dummy_prior)

      # Drop terms
      fixed_melt <- select( fixed_melt, names(e$input) )

      summarise(tibble(fixed_melt), mean = mean(.data[["temp"]]), sd = sd(.data[["temp"]])) %>%
        unlist %>%
        cat("Unscaled: Mean and s.d. of temp:",.,"\n", file = e$log_file)

      # Scale design when making prediction
      cat("Rescaling design\n", file = e$log_file)
      fixed_melt_scaled <- as.data.frame(scale(fixed_melt, center = e$input_centre[[reg]], scale = e$input_scale[[reg]]))

      summarise(tibble(fixed_melt_scaled), mean = mean(.data[["temp"]]), sd = sd(.data[["temp"]])) %>%
        unlist %>%
        cat("Scaled: Mean and s.d. of temp:",.,"\n", file = e$log_file)

      if (e$emul_type == "DK") pred_tdep <- DiceKriging::predict( emulator, newdata = fixed_melt_scaled, type = "UK", checkNames = TRUE )
      if (e$emul_type == "RG") {
        fixed_melt_scaled <- as.matrix(fixed_melt_scaled)
        trend.test.rgasp <- cbind(rep(1,dim(fixed_melt_scaled)[1]), fixed_melt_scaled )
        pred_tdep <- RobustGaSP::predict( emulator, fixed_melt_scaled, testing_trend = trend.test.rgasp)
      }

      # Emulator mean +/- s.d. first
      if (melt_val == melt_val_list[1]) {
        lty = 1
        col = col1
      } else {
        lty = 3
        col = col2
      }

      # Plot emulator predictions
      lines( temp_vec, pred_tdep$mean, lty = lty, lwd = 0.5)
      lines( temp_vec, pred_tdep$mean + 2*pred_tdep$sd, lty = lty, lwd = 0.5)
      lines( temp_vec, pred_tdep$mean - 2*pred_tdep$sd, lty = lty, lwd = 0.5)
      polygon(c(temp_vec, rev(temp_vec)),
              c(pred_tdep$mean + 2*pred_tdep$sd,
                rev(pred_tdep$mean - 2*pred_tdep$sd)), border = col,
              col = col)

    } # melt prediction loop

  } # ice sheets

  # MODEL DATA: all for Glaciers
  for ( ss in 1:dim(sim)[1]) { # xxx improve: vectorise?

    # Select only fixed melt same values, and collapse off, for ice sheets
    if ( is %in% c("GrIS", "AIS")) {
      pch = NA
      if ( abs( sim[ss, "melt"] - e$melt_values[[is]][[ melt_val_list[1] ]] ) < 0.01 ) pch <- 16
      if ( abs( sim[ss, "melt"] -  e$melt_values[[is]][[ melt_val_list[2] ]] ) < 0.01 ) pch <- 1 # open
      if (is == "AIS") if ( sim[ss, "collapse"] == 1) pch = NA # don't show collapse = on
    } else {
      # Plot all glacier simulations
      pch <- 16
    }

    # Plot
    if ( !is.na(pch) ) {
      sc <- sim[ss, "scenario"] %>%
        substr(., nchar(.)-1, nchar(.))
      points( sim[ ss, "temp" ], sim[ ss, yy ],
              pch = pch, cex = 1.8, col = e$scen_col_trans[[sc]] )

    }
  }

  # Legend
  text(-0.12,
       e$range_sle[[reg]][1] + 0.95*(e$range_sle[[reg]][2]-e$range_sle[[reg]][1]),
       reg_name, pos = 4, cex = 1.8 )
  frac <- 0.88
  for (sc in names(e$scen_col_trans)) {
    scenario <- e$rcp_name_list[[sc]]
    if (length(scenario) > 0 && scenario %in% e$scen_name_list[["CMIP5"]]) {
      # Ice sheets don't have middle RCPs and do have equivalent SSPs
      if (is %in% c("GrIS", "AIS")) {
        if (scenario %in% c("RCP4.5", "RCP6.0")) next
        if (scenario == "RCP2.6") scenario <- paste(scenario, "/ SSP1-26")
        if (scenario == "RCP8.5") scenario <- paste(scenario, "/ SSP5-85")
      }
      text( 0.1,
            e$range_sle[[reg]][1] + frac*(e$range_sle[[reg]][2]-e$range_sle[[reg]][1]),
            scenario, pos = 4, cex = 1.6)
      points( 0.0,
              e$range_sle[[reg]][1] + frac*(e$range_sle[[reg]][2]-e$range_sle[[reg]][1]),
              pch = 20, cex = 1.8, col = e$scen_col_trans[[sc]])
      frac <- frac - 0.05
    }
  }

  # Add covariance to legend
  leg_cov <- e$covariance_fn[[reg]]
  if (leg_cov == "pow_exp") leg_cov = paste(leg_cov, "p =", e$alpha_powexp[[reg]])
  #  text(0, frac * e$range_sle[[reg]][2], leg_cov, pos = 4, cex = 1.2)

  subfig = ""
  if (reg == "ALL") subfig <- "a"
  if (reg == "WAIS") subfig <- "b"
  if (reg == "EAIS") subfig <- "c"
  if (reg == "region_5") subfig <- "d"
  if (reg == "PEN") subfig <- "e"
  if (reg == "region_19") subfig <- "f"
  if (reg == "region_1") subfig <- "g"
  if (reg == "region_3") subfig <- "h"
  if (reg == "region_4") subfig <- "i"
  if (reg == "region_9") subfig <- "j"
  if (reg == "region_11") subfig <- "k"
  if (reg == "region_12") subfig <- "l"
  if (reg == "region_13") subfig <- "m"
  if (reg == "region_14") subfig <- "n"
  if (reg == "region_15") subfig <- "o"

  text( 0.92 * max(temp_vec),
        e$range_sle[[reg]][1] + 0.06 * (e$range_sle[[reg]][2] - e$range_sle[[reg]][1]),
        pos = 4, font = 2, subfig, cex = 1.8)


  dev.off()

  # T dep, melt --------------------------------------
  # Ice sheets only: T-dependence with melt uncertainty

  if (is %in% c("GrIS", "AIS") && N_melt > 1) {

    temp_vec <- seq(0, 7) # Only need to do to 7

    # Factorial design of temp vector and melt prior samples
    vary_melt <- sample( unlist(e$melt_prior[[is]]), N_melt, replace = TRUE ) %>%
      expand.grid(temp_vec, .) %>%
      as.data.frame() %>%
      rename(temp = Var1, melt = Var2)

    # Add collapse column of same length
    if (is == "AIS") {
      vary_melt <- cbind(vary_melt, sample(c(0,1), length(temp_vec) * N_melt, replace = TRUE))
      colnames(vary_melt)[3] <- "collapse"
    }

    # Prior samples over all models
    # Longer sample this time
    if (e$add_dummy %in% c("model", "group", "melt")) {

      if (e$add_dummy %in% c("model", "group")) {
        dummy_prior <- as.data.frame( t(replicate( length(temp_vec) * N_melt, sample( c(1,rep(0,nd)), nd, replace = FALSE ))) )
        names(dummy_prior) <- ice_models
      }
      if (e$add_dummy == "melt") {
        if ( abs(e$dummy_melt_pred - 0.5) < 0.01 ) dummy_melt <- sample(c(0,1), length(temp_vec) * N_melt, replace = TRUE)
        else dummy_melt <- rep(e$dummy_melt_pred, length(temp_vec) * N_melt)
        dummy_prior <- data.frame(melt0 = dummy_melt)
      }

      vary_melt <- cbind( vary_melt, dummy_prior )

    }

    # Drop terms
    vary_melt <- select(vary_melt, names(e$input))

    vary_melt_scaled <- as.data.frame(scale(vary_melt, center = e$input_centre[[reg]], scale = e$input_scale[[reg]]))

    # Predict
    if (e$emul_type == "DK") pred_tdep <- DiceKriging::predict( emulator, newdata = vary_melt_scaled, type = "UK", checkNames = TRUE )
    if (e$emul_type == "RG") {
      vary_melt_mat <- as.matrix(vary_melt_scaled)
      trend.test.rgasp <- cbind(rep(1,dim(vary_melt_mat)[1]), vary_melt_mat )
      pred_tdep <- RobustGaSP::predict( emulator, vary_melt_mat, testing_trend = trend.test.rgasp)
    }

    # Plot
    boxplot( sle ~ temp, data = data.frame(temp = select(vary_melt,"temp"), sle = pred_tdep$mean),
             boxwex = 0.3, outpch = NA, col = grey(0.05, alpha = 0.05),
             border = grey(0.8),
             xlim = c(-0.5, e$max_temp[[yy_num]]), ylim = e$range_sle[[reg]], at = temp_vec,
             xlab = tlab, ylab = slab,
             main = paste0(reg_name, " ", yy_num, ": melt uncertainty" ) )
    abline( h = 0)

    # MODEL DATA
    sim_plot <- sim

    # Loop for nice RCP colours and melt symbols xxx improve: vectorise?
    for ( ss in 1:dim(sim_plot)[1]) {
      pch = 16 # closed circle
      if (is == "AIS") {
        if (sim_plot[ss, "melt"] > 50000/e$sc) pch = 1 # open circle for PIGL values
        if (sim_plot[ss, "collapse"] == 1) pch = 8 # star for collapse
      }
      sc <- sim_plot[ss, "scenario"] %>%
        substr(., nchar(.)-1, nchar(.))
      points( sim_plot[ ss, "temp" ], sim_plot[ ss, yy ],
              pch = pch, cex = 0.8, col = e$scen_col_trans[[sc]])
    }

    # Legend
    frac <- 0.95
    for (sc in names(e$scen_col_trans)) {

      scenario <- e$rcp_name_list[[sc]]

      if (length(scenario) > 0 && scenario %in% e$scen_name_list[["CMIP5"]]) {

        # Ice sheets don't have middle RCPs and do have equivalent SSPs
        if (is %in% c("GrIS", "AIS")) {
          if (scenario %in% c("RCP4.5", "RCP6.0")) next
          if (scenario == "RCP2.6") scenario <- paste(scenario, "/ SSP1-26")
          if (scenario == "RCP8.5") scenario <- paste(scenario, "/ SSP5-85")
        }
        text( 0, frac * e$range_sle[[reg]][2], scenario, pos = 4, cex = 1.2)
        points(0, frac * e$range_sle[[reg]][2], pch = 20, cex = 0.9, col = e$scen_col_trans[[sc]])
        frac <- frac - 0.04
      }
    }

    # T dep, melt, em --------------------------------------

    # Ice sheets only: T-dependence with melt and emulator uncertainties
    print("ED Figure 2")

    # Create data frame for e$smid density estimates per T value
    # em_tt_mean_all -> pred_tdep_em
    pred_tdep_em <- expand.grid(temp_vec, rep(NA, length(e$smid)) ) %>%
      as.data.frame() %>%
      rename(temp = Var1, dens = Var2)

    # Estimate pdf using mid-point rule separately for each T value
    for (tt in temp_vec) {

      # Get temp_vec value: ta_melt_grid$gsat -> vary_melt$temp
      tind <- abs(vary_melt$temp - tt) < 0.001

      # xxx improve: maybe shouldn't overwrite name?
      # xxx improve: use tidyverse filter()?
      pred_tdep_mean <- pred_tdep$mean[ tind ]
      pred_tdep_sd <- pred_tdep$sd[ tind ]
      stopifnot( length(pred_tdep_mean) == N_melt && length(pred_tdep_sd) == N_melt)

      # Start with zero density
      melt_int <- rep(0, length(e$smid))

      # Sum of individual pdfs (there are N_melt predictions for each temp_vec value)
      # xxx improve: vectorise?
      for ( mm in 1:N_melt) {

        # Add density for this prediction (normal dist with emulator mean and s.d.)
        melt_int <- melt_int + dnorm(e$smid, mean = pred_tdep_mean[ mm ],
                                     sd = pred_tdep_sd[ mm ] )
      }

      # Normalise by sum and save
      # xxx improve: filter()?
      pred_tdep_em[ abs(pred_tdep_em$temp - tt) < 0.001, "dens" ] <- melt_int / sum(melt_int)

    } # temperature values

    pdf( file = paste0( e$outdir, "/ED_Fig2_", reg,".pdf"), width = 9, height = 8)

    par(mar = c(5, 6, 1, 1))

    # N.B. Not boxplot() like previous plot because pdf, not sample
    plot( 1:3, 1:3, type = "n",
          xlim = c(-0.5, e$max_temp[[yy_num]]), ylim = e$range_sle[[reg]], xaxs = "i", yaxs = "i",
          xlab = tlab, ylab = slab,
          main = NULL, cex = 1.5, cex.main = 1.5, cex.axis = 1.4, cex.lab = 1.5 )
    abline( h = 0)
    tinc <- 0.2 # box plot width

    for (tt in temp_vec) {

      # Estimate quantiles from pdf
      qq <- pred_tdep_em[ abs(pred_tdep_em$temp - tt) < 0.001, "dens" ] %>%
        cumsum(.) %>%
        approx( e$smid, xout = c(0.05, 0.25, 0.50, 0.75, 0.95) )
      qq <- qq$y

      # Box and whisker
      lines( rep(tt, 2), qq[c(1,2)], col = grey(0.8))
      lines( rep(tt, 2), qq[c(4,5)], col = grey(0.8))
      rect( tt - tinc, qq[2], tt + tinc, qq[4], lwd = 1.5,
            col = grey(0.05, alpha = 0.05), border = grey(0.8))
      lines( tt + tinc*c(-1,1), rep(qq[3], 2), col = grey(0.8), lwd = 2)
      lines( tt + tinc*c(-0.5,0.5), rep(qq[1], 2), col = grey(0.8))
      lines( tt + tinc*c(-0.5,0.5), rep(qq[5], 2), col = grey(0.8))

    }

    # Loop for nice RCP colours and melt symbols
    sim_plot <- sim
    for ( ss in 1:dim(sim_plot)[1]) {
      pch = 16 # closed
      if (is == "AIS") {
        if (sim_plot[ss, "melt"] > 50000/e$sc) pch = 1 # open circle for 3 PIGL values # xxx improve: match melt_values
        if (sim_plot[ss, "collapse"] == 1) pch = 8 # star for collapse
      }
      sc <- sim_plot[ss, "scenario"] %>%
        substr(., nchar(.)-1, nchar(.))
      points( sim_plot[ ss, "temp" ], sim_plot[ ss, yy ],
              pch = pch, cex = 0.8, col = e$scen_col_trans[[sc]])
    }

    # Legend
    text(0, 0.90*e$range_sle[[reg]][2], reg_name, pos = 4, cex = 1.4 )
    frac <- 0.80
    for (sc in names(e$scen_col_trans)) {
      scenario <- e$rcp_name_list[[sc]]
      if (length(scenario) > 0 && scenario %in% e$scen_name_list[["CMIP5"]]) {
        # Ice sheets don't have middle RCPs and do have equivalent SSPs
        if (scenario %in% c("RCP4.5", "RCP6.0")) next
        if (scenario == "RCP2.6") scenario <- paste(scenario, "/ SSP1-26")
        if (scenario == "RCP8.5") scenario <- paste(scenario, "/ SSP5-85")
        text( 0, frac * e$range_sle[[reg]][2], scenario, pos = 4, cex = 1.2)
        points(0, frac * e$range_sle[[reg]][2], pch = 20, cex = 0.9, col = e$scen_col_trans[[sc]])
        frac <- frac - 0.05
      }
    }

    # Box and whisker for legend
    tinc <- 0.02 * e$range_sle[[reg]][2]
    lines( c(0, 0.25), rep(frac * e$range_sle[[reg]][2], 2), col = grey(0.8))
    lines( c(0.75, 1), rep(frac * e$range_sle[[reg]][2], 2), col = grey(0.8))
    rect( 0.25, frac * e$range_sle[[reg]][2] - tinc,
          0.75, frac * e$range_sle[[reg]][2] + tinc, lwd = 1.5,
          col = grey(0.05, alpha = 0.05), border = grey(0.8))
    lines( rep(0.5, 2), frac * e$range_sle[[reg]][2] + tinc*c(-1,1), col = grey(0.8), lwd = 2)
    lines( rep(0, 2), frac * e$range_sle[[reg]][2] + tinc*c(-0.5,0.5), col = grey(0.8))
    lines( rep(1, 2), frac * e$range_sle[[reg]][2] + tinc*c(-0.5,0.5),  col = grey(0.8))

    text( 1, frac * e$range_sle[[reg]][2], "Emulator", pos = 4, cex = 1.2)

    dev.off()

  } # ice sheets and N_melt > 1

  # MAIN: Melt dep --------------------------------------

  # Greenland, WAIS and EAIS
  if ("melt" %in% names(e$input)) {

    # PLOT VS MELT AND COLLAPSE PARAMETERS

    # GrIS: 8.5	MIROC5; exps 5, 9, 10
    # AIS: 8.5	NorESM1-M collapse = FALSE
    # exps (1) 5, 9, 10, 13, D51, D52

    # Single GCM, scenario and collapse so only comparing melt perturbations (if all 3-4 there)
    if (is == "GrIS") {
      sim_plot <- filter(sim, ice_source == "GrIS" & scenario == "RCP85" & GCM == "MIROC5" )
      xlab <- expression("Greenland glacier retreat parameter," ~ kappa ~
                           "(km ("*m^3 ~ s^-1 *")"^-0.4 ~ degree*C^-1 * ")")
    }
    if (is == "AIS") {
      sim_plot <- filter(sim, ice_source == "AIS" & scenario == "RCP85" & GCM == "NorESM1-M" & collapse == 0)
      xlab <- expression(paste("Antarctic basal melt parameter, ", gamma, " (x",10^3," m ",a^-1,")"))
    }

    # Predict for different melt and same fixed RCP+GCM temperature
    temp_med <- unlist(unique(select(sim_plot, temp)))
    stopifnot(length(temp_med) == 1)

    melt_vec <- seq(from = e$range_melt[[is]][1], to = e$range_melt[[is]][2], length = 40)

    # Combine columns and reject any not used
    fixed_temp <- data.frame( temp = rep(temp_med, length(melt_vec)),
                              melt = melt_vec,
                              collapse = rep(0, length(melt_vec)))

    # Prior samples over all models
    if (e$add_dummy %in% c("model", "group", "melt")) {
      if (e$add_dummy %in% c("model", "group")) {
        dummy_prior <- as.data.frame( t(replicate( length(melt_vec), sample( c(1,rep(0,nd)), nd, replace = FALSE ))) )
        names(dummy_prior) <- ice_models
      }
      if (e$add_dummy == "melt") {
        if ( abs(e$dummy_melt_pred - 0.5) < 0.01 ) dummy_melt <- sample(c(0,1), length(melt_vec), replace = TRUE)
        else dummy_melt <- rep(e$dummy_melt_pred, length(melt_vec))
        dummy_prior <- data.frame(melt0 = dummy_melt)
      }

      fixed_temp <- cbind(fixed_temp, dummy_prior)
    }

    # Drop terms
    fixed_temp <- select( fixed_temp, names(e$input) )

    # Scale to mean = 0, s.d. = 1
    fixed_temp_scaled <- as.data.frame(scale(fixed_temp, center = e$input_centre[[reg]], scale = e$input_scale[[reg]]))

    if (e$emul_type == "DK") pred_mdep <- DiceKriging::predict( emulator, newdata = fixed_temp_scaled, type = "UK", checkNames = TRUE )

    if (e$emul_type == "RG") {
      fixed_temp <- as.matrix(fixed_temp_scaled)
      trend.test.rgasp <- cbind(rep(1,dim(fixed_temp)[1]), fixed_temp)
      pred_mdep <- RobustGaSP::predict( emulator, fixed_temp, testing_trend = trend.test.rgasp)
    }


    pdf( file = paste0( e$outdir, "/MAIN_meltdep_", reg,".pdf" ), width = 9, height = 8)
    par(mar = c(5, 6, 1, 2))

    plot( 1:3, 1:3, type = "n", xlim = e$range_melt[[is]],
          ylim = e$range_sle[[reg]], yaxs = "i", xaxs = "i",
          xlab = xlab, ylab = paste("Sea level contribution at", yy_num, "(cm SLE)"),
          main = NULL, cex = 1.5, cex.main = 1.5, cex.axis = 1.5, cex.lab = 1.7 )
    abline( h = 0)
    if (is != "GrIS") abline(v=0)

    if (is == "GrIS") {
      text(-0.95, 0.95*e$range_sle[[reg]][2], reg_name, pos = 4, cex = 1.6 )
      frac <- 0.88
      text( -0.95, frac * e$range_sle[[reg]][2], "RCP8.5: MIROC5", pos = 4, cex = 1.4)
      frac <- frac - 0.06
      points(-0.91, frac * e$range_sle[[reg]][2], pch = 17, cex = 1.6, col = e$scen_col_trans[["85"]])
      text( -0.9, frac * e$range_sle[[reg]][2], "IMAUICE", pos = 4, cex = 1.4)
      frac <- frac - 0.06
      points(-0.91, frac * e$range_sle[[reg]][2], pch = 16, cex = 1.6, col = e$scen_col_trans[["85"]])
      text( -0.9, frac * e$range_sle[[reg]][2], "Other models", pos = 4, cex = 1.4)
      frac <- frac - 0.06
      points(-0.91, frac * e$range_sle[[reg]][2], pch = 3, cex = 1.6, col = e$scen_col_trans[["85"]])
      text( -0.9, frac * e$range_sle[[reg]][2], "Not using parameterisation", pos = 4, cex = 1.4)
      #    frac <- frac - 0.05
      #    text(-1.2, frac * e$range_sle[[reg]][2], leg_cov, pos = 4, cex = 1.2)
    }

    if (is == "AIS") {
      text( 20,
            e$range_sle[[reg]][1] + 0.95*(e$range_sle[[reg]][2]-e$range_sle[[reg]][1]), reg_name, pos = 4, cex = 1.6 )
      frac <- 0.9
      text( 20,
            e$range_sle[[reg]][1] + frac*(e$range_sle[[reg]][2]-e$range_sle[[reg]][1]),
            "RCP8.5: NorESM1-M, no shelf collapse", pos = 4, cex = 1.4)
      frac <- frac - 0.04
      points( 35,
              e$range_sle[[reg]][1] + frac*(e$range_sle[[reg]][2]-e$range_sle[[reg]][1]),
              pch = 16, cex = 1.6, col = e$scen_col_trans[["85"]])
      text( 40,
            e$range_sle[[reg]][1] + frac*(e$range_sle[[reg]][2]-e$range_sle[[reg]][1]),
            "Models", pos = 4, cex = 1.4)
      frac <- frac - 0.04
      points( 35,
              e$range_sle[[reg]][1] + frac*(e$range_sle[[reg]][2]-e$range_sle[[reg]][1]),
              pch = 3, cex = 1.6, col = e$scen_col_trans[["85"]])
      text( 40,
            e$range_sle[[reg]][1] + frac*(e$range_sle[[reg]][2]-e$range_sle[[reg]][1]),
            "Models not using parameterisation", pos = 4, cex = 1.4)
      #    frac <- frac - 0.05
      #    text(0, frac * e$range_sle[[reg]][2], leg_cov, pos = 4, cex = 1.2)
    }


    # Emulator mean +/- s.d. first
    lines( melt_vec, pred_mdep$mean, lwd = 0.5)
    lines( melt_vec, pred_mdep$mean + 2*pred_mdep$sd, lwd = 0.5)
    lines( melt_vec, pred_mdep$mean - 2*pred_mdep$sd, lwd = 0.5)
    polygon(c(melt_vec, rev(melt_vec)),
            c(pred_mdep$mean + 2*pred_mdep$sd,
              rev(pred_mdep$mean - 2*pred_mdep$sd)),
            border = grey(0.05, alpha = 0.05),
            col = grey(0.05, alpha = 0.05))

    # Check at least 3 (Greenland) or 4 (Antarctica) melt values per model
    if (is == "GrIS") nm <- 3
    if (is == "AIS") nm <- 4

    for ( ss in 1:dim(sim_plot)[1]) {

      this_mod <- filter(sim_plot, group == sim_plot[ss, "group"] & model == sim_plot[ss, "model"])
      if (length(unique(this_mod[,"melt"])) < nm) next

      pch = 16 # closed circle
      col = e$scen_col_trans[[sc]] # semi-transparent

      if ( (abs(sim_plot[ ss, "melt"] - e$open_melt[[is]]) < 1e-6)
           || (is == "GrIS" && sim_plot[ss, "group"] == "BGC") ) pch = 3 # cross for open melt (not equal to any melt)

      # improve: use exp01 to select open melt for AIS??
      if (is == "GrIS" && sim_plot[ss, "model"] == "IMAUICE1") {
        col <- e$scen_col[[sc]]
        pch = 17
      }
      sc <- sim_plot[ss, "scenario"] %>%
        substr(., nchar(.)-1, nchar(.))

      points( sim_plot[ss, "melt"], sim_plot[ss , yy], pch = pch, cex = 1.6, col = col )

    }

    if (reg == "ALL") subfig <- "a"
    if (reg == "WAIS") subfig <- "b"
    if (reg == "EAIS") subfig <- "c"
    if (reg == "PEN") subfig <- "d"

    if (reg == "ALL") {
      text( 0.3,
            e$range_sle[[reg]][1] + 0.06 * (e$range_sle[[reg]][2] - e$range_sle[[reg]][1]),
            pos = 4, font = 2, subfig, cex = 1.6)

    } else {
      text( 650,
            e$range_sle[[reg]][1] + 0.06 * (e$range_sle[[reg]][2] - e$range_sle[[reg]][1]),
            pos = 4, font = 2, subfig, cex = 1.6)
    }
    dev.off()

  } # only ice sheet regions with melt

}
