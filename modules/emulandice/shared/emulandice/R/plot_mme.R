plot_mme <- function(is, reg) {

  #' Plots and outputs from raw simulation data
  #' @param is Ice source
  #' @param reg Ice source region

  # Output written to stats.txt

  to_plot <- filter(e$sle_data, ice_source == is & region == reg)
  reg_name <- e$region_name_list[[is]][ which(e$region_list[[is]] == reg ) ]
  print("", quote = FALSE)
  print("==============================================", quote = FALSE)
  print(paste("PLOT MME:", reg_name))
  print("==============================================", quote = FALSE)

  # MAIN: scenario dep --------------------------------------

  if (is == "GrIS" || (is == "AIS" && reg != "ALL") ) {

    print("", quote = FALSE)
    print("____________________________________", quote = FALSE)
    print("Scenario dependence quadrant plots", quote = FALSE)

    if (is == "GrIS") ylim <- c(-10,25)

    if (is == "AIS") {
      if (reg == "WAIS") ylim <- c(-20,25)
      if (reg == "EAIS") ylim <- c(-8,8)
      if (reg == "PEN") ylim <- c(-3,3)
    }

    xlim <- ylim

    pdf(file = paste0( e$outdir,"/MAIN_scen_dep_", reg, ".pdf"), width = 9, height = 8)
    par(mar = c(5, 5, 1, 1))

    plot(1:2, 1:2, type = "n", main = NULL, xaxs = "i", yaxs = "i",
         xlim = xlim, xlab = "Sea level contribution at 2100 for RCP2.6/SSP1-26 (cm SLE)",
         ylim = ylim, ylab = "Sea level contribution at 2100 for RCP8.5/SSP5-85 (cm SLE)",
         cex.main = 1.5, cex.axis = 1.5, cex.lab = 1.6 ) # "Forcing at 2100 (Wm-2)")
    abline( h=0, col = "grey", lwd = 0.5 )
    abline( v=0, col = "grey", lwd = 0.5 )

    col_red <- rgb(238, 0, 0, maxColorValue = 255, alpha = 51)
    col_green <- rgb(34, 139, 34, maxColorValue = 255, alpha = 51)

    # Red triangle
    polygon( c(0, xlim[2], 0), c(0, ylim[2],ylim[2]), col = col_red, border = NA)
    if (reg == "ALL") xt <- xlim[2]*0.5
    else xt <- xlim[2]*0.4
    text( xt, ylim[2] - 0.08*(ylim[2] - ylim[1]),
          "MASS LOSS DOMINATES\nMore mass loss\nunder high emissions", cex = 1.5)

    # Green triangle
    polygon( c(0, xlim[1], 0), c(0, ylim[1],ylim[1]), col = col_green, border = NA )
    if (reg != "ALL") text( xlim[1]  + 0.3*(xlim[2]-xlim[1]), ylim[1] + 0.08*(ylim[2]-ylim[1]),
                            "ACCUMULATION DOMINATES\nMore mass gain\nunder high emissions", cex = 1.5 )

    # Other labels
    if (reg != "ALL") {
      text( xlim[2]*0.6, ylim[1] + 0.6*(ylim[2]-ylim[1]),
            "More mass loss\nunder low emissions", cex = 1.5 )
      text( xlim[2]*0.5, ylim[1] + 0.25*(ylim[2]-ylim[1]),
            "Mass loss under low emissions,\nmass gain under high emissions", cex = 1.5 )
    }

    col1 <- rgb(6,104,139, maxColorValue = 255)
    col2 <- rgb(243,122,107, maxColorValue = 255)
    col3 <- rgb(119,12,46, maxColorValue = 255) # rgb(254,206,95,

    # Region name

    ypos <- ylim[2] - 0.07*(ylim[2] - ylim[1])
    xpos <- xlim[1] + 0.02*(xlim[2] - xlim[1])

    text(xpos, ypos, pos = 4, reg_name, cex = 2)

    xpos <- xpos + 0.02*(xlim[2] - xlim[1])
    ypos <- ypos - 0.07*(ylim[2] - ylim[1])

    if (is == "GrIS") {
      points(xpos, ypos, pch = 19, col = col1, cex = 1.7)
      text(xpos + 0.01*(xlim[2] - xlim[1]), ypos, pos = 4, "MIROC5", cex = 1.5)
      ypos <- ypos - 0.05*(ylim[2] - ylim[1])
      points(xpos, ypos, pch = 19, col = col3, cex = 1.7)
      text(xpos + 0.01*(xlim[2] - xlim[1]), ypos, pos = 4, "CNRM-CM6-1", cex = 1.5)

    }
    if (is == "AIS") {
      points( xpos, ypos, pch = 19, col = col1, cex = 1.7) # col = rgb(136, 20, 72, maxColorValue = 255, alpha = 153))
      text( xpos + 0.01*(xlim[2] - xlim[1]), ypos, pos = 4, "NorESM1-M", cex = 1.5)

      ypos <- ypos - 0.05*(ylim[2] - ylim[1])
      points(xpos, ypos, pch = 19, col = col3, cex = 1.7) #col = rgb(56, 240, 172, maxColorValue = 255, alpha = 153))
      text(xpos + 0.01*(xlim[2] - xlim[1]), ypos, pos = 4, "CNRM-CM6-1", cex = 1.5)

      ypos <- ypos - 0.05*(ylim[2] - ylim[1])
      points(xpos, ypos, pch = 19, col = col2, cex = 1.7) # rgb(57, 146, 131, maxColorValue = 255, alpha = 153))
      text( xpos + 0.01*(xlim[2] - xlim[1]), ypos, pos = 4, "IPSL-CM5A-MR",cex = 1.5)

    }

    ypos <- ypos - 0.07*(ylim[2] - ylim[1])
    points(xpos, ypos, pch = 1, cex = 1.7)
    ypos <- ypos - 0.03*(ylim[2] - ylim[1])
    text( xpos + 0.01*(xlim[2] - xlim[1]), ypos, pos = 4, "Not using\nparameterisation", cex = 1.5)

    if (reg == "ALL") subfig <- "a"
    if (reg == "WAIS") subfig <- "b"
    if (reg == "EAIS") subfig <- "c"
    if (reg == "PEN") subfig <- "d"

    text( xlim[1] + 0.92*(xlim[2] - xlim[1]),
          ylim[1] + 0.07*(ylim[2] - ylim[1]),
          pos = 4, font = 2, subfig, cex = 2.2)


    print("", quote = FALSE)
    print("Plot data:")
    print("", quote = FALSE)

    for (ff in 1:dim(to_plot)[1]) {

      if ( to_plot[ff, "scenario"] %in% c("RCP26", "SSP126") ) {

        # Get unique run identifiers
        run_tag <- unlist(select(to_plot[ff,], "GCM", "group", "model", "melt", "collapse"))

        # Get both scenarios
        fd <- to_plot %>%
          filter(GCM == run_tag[1] & group == run_tag[2] & model == run_tag[3] &
                   abs(melt - as.numeric(run_tag[4])) < 1e-6 )
        if (is == "AIS") fd <- filter(fd, collapse == run_tag[5])

        # If find two
        if(dim(fd)[1] == 2) {

          if (run_tag[1] == "CNRM-CM6-1") col_gcm <- col3
          if (run_tag[1] == "MIROC5") col_gcm <- col1
          if (run_tag[1] == "IPSL-CM5A-MR") col_gcm <- col2
          if (run_tag[1] == "NorESM1-M") col_gcm <- col1

          # Standard melt
          pch_melt <- 19 # filled circle

          # Open melt
          if ( abs(as.numeric(run_tag[4]) - e$open_melt[[is]]) < 1e-6 ) pch_melt <- 1 # open circle

          points( fd[fd$scenario %in% c("RCP26", "SSP126"), "y2100"],
                  fd[fd$scenario %in% c("RCP85", "SSP585"), "y2100"],
                  pch = pch_melt, col = col_gcm, cex = 1.7 )
          print(paste(unlist(run_tag), collapse = " "))
        }
      }
    }

    dev.off()

  } # if ice sheet


  if (is == "AIS") {

    print("", quote = FALSE)
    print("", quote = FALSE)
    print("____________________________________________", quote = FALSE)
    print("Antarctic gamma-dependence by GCM plots:", quote = FALSE)

    # ED: gamma dep --------------------------------------

    # Melt-dependence of 4 models that ran all
    # Separate plot for each GCM
    # Vertical dashed is imputed melt

    # SICOPOLIS MeanAnt_low OK PEN?
    # Note JPL1 ISSM sparse
    # Note NCAR CISM MIROC PIG_med nearly identical to SICOPOLIS, but not quite
    # And indicates open is lower half / middle of PIG -> suggests tuned to observations

    ylim <- c(-10, 30)
    if (reg == "PEN") ylim <- c(-0.5, 1.5)

    for (gcm in c("NorESM1-M", "MIROC-ESM-CHEM", "CCSM4")) {

      pdf(file = paste0( e$outdir,"/ED_gamma_dep_", gcm, "_", reg, ".pdf"), width = 9, height = 8)
      par(mar = c(5, 5, 1.5, 1))

      xlab <- expression(paste("Antarctic basal melt parameter, ", gamma, " (x",10^3," m ",a^-1,")"))

      plot(1:3, 1:3, xlim = c(0, 500000)/e$sc, ylim = ylim, yaxs = "i", type = "n", main = NULL,
           xlab = xlab, ylab = "Sea level contribution at 2100 (cm SLE)",
           cex.main = 1.5, cex.axis = 1.5, cex.lab = 1.8 )
      abline(h = 0)
      abline(v = e$open_melt[["AIS"]], lty = 5, lwd = 0.5)

      # Select RCP8.5, collapse off and GCM
      fd <- to_plot %>%
        filter( scenario == "RCP85" & GCM == gcm & collapse == 0 )

      # Plot each ice sheet model and add to legend
      fd_model <- filter( fd, group == "ILTS_PIK" & model == "SICOPOLIS")
      points( select(fd_model, "melt", "y2100"), pch = 16, col = e$scen_col_trans[["85"]], cex = 2 )
      print("", quote = FALSE)
      print(select(fd_model, -melt0, -resolution, -y2015, -y2099))
      points(10, 0.75*ylim[2], pch = 16, col = e$scen_col_trans[["85"]], cex = 2)
      text(14, 0.75*ylim[2], pos = 4,
           paste(unlist(unique(select(fd_model, group))), unlist(unique(select(fd_model, model))), sep="/"), cex = 1.7)

      fd_model <- filter( fd, group == "NCAR" & model == "CISM")
      points( select(fd_model, "melt", "y2100"), pch = 3, col = e$scen_col_trans[["85"]], cex = 2.2)
      print("", quote = FALSE)
      print(select(fd_model, -melt0, -resolution, -y2015, -y2099))
      points(10, 0.68*ylim[2], pch = 3, col = e$scen_col_trans[["85"]], cex = 2)
      text(14, 0.68*ylim[2], pos = 4,
           paste(unlist(unique(select(fd_model, group))), unlist(unique(select(fd_model, model))), sep="/"), cex = 1.7)

      fd_model <- filter( fd, group == "LSCE" & model == "GRISLI")
      points( select(fd_model, "melt", "y2100"), pch = 8, col = e$scen_col_trans[["85"]], cex = 2)
      print("", quote = FALSE)
      print(select(fd_model, -melt0, -resolution, -y2015, -y2099))
      points(10, 0.61*ylim[2], pch = 8, col = e$scen_col_trans[["85"]], cex = 2)
      text(14, 0.61*ylim[2], pos = 4,
           paste(unlist(unique(select(fd_model, group))), unlist(unique(select(fd_model, model))), sep="/"), cex = 1.7)

      fd_model <- filter( fd, group == "JPL1" & model == "ISSM")
      points( select(fd_model, "melt", "y2100"), pch = 1, col = e$scen_col_trans[["85"]], cex = 2)
      print("", quote = FALSE)
      print(select(fd_model, -melt0, -resolution, -y2015, -y2099))
      points(10, 0.54*ylim[2], pch = 1, col = e$scen_col_trans[["85"]], cex = 2)
      text(14, 0.54*ylim[2], pos = 4,
           paste(unlist(unique(select(fd_model, group))), unlist(unique(select(fd_model, model))), sep="/"), cex = 1.7)

      fd_model <- filter( fd, group == "CPOM" & model == "BISICLES")
      points( select(fd_model, "melt", "y2100"), pch = 17, col = e$scen_col_trans[["85"]], cex = 2)
      print("", quote = FALSE)
      print(select(fd_model, -melt0, -resolution, -y2015, -y2099))
      points(10, 0.47*ylim[2], pch = 17, col = e$scen_col_trans[["85"]], cex = 2)
      text(14, 0.47*ylim[2], pos = 4,
           paste(unlist(unique(select(fd_model, group))), unlist(unique(select(fd_model, model))), sep="/"), cex = 1.7)

      # Legend
      text( 0, 0.92*ylim[2], pos = 4, reg_name, cex = 2)
      text( 0, 0.85*ylim[2], pos = 4, paste("Climate model:", gcm), cex = 1.7)


      if (gcm == "NorESM1-M") subfig <- ifelse(reg == "WAIS", "a", "b")
      if (gcm == "MIROC-ESM-CHEM") subfig <- ifelse(reg == "WAIS", "c", "d")
      if (gcm == "CCSM4") subfig <- ifelse(reg == "WAIS", "e", "f")

      text( 480, -8, pos = 4, font = 2, subfig, cex = 2.2)

      dev.off()

    }


    # gamma-RCP  --------------------------------------

    print("", quote = FALSE)
    print("", quote = FALSE)
    print("____________________________________________", quote = FALSE)
    print("Gamma-RCP dependence", quote = FALSE)
    print("", quote = FALSE)

    if (reg == "WAIS") ylim <- c(-2, 20)
    if (reg == "EAIS") ylim <- c(-5, 10)
    if (reg == "PEN") ylim <- c(-0.5, 1.5)

    pdf(file = paste0( e$outdir,"/ED_gamma_RCP_", reg, ".pdf"), width = 9, height = 8)
    par(mar = c(5, 5, 2, 2))

    plot(1:3, 1:3, xlim = c(0, 500000)/e$sc, ylim = ylim, yaxs = "i", type = "n", main = NULL,
         xlab = expression(paste("Antarctic basal melt parameter, ", gamma[0], " (x ",10^-3,")")),
         ylab = "Sea level contribution at 2100 (cm SLE)",
         cex.main = 1.5, cex.axis = 1.5, cex.lab = 1.5 )
    abline(h = 0)
    text( 0, 0.92*ylim[2], pos = 4, "Effect of melt by RCP for NorESM1-M", cex = 1.2)
    text( 0, 0.85*ylim[2], pos = 4, reg_name, cex = 1.5)

    points( c(10, 20), rep(0.77*ylim[2], 2), pch = 24, cex = 1.1,
            col = c(e$scen_col[["85"]], e$scen_col[["26"]]), bg = c(e$scen_col_trans[["85"]], e$scen_col_trans[["26"]]))
    text( 23, 0.77*ylim[2], pos = 4, "BISICLES: RCP8.5, RCP2.6", cex = 1.2)
    points( c(10, 20), rep(0.72*ylim[2], 2), pch = 21, cex = 1.2,
            col = c(e$scen_col[["85"]], e$scen_col[["26"]]), bg = c(e$scen_col_trans[["85"]], e$scen_col_trans[["26"]]))
    text( 23, 0.72*ylim[2], pos = 4, "ISSM: RCP8.5, RCP2.6", cex = 1.2)

    # Get GCM and melt values with both RCPs present for both models
    fd <- filter( to_plot, GCM == "NorESM1-M",
                  abs(melt - e$melt_values[["AIS"]]["med"]) < 1e-6 |
                    abs(melt - e$melt_values[["AIS"]]["PIGL_med"]) < 1e-6 |
                    abs(melt - e$melt_values[["AIS"]]["PIGL_high"]) < 1e-6  )

    for (mod in c("ISSM", "BISICLES"))  {

      print(mod)
      if (mod == "ISSM") {
        fd_model <- filter(fd, group == "JPL1" & model == "ISSM")
        pch = 21
      }
      if (mod == "BISICLES") {
        fd_model <- filter(fd, group == "CPOM" & model == "BISICLES")
        ypos <- 0.7*ylim[2]
        pch = 24
      }

      print(select(fd_model, -melt0, -resolution, -y2015))

      diff1 <- select(filter(fd_model, scenario == "RCP85",
                             abs(melt - e$melt_values[["AIS"]]["med"]) < 1e-6 ), "y2100") -
        select(filter(fd_model, scenario == "RCP85",
                      abs(melt - e$melt_values[["AIS"]]["PIGL_high"]) < 1e-6 ), "y2100")

      diff2 <- select(filter(fd_model, scenario == "RCP26",
                             abs(melt - e$melt_values[["AIS"]]["med"]) < 1e-6 ), "y2100") -
        select(filter(fd_model, scenario == "RCP26",
                      abs(melt - e$melt_values[["AIS"]]["PIGL_high"]) < 1e-6 ), "y2100")

      print( paste("RCP26/RCP85 ratio:", signif(diff2/diff1, 2) ) )

      points( select(filter(fd_model, scenario == "RCP85"), melt, y2100), pch = pch, cex = 1.5,
              col = e$scen_col[["85"]], bg = e$scen_col_trans[["85"]])
      points( select(filter(fd_model, scenario == "RCP26"), melt, y2100), pch = pch, cex = 1.5,
              col = e$scen_col[["26"]], bg = e$scen_col_trans[["26"]])

    }

    dev.off()

    # shelf-gamma for RCP8.5 --------------------------------------

    print("", quote = FALSE)
    print("_____________________________", quote = FALSE)
    print("Shelf-gamma interactions")
    print("", quote = FALSE)

    # No need to plot - just 4 points
    #pdf(file = paste0( e$outdir,"/ED_collapse_gamma_", reg, ".pdf"), width = 9, height = 8)
    par(mar = c(5, 5, 2, 2))

    if (reg == "WAIS") ylim <- c(0, 3)
    if (reg == "EAIS") ylim <- c(0, 3)
    if (reg == "PEN") ylim <- c(0, 3)

    plot(1:3, 1:3, xlim = c(0, 500000)/e$sc, ylim = ylim, yaxs = "i", type = "n",
         main = NULL,
         xlab = expression(paste("Antarctic basal melt parameter, ", gamma[0], " (x ",10^-3,")")),
         ylab = "Additional sea level contribution at 2100 (cm SLE)", cex.main = 1.5, cex.axis = 1.5, cex.lab = 1.5 )
    abline(h = 0)

    # Legend
    text(0, 0.9*ylim[2], pos = 4, paste("Effect of ice shelf collapse under RCP8.5"), cex = 1.2)
    text(0, 0.85*ylim[2], pos = 4, reg_name, cex = 1.2)

    fd <- to_plot %>%
      filter( scenario == "RCP85" & GCM == "CCSM4")

    for (mod in c("ISSM", "BISICLES"))  {

      print("", quote = FALSE)
      print(mod, quote = FALSE)
      if (mod == "ISSM") {
        fd_model <- filter(fd, group == "JPL1" & model == "ISSM")
        ypos <- 0.75*ylim[2]
        pch = 19
      }
      if (mod == "BISICLES") {
        fd_model <- filter(fd, group == "CPOM" & model == "BISICLES")
        ypos <- 0.7*ylim[2]
        pch = 17
      }

      melt <- unlist(select(filter(fd_model, collapse == 1), "melt"))

      for (mv in melt) {

        fd_mv <- filter(fd_model, abs(melt - mv) < 1e-3)

        if (dim(fd_mv)[1] == 2) {
          diff <- select(filter(fd_model, abs(melt - mv) < 1e-6, collapse == 1), "y2100") - select(filter(fd_model, abs(melt - mv) < 1e-6, collapse == 0), "y2100")
          print(sprintf("Melt: %.2f. Additional contribution: %.2f", mv, diff))
          points( mv, diff, pch = pch, col = e$scen_col[["85"]], cex = 1.5)
        }

      }

      points(15, ypos, pch = pch, col = e$scen_col[["85"]], cex = 1.2)
      text(15, ypos, pos = 4, paste(mod))

    }

    #dev.off()

    # ED: collapse by GCM --------------------------------------

    pdf(file = paste0( e$outdir,"/ED_collapse_GCM_", reg, ".pdf"), width = 9, height = 10)
    par(mar = c(14, 5, 3, 2))

    if (reg == "WAIS") ylim <- c(-1, 2)
    if (reg == "EAIS") ylim <- c(-1, 8)
    if (reg == "PEN") ylim <- c(-1, 6)

    gcm_list <- c("IPSL-CM5A-MR", "MIROC-ESM-CHEM", "HadGEM2-ES", "NorESM1-M", "CCSM4", "CSIRO-Mk3-6-0")

    plot(1:3, 1:3, xlim = c(0.5,6.5), ylim = ylim, yaxs = "i", type = "n", main = NULL,
         xaxt = "n", xlab ="", ylab = "Additional sea level contribution at 2100 (cm SLE)",
         cex.main = 1.5, cex.axis = 1.5, cex.lab = 1.8)
    axis( side = 1, labels = gcm_list, at = 1:6, las = 2, cex.axis = 1.6,
          xlab = "Climate model forcing")
    abline(h = 0)
    text(0.5, ylim[1] + 0.92*(ylim[2]-ylim[1]), pos = 4, reg_name, cex = 2)
    text(0.5, ylim[1] + 0.85*(ylim[2]-ylim[1]), pos = 4, "Effect of ice shelf collapse by climate model", cex = 1.8)

    fd <- to_plot %>%
      filter( scenario == "RCP85" & abs(melt - e$melt_values[["AIS"]]["med"]) < 0.01 )

    # Approx 2cm increase with shelf collapse
    gcm_ind <- 0
    finish_leg <- FALSE

    for (gcm in gcm_list) {

      fd_model <- filter(fd, GCM == gcm)
      gcm_ind <- gcm_ind + 1

      model_list <- c( "SICOPOLIS", "GRISLI", # These are the two that did all GCMs
                       "PISM1", "MALI", "IMAUICE1", "IMAUICE2", "JPL1_ISSM", "UCIJPL_ISSM",
                       "fETISh_16km", "fETISh_32km")

      pch_cex <- 2

      for (mm in model_list) {

        fd_model2 <- filter(fd_model, model == mm)
        if (mm == "JPL1_ISSM") fd_model2 <- filter(fd_model, group == "JPL1" & model == "ISSM")
        if (mm == "UCIJPL_ISSM") fd_model2 <- filter(fd_model, group == "UCIJPL" & model == "ISSM")

        pch_ind <- 19
        col <- e$scen_col[["85"]]

        if (mm == "SICOPOLIS") {
          yleg <- 0.75
          pch_ind <- 24
        }
        if (mm == "GRISLI") {
          pch_ind <- 22
          yleg <- 0.69
        }
        if ( ! mm %in% c("SICOPOLIS", "GRISLI")) col <- e$scen_col_trans2[["85"]]

        if ( mm %in% c("SICOPOLIS", "GRISLI")) {
          points( 0.7, ylim[1] + yleg*(ylim[2]-ylim[1]), pch = pch_ind, cex = pch_cex, col = col, bg = col )
          text( 0.73, ylim[1] + yleg*(ylim[2]-ylim[1]), pos = 4, paste(unlist(select(fd_model2, group))[1], mm, sep = "/"), cex = 1.5 )
        } else {
          if ( ! finish_leg) {
            yleg <- 0.63
            points( 0.7, ylim[1] + yleg*(ylim[2]-ylim[1]), pch = pch_ind, cex = pch_cex, col = col, bg = col )
            text( 0.73, ylim[1] + yleg*(ylim[2]-ylim[1]), pos = 4, "Other models", cex = 1.5 )
            finish_leg <- TRUE
          }
        }

        if (dim(fd_model2)[1] == 2) {
          diff <- select(filter(fd_model2, collapse == 1), "y2100") - select(filter(fd_model2, collapse == 0),  "y2100")
          # print(diff)
          points( gcm_ind, diff, pch = pch_ind, cex = pch_cex, col = col, bg = col)
        }
      }

    }

    if (reg == "WAIS") subfig <- "a"
    if (reg == "EAIS") subfig <- "b"
    if (reg == "PEN") subfig <- "c"
    text( 6.2, ylim[1] + 0.92*(ylim[2] - ylim[1]), pos = 4, font = 2, subfig, cex = 2.2)

    dev.off()

    #      points( select(filter(fd_model, collapse == 1), "collapse", "y2100"), pch = 8, col = e$scen_col_trans[["85"]])
    #     points( select(filter(fd_model, collapse == 0), "collapse", "y2100"), pch = 3, col = e$scen_col_trans[["85"]])
    #    lines( x = unlist(select(filter(fd_model, GCM == "NorESM1-M"), "collapse")),
    #          y = unlist(select(filter(fd_model, GCM == "NorESM1-M"), "y2100") ))
    #  points( select(filter(fd_model, collapse == 1 & GCM == "NorESM1-M"), "collapse", "y2100"), pch = 8, col = "black")
    # points( select(filter(fd_model, collapse == 0 & GCM == "NorESM1-M"), "collapse", "y2100"), pch = 3, col = "black")

    # Little change/small decrease
    #fd_model <- filter(fd, model == "GRISLI")
    #points( select(filter(fd_model, collapse == 1), "collapse", "y2100"), pch = 19, col = e$scen_col_trans[["85"]])
    #points( select(filter(fd_model, collapse == 0), "collapse", "y2100"), pch = 1, col = e$scen_col_trans[["85"]])
    #lines( x = unlist(select(filter(fd_model, GCM == "NorESM1-M"), "collapse")),
    #      y = unlist(select(filter(fd_model, GCM == "NorESM1-M"), "y2100") ))
    #  points( select(filter(fd_model, collapse == 1 & GCM == "NorESM1-M"), "collapse", "y2100"), pch = 19, col = "black")
    #  points( select(filter(fd_model, collapse == 0 & GCM == "NorESM1-M"), "collapse", "y2100"), pch = 1, col = "black")



  } # AIS
}


