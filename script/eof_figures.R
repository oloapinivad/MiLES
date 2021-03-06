######################################################
#--------Routines for EOFs plotting for MiLES--------#
#-------------P. Davini (May 2017)-------------------#
######################################################

# DECLARING THE FUNCTION: EXECUTION IS AT THE BOTTOM OF THE SCRIPT

miles.eof.figures <- function(project, dataset, expid, ens, year1, year2,
                              project_ref, dataset_ref, expid_ref, ens_ref, year1_ref, year2_ref,
                              season, FIGDIR, FILESDIR, REFDIR, CFGSCRIPT, tele) {

  # R configuration file
  source(CFGSCRIPT)

  # use filebuilding script to access to file
  nomefile_exp <- file.builder(FILESDIR, paste0("EOFs/", tele), "EOFs", project, dataset, expid, ens, year1, year2, season)

  # check for REFDIR==FILESDIR, i.e. if we are using the climatology provided by MiLES or another dataset MiLES-generated
  if (REFDIR != FILESDIR) {
    nomefile_ref <- paste0(file.path(REFDIR, paste0("EOFs/", tele)), "/EOFs_", dataset_ref, "_", year1_ref, "_", year2_ref, "_", season, ".nc")
  } else {
    # use file.builder to create the path of the blocking files
    nomefile_ref <- file.builder(FILESDIR, paste0("EOFs/", tele), "EOFs", project_ref, dataset_ref, expid_ref, ens_ref, year1_ref, year2_ref, season)
  }

  # EOFs to plot (depends on how many computed by CDO!)
  neofs <- 4

  ##########################################################
  #-----------------Loading datasets-----------------------#
  ##########################################################

  # loading anomalies and variances of experiment
  variance_exp <- ncdf.opener(nomefile_exp, namevar = "Variances", rotate = "no") * 100 # convert to percentage
  regressions_exp <- ncdf.opener(nomefile_exp, namevar = "Regressions", rotate = "no")

  # loading reference field
  variance_ref <- ncdf.opener(nomefile_ref, namevar = "Variances", rotate = "no") * 100 # convert to percentage
  regressions_ref <- ncdf.opener(nomefile_ref, namevar = "Regressions", rotate = "no")


  ##########################################################
  #-----------------Produce figures------------------------#
  ##########################################################

  # plot properties
  info_exp <- info.builder(dataset, expid, ens, year1, year2, season)
  info_ref <- info.builder(dataset_ref, expid_ref, ens_ref, year1_ref, year2_ref, season)
  lev_field <- seq(-150, 150, 20)
  lev_diff <- seq(-95, 95, 10)

  # loop on number of EOFs
  for (neof in 1:neofs) {
    linear_exp <- regressions_exp[, , neof]
    linear_ref <- regressions_ref[, , neof]

    # check and flip signs (to be in agreement with reference field)
    if (cor(c(linear_ref), c(linear_exp)) < 0) {
      linear_exp <- (-linear_exp)
    }

    #-----plotting-------#

    # plot properties
    region <- tele # if it is a box of lonlat
    if (tele == "NAO") {
      region <- "North Atlantic"
    }
    if (tele == "AO") {
      region <- "Northern Hemisphere"
    }
    if (tele == "PNA") {
      region <- "North Pacific"
    }
    title_name <- paste0(region, " EOF", neof)

    # define figure
    figname <- fig.builder(FIGDIR, paste0("EOFs/", tele), paste0("EOF", neof), project, dataset, expid, ens, year1, year2, season, output_file_type)
    print(figname)

    # Chose output format for figure - by JvH
    open.plot.device(figname, output_file_type, CFGSCRIPT)

    # where to plot variances values
    if (map_projection == "no") {
      varpoints <- c(120, 85)
    } else {
      varpoints <- c(0, 0.7)
    }

    # plot properties
    par(plotpar)

    im <- plot.prepare(ics, ipsilon, linear_exp, proj = map_projection, lat_lim = lat_lim)
    filled.contour3(im$x, im$y, im$z, xlab = im$xlab, ylab = im$ylab, main = paste(info_exp), levels = lev_field, color.palette = palette3, xlim = im$xlim, ylim = im$ylim, axes = im$axes)
    mtext(title_name, side = 3, line = .5, outer = TRUE, cex = 2, font = 2)
    proj.addland(ics, ipsilon, proj = map_projection)
    text(varpoints[1], varpoints[2], paste("Variance Explained: ", round(variance_exp[neof], 2), "%", sep = ""), cex = 2)

    im <- plot.prepare(ics, ipsilon, linear_ref, proj = map_projection, lat_lim = lat_lim)
    filled.contour3(im$x, im$y, im$z, xlab = im$xlab, ylab = im$ylab, main = paste(info_ref), levels = lev_field, color.palette = palette3, xlim = im$xlim, ylim = im$ylim, axes = im$axes)
    mtext(title_name, side = 3, line = .5, outer = TRUE, cex = 2, font = 2)
    proj.addland(ics, ipsilon, proj = map_projection)
    image.scale3(volcano, levels = lev_field, color.palette = palette3, colorbar.label = "m", cex.colorbar = imgscl_colorbar, cex.label = imgscl_label, colorbar.width = 1 * af, line.label = imgscl_line)
    text(varpoints[1], varpoints[2], paste("Variance Explained: ", round(variance_ref[neof], 2), "%", sep = ""), cex = 2)

    # delta field plot
    im <- plot.prepare(ics, ipsilon, linear_exp - linear_ref, proj = map_projection, lat_lim = lat_lim)
    filled.contour3(im$x, im$y, im$z, xlab = im$xlab, ylab = im$ylab, main = paste("Difference"), levels = lev_diff, color.palette = palette2, xlim = im$xlim, ylim = im$ylim, axes = im$axes)
    proj.addland(ics, ipsilon, proj = map_projection)
    image.scale3(volcano, levels = lev_diff, color.palette = palette2, colorbar.label = "m", cex.colorbar = imgscl_colorbar, cex.label = imgscl_label, colorbar.width = 1 * af, line.label = imgscl_line)


    dev.off()
  }
}

# blank line
cat("\n\n\n")

# REAL EXECUTION OF THE SCRIPT
# read command line
args <- commandArgs(TRUE)

# number of required arguments from command line
name_args <- c("project", "dataset", "expid", "ens", "year1", "year2", "project_ref", "dataset_ref", "expid_ref", "ens_ref", "year1_ref", "year2_ref", "season", "FIGDIR", "FILESDIR", "REFDIR", "CFGSCRIPT", "PROGDIR", "tele")
req_args <- length(name_args)

# if there arguments, check them required args and assign
if (length(args) != 0) {
  if (length(args) != req_args) {
    # stop if something is wrong
    print(paste(length(args), "arguments received: please specify the following", req_args, "arguments:"))
    print(name_args)
    stop("ERROR!")
  } else {
    # when the number of arguments is ok run the function()
    for (k in 1:req_args) {
      if (args[k] == "") {
        args[k] <- NA
      }
      assign(name_args[k], args[k])
    }
    source(paste0(PROGDIR, "/script/basis_functions.R"))
    miles.eof.figures(
      project, dataset, expid, ens, year1, year2,
      project_ref, dataset_ref, expid_ref, ens_ref, year1_ref, year2_ref,
      season, FIGDIR, FILESDIR, REFDIR, CFGSCRIPT, tele
    )
  }
}
