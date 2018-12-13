######################################################
#-----EOFs routines computation for MiLES--------#
#-------------P. Davini (Feb 2018)-------------------#
######################################################
miles.eofs.fast <- function(project, dataset, expid, ens, year1, year2, season, tele, z500filename, FILESDIR, doforce) {

  # standard defined 4 EOFs
  neofs <- 4

  # t0
  t0 <- proc.time()

  # setting up time domain
  years <- year1:year2
  timeseason <- season2timeseason(season)

  # define folders using file.builder function (takes care of ensembles)
  savefile1 <- file.builder(FILESDIR, paste0("EOFs/", tele), "EOFs", project, dataset, expid, ens, year1, year2, season)

  # select teleconnection region
  if (tele == "NAO") {
    xlim <- c(-90, 40)
    ylim <- c(20, 85)
    rotation <- "full"
  } else if (tele == "AO") {
    xlim <- c(-180, 180)
    ylim <- c(20, 85)
    rotation <- "full"
  } else if (tele == "PNA") {
    xlim <- c(140, 280)
    ylim <- c(20, 85)
    rotation <- "no" # 140E-80W: use trick of rotation for cross-dateline
  } else {
    # use non standard region, detect region with strsplit
    splitter <- as.numeric(strsplit(tele, "_")[[1]])
    if (length(splitter) == 4) {
      xlim <- c(splitter[1], splitter[2])
      ylim <- c(splitter[3], splitter[4])
      if (xlim[2] > 180) {
        rotation <- "no"
      } else {
        rotation <- "full"
      }
    } else {
      stop("Wrong teleconnection region!")
    }
  }

  # check if data is already there to avoid re-run
  if (file.exists(savefile1)) {
    print("Actually requested EOFs data is already there!")
    print(savefile1)
    if (as.logical(doforce)) {
      print("Running with doforce=true... re-run!")
    } else {
      print("Skipping... activate doforce=true if you want to re-run it")
      q()
    }
  }

  # new file opening
  fieldlist <- ncdf.opener.universal(z500filename, namevar = "zg", tmonths = timeseason, tyears = years, rotate = rotation, exportlonlat = F)
  print(str(fieldlist))

  # assign variables
  ics <- fieldlist$lon
  ipsilon <- fieldlist$lat

  # extract calendar and time unit from the original file
  timeaxis <- fieldlist$time

  # declare variable and clean
  Z500 <- fieldlist$field
  rm(fieldlist)

  # time array to simplify time filtering
  etime <- power.date.new(timeaxis, verbose = T)
  totdays <- length(timeaxis)

  # monthly averaging
  print("monthly mean...")

  # new faster monthly mean function
  Z500monthly <- monthly.mean(ics, ipsilon, Z500, etime)

  # climatology
  print("climatological mean...")
  Z500clim <- apply(Z500monthly, c(1, 2), ave, rep(timeseason, length(years)))
  Z500clim <- aperm(Z500clim, c(2, 3, 1))

  # monthly anomalies
  print("anomalies...")
  Z500anom <- Z500monthly - Z500clim

  # compute EOFs
  print("EOFs...")
  EOFS <- eofs(ics, ipsilon, Z500anom, neof = neofs, xlim, ylim, method = "SVD", do_standardize = T, do_regression = T)
  # COEFF=eofs.coeff(ics,ipsilon,Z500anom,EOFS,do_standardize=T) #do we really need this?

  # flip signs of patterns and regressions for NAO and AO
  print("checking signs...")
  for (i in 1:neofs) {
    posreg <- NULL

    # define regions for sign control: boxes where values should be positive
    if (tele == "NAO") {
      if (i == 1) {
        posreg <- c(-30, 30, 40, 50)
      } # NAO
      if (i == 2) {
        posreg <- c(-60, 0, 40, 60)
      } # East Atlantic Pattern
      if (i == 3) {
        posreg <- c(-30, 30, 50, 70)
      } # Scandinavian Blocking
    }

    if (tele == "AO") {
      if (i == 1) {
        posreg <- c(-180, 180, 20, 50)
      } # Arctic Oscillation
      if (i == 2) {
        posreg <- c(-120, -60, 40, 60)
      } # PNA
    }

    # if definition of region exists
    if (!is.null(posreg)) {
      # convert into indices
      xbox <- whicher(EOFS$pattern$x, posreg[1]):whicher(EOFS$pattern$x, posreg[2])
      ybox <- whicher(EOFS$pattern$y, posreg[3]):whicher(EOFS$pattern$y, posreg[4])
      valuereg <- mean(EOFS$pattern$z[xbox, ybox, i])

      # if negative in the box, flip all signs!
      if (valuereg < 0) {
        EOFS$pattern$z[, , i] <- -EOFS$pattern$z[, , i]
        # COEFF[,i]=-COEFF[,i]
        EOFS$regression <- -EOFS$regression
      }
    }
  }

  # expand EOF pattern to save it
  expanded_pattern <- EOFS$regression * NA
  expanded_pattern[whicher(ics, xlim[1]):whicher(ics, xlim[2]), whicher(ipsilon, ylim[1]):whicher(ipsilon, ylim[2]), ] <- EOFS$pattern$z

  t1 <- proc.time() - t0
  print(t1)


  ##########################################################
  #------------------------Save to NetCDF------------------#
  ##########################################################

  # saving output to netcdf files
  print("saving NetCDF climatologies...")
  print(savefile1)

  # dims
  dims <- ncdf.defdims(ics, ipsilon, timeaxis)

  # extra dimensions definition
  ef <- ncdim_def("eof", units = "", 1:neofs, longname = "EOF Number")

  # defining vars
  unit <- "m"
  longvar <- "EOFs Loading Pattern"
  pattern_ncdf <- ncvar_def("Patterns", unit, list(dims$x, dims$y, ef), -999, longname = longvar, prec = "single", compression = 1)

  unit <- "m"
  longvar <- "EOFs Linear Regressions"
  regression_ncdf <- ncvar_def("Regressions", unit, list(dims$x, dims$y, ef), -999, longname = longvar, prec = "single", compression = 1)

  unit <- paste0("0-", neofs)
  longvar <- "PCs timeseries"
  pc_ncdf <- ncvar_def("PCs", unit, list(ef, dims$tmon), -999, longname = longvar, prec = "single", compression = 1)

  unit <- "%"
  longvar <- "EOFs variance"
  variance_ncdf <- ncvar_def("Variances", unit, list(ef), -999, longname = longvar, prec = "single", compression = 1)

  # saving files (new list based method with ncdf.writer())
  nc_var <- list(pattern_ncdf, pc_ncdf, variance_ncdf, regression_ncdf)
  nc_field <- list(expanded_pattern, EOFS$coeff, EOFS$variance, EOFS$regression)
  names(nc_field) <- names(nc_var) <- c("Patterns", "PCs", "Variances", "Regressions")
  ncdf.writer(savefile1, nc_var, nc_field)


}

# blank line
cat("\n\n\n")

# REAL EXECUTION OF THE SCRIPT
# read command line
args <- commandArgs(TRUE)

# number of required arguments from command line
name_args <- c("project", "dataset", "expid", "ens", "year1", "year2", "season", "tele", "z500filename", "FILESDIR", "PROGDIR", "doforce")
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
    source(file.path(PROGDIR, "script/basis_functions.R"))
    miles.eofs.fast(project, dataset, expid, ens, year1, year2, season, tele, z500filename, FILESDIR, doforce)
  }
}
