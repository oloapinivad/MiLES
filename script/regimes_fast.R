######################################################
#------Regimes routines computation for MiLES--------#
#-------------P. Davini (May 2017)-------------------#
######################################################

miles.regimes.fast <- function(project, dataset, expid, ens, year1, year2, season, z500filename, FILESDIR, nclusters = nclusters, doforce) {

  # t0
  t0 <- proc.time()

  # region boundaries for North Atlantic
  if (nclusters != 4 | season != "DJF") {
    stop("Beta version: unsupported season and/or number of clusters")
  }

  # test function to smooth seasonal cycle: it does not work fine yet, keep it false
  smoothing <- T
  xlim <- c(-80, 40)
  ylim <- c(30, 87.5)

  # define file where save data
  savefile1 <- file.builder(FILESDIR, "Regimes", "RegimesPattern", project, dataset, expid, ens, year1, year2, season)

  # check if data is already there to avoid re-run
  if (file.exists(savefile1)) {
    print("Actually requested weather regimes data is already there!")
    print(savefile1)
    if (as.logical(doforce)) {
      print("Running with doforce=true... re-run!")
    } else {
      print("Skipping... activate doforce=true if you want to re-run it")
      q()
    }
  }

  # setting up time domain
  years <- year1:year2
  timeseason <- season2timeseason(season)

  # decide if we want to include this
  # increase the number of files to load for edge-day smoothing
  # if (smoothing) {
  # 	timeseason0=timeseason
  # 	if (season=="DJF") {
  # 		timeseason=sort(c(timeseason,11,3))
  # 	} else {
  # 	timeseason=sort(c(timeseason[1]-1,timeseason,timeseason[length(timeseason)]+1))
  # 	}
  # }

  # file opening
  fieldlist <- ncdf.opener.universal(z500filename, namevar = "zg", tmonths = timeseason, tyears = years, rotate = "full", exportlonlat = F)
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

  print("Compute anomalies based on daily mean")
  # smoothing flag and daily anomalies
  if (smoothing) {
    Z500anom <- daily.anom.run.mean5(ics, ipsilon, Z500, etime)
  } else {
    Z500anom <- daily.anom.mean(ics, ipsilon, Z500, etime)
  }

  # compute weather regimes: new regimes2 function with minimum variance evaluation
  weather_regimes <- regimes2(ics, ipsilon, Z500anom, ncluster = nclusters, ntime = 1000, minvar = 0.8, xlim, ylim, alg = "Hartigan-Wong")

  # Cluster assignation: based on the position of the absolute maximum/minimum
  # negative value for NAO-, maximum for the other 3 regimes
  compose <- weather_regimes$regimes
  names <- paste("Regimes", 1:nclusters)
  position <- rbind(c(-45, 65), c(-35, 50), c(10, 60), c(-20, 60))
  rownames(position) <- c("NAO-", "Atlantic Ridge", "Scandinavian Blocking", "NAO+")

  # minimum distance in degrees to assign a regime name
  min_dist_in_deg <- 20

  # loop
  for (i in 1:nclusters) {

    # find position of max and minimum values
    MM <- which(compose[, , i] == max(compose[, , i], na.rm = T), arr.ind = T)
    mm <- which(compose[, , i] == min(compose[, , i], na.rm = T), arr.ind = T)

    # use maximum or minimum (use special vector to alterate distance when needed)
    if (max(compose[, , i], na.rm = T) > abs(min(compose[, , i], na.rm = T))) {
      distmatrix <- rbind(c(ics[MM[1]], ipsilon[MM[2]]), position + c(0, 0, 0, 1000))
    } else {
      distmatrix <- rbind(c(ics[mm[1]], ipsilon[mm[2]]), position + c(1000, 1000, 1000, 0))
    }

    # compute distances and names assignation
    distMM <- dist(distmatrix)[1:nclusters]
    print(distMM)

    # minimum distance for correct assignation of 15 deg
    if (min(distMM) < min_dist_in_deg) {
      names[i] <- rownames(position)[which.min(distMM)]

      # avoid double assignation
      if (i > 1 & any(names[i] == names[1:max(c(1, i - 1))])) {
        print("Warning: double assignation of the same regime. Avoiding last assignation...")
        names[i] <- paste("Regime", i)
      }
    }
    print(names[i])
  }

  t1 <- proc.time() - t0
  print(t1)

  ##########################################################
  #------------------------Save to NetCDF------------------#
  ##########################################################

  # saving output to netcdf files
  print("saving NetCDF climatologies...")

  # dimension definition (using default 1850-01-01 reftime)
  dims <- ncdf.defdims(ics, ipsilon, timeaxis)

  # extra dimensions definition
  cl <- ncdim_def("clust", units = "" , 1:nclusters, longname = "Cluster index")

  # var definition
  unit <- "m"
  longvar <- "Weather Regimes Pattern"
  pattern_ncdf <- ncvar_def("Regimes", unit, list(dims$x, dims$y, cl), -999, longname = longvar, prec = "single", compression = 1)

  unit <- paste0("0-", nclusters)
  longvar <- "Weather Regimes Cluster Index"
  cluster_ncdf <- ncvar_def("Indices", unit, list(dims$t), -999, longname = longvar, prec = "single", compression = 1)

  unit <- "%"
  longvar <- "Weather Regimes Frequencies"
  frequencies_ncdf <- ncvar_def("Frequencies", unit, list(cl), -999, longname = longvar, prec = "single", compression = 1)

  # testnames
  dimnchar <- ncdim_def("nchar", "", 1:max(nchar(names)), create_dimvar = FALSE)
  names_ncdf <- ncvar_def("Names", "", list(dimnchar, cl), prec = "char")

  # saving files (new list based method with ncdf.writer())
  nc_var <- list(pattern_ncdf, cluster_ncdf, frequencies_ncdf, names_ncdf)
  nc_field <- list(weather_regimes$regimes, weather_regimes$cluster, weather_regimes$frequencies, names)
  names(nc_field) <- names(nc_var) <- c("Regimes", "Indices", "Frequencies", "Names")
  ncdf.writer(savefile1, nc_var, nc_field)

}

# blank line
cat("\n\n\n")

# REAL EXECUTION OF THE SCRIPT
# read command line
args <- commandArgs(TRUE)

# number of required arguments from command line
name_args <- c("project", "dataset", "expid", "ens", "year1", "year2", "season", "z500filename", "FILESDIR", "PROGDIR", "nclusters", "doforce")

# if there arguments, check them required args and assign
if (length(args) != 0) {
  req_args <- length(name_args)
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
    miles.regimes.fast(project, dataset, expid, ens, year1, year2, season, z500filename, FILESDIR, nclusters, doforce)
  }
}
