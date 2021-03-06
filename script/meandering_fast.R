######################################################
#-----Meandering routines computation for MiLES------#
#--------G. Di Capua & P. Davini (Apr 2018)----------#
######################################################
miles.meandering <- function(project, dataset, expid, ens, year1, year2, season, z500filename, FILESDIR, doforce) {

  # t0
  t0 <- proc.time()

  # setting up time domain
  years <- year1:year2
  timeseason <- season2timeseason(season)

  # define folders using file.builder function (takes care of ensembles)
  savefile1 <- file.builder(FILESDIR, "MI", "MIClim", project, dataset, expid, ens, year1, year2, season)
  savefile2 <- file.builder(FILESDIR, "MI", "MIFull", project, dataset, expid, ens, year1, year2, season)

  # check if data is already there to avoid re-run
  if (file.exists(savefile1)) {
    print("Actually requested Meandering Index data is already there!")
    if (as.logical(doforce)) {
      print("Running with doforce=true... re-run!")
    } else {
      print("Skipping... activate doforce=true if you want to re-run it")
      q()
    }
  }

  fieldlist <- ncdf.opener.universal(z500filename, namevar = "zg", tmonths = timeseason, tyears = years, rotate = "full", exportlonlat = F)
  print(str(fieldlist))

  # extract calendar and time unit from the original file
  timeaxis <- fieldlist$time
  ics <- fieldlist$lon
  ipsilon <- fieldlist$lat

  # time array to simplify time filtering
  etime <- power.date.new(timeaxis, verbose = T)
  totdays <- length(timeaxis)

  # declare variable and clean
  Z500 <- fieldlist$field
  rm(fieldlist)

  # smoothing with 5-day running mean
  print("5-day running mean")
  runZ500 <- apply(Z500, c(1, 2), run.mean5)
  runZ500 <- aperm(runZ500, c(2, 3, 1))
  runZ500[is.na(runZ500)] <- Z500[is.na(runZ500)]

  # list of isohypses on which evaluate the MI
  isolvls <- seq(4900, 6200, 5)

  # reference latitude (60N following Di Capua et al., 2016)
  ref_lat <- 60

  # lower and upper latitudinal bound for the maximum isohypse
  lower_bound <- 50
  upper_bound <- 75

  # Running the real code
  MI_lat <- MI_value <- 1:totdays * NA

  # library(rbenchmark)
  # print(benchmark(sapply(isolvls,function(x) {MI.fast(ics,ipsilon,Z500[,,1],x,ref_lat,verbose=F)}),
  # 		vapply(isolvls,function(x) {MI.fast(ics,ipsilon,Z500[,,1],x,ref_lat,verbose=F)},list(1,2))))

  for (t in 1:totdays) {
    progression.bar(t, totdays, each = 20)

    # computed MI
    MI_list <- sapply(isolvls, function(x) {
      MI.fast(ics, ipsilon, runZ500[, , t], x, ref_lat, verbose = F)
    })
    # MI_list=vapply(isolvls,function(x) {MI.fast(ics,ipsilon,Z500[,,1],x,ref_lat,verbose=F)},list(1,2))

    # applyng bounds on latitude and longitudes
    subset <- which(unlist(MI_list[2, ]) >= lower_bound & unlist(MI_list[2, ]) <= upper_bound)
    if (length(subset) > 0) {
      MI_value[t] <- max(unlist(MI_list[1, subset]))
      MI_lat[t] <- unlist(MI_list[2, subset])[which.max(unlist(MI_list[1, subset]))]
    } else {
      MI_value[t] <- NA
      MI_lat[t] <- NA
    }
  }


  tf <- proc.time() - t0
  print(tf)

  ##########################################################
  #------------------------Save to NetCDF------------------#
  ##########################################################

  # saving output to netcdf files
  print("saving NetCDF climatologies...")

  # which fieds to plot/save
  savelist <- c("MI", "MI_lat")
  full_savelist <- c("MI", "MI_lat")

  # dimension definition (using default 1850-01-01 reftime)
  dims <- ncdf.defdims(ics, ipsilon, timeaxis)
  x <- ncdim_def("lon", "degrees_east", 0, longname = "longitude")

  # pre-declare list for loading the variable definition and the fields
  nc_field <- nc_var <- sapply(savelist, function(x) NULL)
  nc_fullfield <- nc_fullvar <- sapply(full_savelist, function(x) NULL)
  for (var in savelist)
  {
    # name of the var
    if (var == "MI") {
      longvar <- "Meandering Index"
      unit <- ""
      field <- mean(MI_value)
      full_field <- MI_value
    }
    if (var == "MI_lat") {
      longvar <- "Meandering Index Latitude"
      unit <- "deg"
      field <- mean(MI_lat)
      full_field <- MI_lat
    }

    # variable definitions
    var_ncdf <- ncvar_def(var, unit, list(x, t = dims$tclim), -999, longname = longvar, prec = "single", compression = 1)
    full_var_ncdf <- ncvar_def(var, unit, list(x, t = dims$t), -999, longname = longvar, prec = "single", compression = 1)

    if (var %in% savelist) {
      nc_var[[which(var == savelist)]] <- var_ncdf
      nc_field[[which(var == savelist)]] <- field
    }

    if (var %in% full_savelist) {
      nc_fullvar[[which(var == full_savelist)]] <- full_var_ncdf
      nc_fullfield[[which(var == full_savelist)]] <- full_field
    }
  }

  # save variables
  # Climatologies Netcdf file creation
  ncdf.writer(savefile1, nc_var, nc_field)
  ncdf.writer(savefile2, nc_fullvar, nc_fullfield)


}

# blank lines
cat("\n\n\n")
# REAL EXECUTION OF THE SCRIPT
# read command line
args <- commandArgs(TRUE)

# number of required arguments from command line
name_args <- c("project", "dataset", "expid", "ens", "year1", "year2", "season", "z500filename", "FILESDIR", "PROGDIR", "doforce")
req_args <- length(name_args)

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
    source(file.path(PROGDIR, "script/meandering_functions.R"))
    miles.meandering(project, dataset, expid, ens, year1, year2, season, z500filename, FILESDIR, doforce)
  }
}
