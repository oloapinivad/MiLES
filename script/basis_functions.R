# basis functions


##########################################################
#------------------------Packages------------------------#
##########################################################

# loadin packages
library("maps")
library("ncdf4")
library("PCICt")

# check if fast linear fit is operative (after R 3.1): 3x faster than lm.fit, 36x faster than lm
if (exists(".lm.fit")) {
  lin.fit <- .lm.fit
} else {
  lin.fit <- lm.fit
}

# check R version as numeric
R_version <- as.numeric(R.Version()$major) + as.numeric(R.Version()$minor) / 10

##########################################################
#-----------------Basic functions------------------------#
##########################################################

# constants
g0 <- 9.80655
Earth.Radius <- 6378137
omega <- 7.29211 * 10^(-5)

# normalize a time series
standardize <- function(timeseries) {
  out <- (timeseries - mean(timeseries, na.rm = T)) / sd(timeseries, na.rm = T)
  return(out)
}


# detect ics ipsilon lat-lon
whicher <- function(axis, number) {
  out <- which.min(abs(axis - number))
  return(out)
}

# produce a 2d matrix of area weight
area.weight <- function(ics, ipsilon, root = T) {
  field <- array(NA, dim = c(length(ics), length(ipsilon)))
  if (root == T) {
    for (j in 1:length(ipsilon)) {
      field[, j] <- sqrt(cos(pi / 180 * ipsilon[j]))
    }
  }

  if (root == F) {
    for (j in 1:length(ipsilon)) {
      field[, j] <- cos(pi / 180 * ipsilon[j])
    }
  }

  return(field)
}

# sector details for blocking extra diagnostics and EOFs sectors
sector.details <- function(ics, ipsilon, SECTOR) {
  if (SECTOR == "Euro") {
    lons <- c(-15, 25)
    lats <- c(50, 65)
    namesec <- "Central Europe"
  }
  if (SECTOR == "WestEuro") {
    lons <- c(-25, 15)
    lats <- c(35, 50)
    namesec <- "Western Europe"
  }
  if (SECTOR == "Azores") {
    lons <- c(-70, -10)
    lats <- c(30, 40)
    namesec <- "Central Atlantic"
  }
  if (SECTOR == "Greenland") {
    lons <- c(-65, -15)
    lats <- c(62.5, 72.5)
    namesec <- "Greenland"
  }
  if (SECTOR == "FullPacific") {
    lons <- c(130, -150)
    lats <- c(60, 75)
    namesec <- "North Pacific"
  }
  if (SECTOR == "FullPacific2") {
    lons <- c(130, 210)
    lats <- c(60, 75)
    namesec <- "North Pacific"
  }

  left1 <- which.min(abs(ics - lons[1]))
  right1 <- which.min(abs(ics - lons[2]))
  low1 <- which.min(abs(ipsilon - lats[1]))
  high1 <- which.min(abs(ipsilon - lats[2]))

  latssel <- low1:high1
  if (SECTOR == "FullPacific") {
    lonssel <- c(left1:length(ics), 1:right1)
  } else {
    lonssel <- left1:right1
  }
  out <- list(lons = lons, lonssel = lonssel, lats = lats, latssel = latssel, name = namesec)
  return(out)
}

# weighted correlation
weighted.cor <- function(x, y, w) {
  w.mean.x <- sum(w * x) / sum(w)
  w.mean.y <- sum(w * y) / sum(w)

  w.cov.xy <- sum(w * (x - w.mean.x) * (y - w.mean.y)) / sum(w)
  w.var.y <- sum(w * (y - w.mean.y) * (y - w.mean.y)) / sum(w)
  w.var.x <- sum(w * (x - w.mean.x) * (x - w.mean.x)) / sum(w)

  corr <- w.cov.xy / sqrt(w.var.x * w.var.y)
  return(corr)
}

# weighted standard deviations
weighted.sd <- function(x, w) {
  w.mean <- sum(w * x) / sum(w)
  v1 <- sum(w)
  v2 <- sum(w^2)
  var <- v1 / (v1^2 - v2) * sum(w * (x - w.mean)^2)
  sdd <- sqrt(var)
  return(sdd)
}


# verbose-only printing function
printv <- function(value, verbosity = TRUE) {
  if (verbosity) {
    print(value)
  }
}


##########################################################
#-----------Array manipulations functions----------------#
##########################################################

# Function to generalize through do.call() n-dimensional array subsetting 
# and array indexing. Derived from Stack Overflow issue
# https://stackoverflow.com/questions/14500707/select-along-one-of-n-dimensions-in-array
array_indexing <- function(field, dim, value, drop = FALSE) {

  # Create list representing arguments supplied to [
  # bquote() creates an object corresponding to a missing argument
  indices <- rep(list(bquote()), length(dim(field)))
  indices[[dim]] <- value

  # do.call on the indices
  out <- do.call("[",c(list(field), indices, list(drop = drop)))

  return(out)
}

# n-dimensional generalized evolution of flipper()
# it reverts the ipsilon a required dimension: by default uses the second dimension
# used by ncdf.opener.universal()
flipper <- function(field, dim = 2) {
  ydim <- dim
  ll <- dim(field)[ydim]
  field <- array_indexing(field, ydim, ll:1)
  return(field)
}

# define rotate along longitudes based on array_indexing (Jan 2019)
# universal to every dimension, 50% faster than previous rotation() functon
# can be used also for ad hoc rotation
# used by ncdf.opener.universal() 
rotation <- function(line, rotate) {

  # dimension of the first dimension (the one to be rotated)
  ll <- dim(line)[1]
  dims <- length(dim(line))

  # default options
  if (is.character(rotate)) {
    if (rotate == "full") {
      # 180 degrees rotation of longitude
      move1 <- 1 / 2 * ll
    } else if (rotate == "half") {
      # 90 degree rotation (useful for TM90)
      move1 <- 1 / 4 * ll
    } else if (rotate == "no") {
      return(line)
    }

  # numeric options
  } else {
    if (rotate == 0 | rotate >= ll) {
      print("Nothing to do")
      return(line)
    } else {
      move1 <- rotate + 1 #always add one to avoid crash
    }
  }

  # move2 as the difference of the number of points  
  move2 <- ll - move1

  # create new elements order
  elements <- c((move2 + 2):ll, 1:(move2 + 1))

  # run the selection using array_indexing
  newline <- array_indexing(line, 1, elements)

  return(newline)
}



##########################################################
#---------String manipulations functions-----------------#
##########################################################

# info string creator
info.builder <- function(dataset, expid, ens, year1, year2, season) {

  # loop on descriptors that are concatenated to create info string
  descriptors <- c(dataset, expid, ens, paste0(year1, "-", year2), season)
  info <- NULL
  for (dcode in descriptors) {
    if (!is.na(dcode)) {
      info <- paste(info, dcode)
    }
  }
  return(info)
}

# basic switch to create NetCDF file names and folders (use recursive structure from v0.6)
file.builder <- function(DATADIR, dir_name, file_name, project, dataset, expid, ens, year1, year2, season) {

  # add directory name descriptor
  DATADIR <- file.path(DATADIR, dir_name)

  if (!is.na(project)) {
    DATADIR <- file.path(DATADIR, project)
  }

  # loop on descriptors that are concatenated to create dir and file name
  descriptors <- c(dataset, expid, ens, paste0(year1, "_", year2), season)
  for (dcode in descriptors) {
    if (!is.na(dcode)) {
      DATADIR <- file.path(DATADIR, dcode)
      file_name <- paste0(file_name, "_", dcode)
    }
  }

  # actually dir.exists is in devtools only for R < 3.2, then is included in base package
  if (exists("dir.exists")) {
    if (!dir.exists(DATADIR)) {
      dir.create(DATADIR, recursive = T)
    }
  } else {
    dir.create(DATADIR, recursive = T, showWarnings = F)
  }
  return(file.path(DATADIR, paste0(file_name, ".nc")))
}
# basic switch to create figures names and folders (use recursive structure from v0.6)
fig.builder <- function(FIGDIR, dir_name, file_name, project, dataset, expid, ens, year1, year2, season, output_file_type) {

  # add directory name descriptor
  FIGDIR <- file.path(FIGDIR, dir_name)

  if (!is.na(project)) {
    FIGDIR <- file.path(FIGDIR, project)
  }

  # loop on descriptors that are concatenated to create dir and file name
  descriptors <- c(dataset, expid, ens, paste0(year1, "_", year2), season)
  for (dcode in descriptors) {
    if (!is.na(dcode)) {
      FIGDIR <- file.path(FIGDIR, dcode)
      file_name <- paste0(file_name, "_", dcode)
    }
  }

  # actually dir.exists is in devtools only for R < 3.2, then is included in base package
  if (exists("dir.exists")) {
    if (!dir.exists(FIGDIR)) {
      dir.create(FIGDIR, recursive = T)
    }
  } else {
    dir.create(FIGDIR, recursive = T, showWarnings = F)
  }

  return(file.path(FIGDIR, paste0(file_name, ".", output_file_type)))
}
2
# progression bar
progression.bar <- function(index, total_length, each = 10) {
  if (any(index == round(seq(0, total_length, , each + 1)))) {
    progression <- paste("--->", round(index / total_length * 100), "%")
    print(progression)
  }
}

##########################################################
#--------------Time Based functions----------------------#
##########################################################

# to convert season charname to months number
season2timeseason <- function(season) {

  # special cases
  if (season == "ALL") {
    timeseason <- 1:12
  } else if (season == "JJA") {
    timeseason <- 6:8
  } else if (season == "DJF") {
    timeseason <- c(1:2, 12)
  } else if (season == "MAM") {
    timeseason <- 3:5
  } else if (season == "SON") {
    timeseason <- 9:11
  } else if (season == "DJFM") {
    timeseason <- c(1:3, 12)
  } else if (season == "JJAS") {
    timeseason <- 6:9
  # otherwise look for strings
  } else {
    charseason <- strsplit(season, "_")[[1]]
    print(charseason)
    if (mean(nchar(charseason)) == 3) {
      timeseason <- which(month.abb %in% charseason)
    } else {
      timeseason <- which(month.name %in% charseason)
    }
  }

  if (length(timeseason) == 0 | min(timeseason) < 0 | max(timeseason) > 13) {
    stop("wrong season selected!")
  }
  return(timeseason)
}

# check number of days for each month
number.days.month <- function(datas, calendar) {

  # Dec-18 update to handle 360 and 365 calendars in file opening
  require(PCICt)
  # evaluate the number of days in a defined month of a year
  datas <- as.PCICt(datas, cal = calendar, format = "%Y-%m-%d")
  m <- format(datas, format = "%m")
  while (format(datas, format = "%m") == m) {
    datas <- datas + 86400
  }
  return(as.integer(format(datas - 1, format = "%d")))
}

# new function to create simple list with date values - Oct-18
# it needs a date or PCICt object, and returns also the season subdivision
power.date.new <- function(datas, verbose = FALSE) {

  # create a "season" for continuous time, used by persistance tracking
  startpoints <- c(0, which(diff(datas) > 1))
  #print(startpoints)
  deltapoints <- diff(c(startpoints, length(datas)))
  seas <- inverse.rle(list(lengths = deltapoints, values = seq(1, length(startpoints))))

  etime <- list(
    day = as.numeric(format(datas, "%d")), month = as.numeric(format(datas, "%m")),
    year = as.numeric(format(datas, "%Y")), data = datas, season = seas
  )

  printv("Time Array Built", verbose)
  printv(paste("Length:", length(seas)), verbose)
  printv(paste("From", datas[1], "to", datas[length(seas)]), verbose)
  return(etime)
}

# create a PCICt object from time, units and calendar type of a NetCDF file (Jan-19)
create.timeline <- function(time, units, caldata = "standard", verbose = F) {

  # assume default
  timeline <- NA

  require(PCICt)

  # If the calendar is missing, assume a standard one
  if (caldata == 0) {
     warning("No calendar found, assuming standard calendar!")
     caldata <- "standard"
  }

  # new method including both absolute and releative time axis (Oct 2018)
  # if "as" is present, this is an absolute time axis
  if (grepl("as", units, fixed = TRUE)) {
    printv("Absolute time axis!", verbose)
    timeline <- as.PCICt(as.character(time), format = "%Y%m%d", cal = caldata)
    # if a "since" is present, this is a relative time axis
  } else if (grepl("since", units, fixed = TRUE)) {
    printv("Relative time axis!", verbose)
    origin <- substr(gsub("[a-zA-Z ]", "", units), 1, 10)
    origin.pcict <- as.PCICt(origin, cal = caldata, format = "%Y-%m-%d")

    # distinguish between day and seconds based axis
    if (grepl("day", substr(units, 1, 3), fixed = TRUE)) {
      timeline <- origin.pcict + (floor(time) * 86400)
    } else if (grepl("sec", substr(units, 1, 3), fixed = TRUE)) {
      timeline <- origin.pcict + floor(time)
    } else if (grepl("hours", substr(units, 1, 5), fixed = TRUE)) {
      timeline <- origin.pcict + floor(time) * 3600
    } else {
      warning("Not recognised relative time axis!")
    }
  } else {
    warning("Not recognized time units!")
  }

  # warning for not-recognised calendar, stop if timeflag is on, replace time with input time
  if (any(is.na(timeline))) {
    print("Calendar from NetCDF is unsupported or not present")
    timeline <- time
  }

  return(timeline)
}


##########################################################
#--------------NetCDF loading function-------------------#
##########################################################

# universal function to open a single var 4D (x,y,z,time) ncdf files: it includes rotation, y-axis flipping, possible time selection and CDO-based interpolation
# automatically rotate matrix to place greenwich at the center (flag "rotate") and flip the latitudes in order to have increasing values
# if required (flag "interp2grid") additional interpolation with CDO can be used. "grid" can be used to specify the target grid name
# time selection based on package PCICt must be specifed with both "tmonths" and "tyears" flags
# level selection can be done with "tlev"
# it returns a list including its own dimensions
# last update in Jan-19
ncdf.opener.universal <- function(namefile, namevar = NULL, namelon = NULL, namelat = NULL, namelev = NULL, 
                                  tmonths = 1:12, tyears = NULL, tlev = NULL,
                                  rotate = "full", interp2grid = F, grid = "r144x73", remap_method = "remapbil",
                                  exportlonlat = FALSE, verbose = FALSE) {

  # load package
  require(ncdf4)

  # check if timeflag (by tyears) is activated or full file must be loaded
  if (is.null(tyears)) {
    timeflag <- FALSE
    printv("No years have been specified, loading all the data", verbose)
  } else {
    timeflag <- TRUE
    printv("tyears is set!", verbose)
  }

  # interpolation made with CDO: second order conservative remapping
  if (interp2grid) {
    print(paste("Remapping with CDO on", grid, "grid"))
    filename <- basename(normalizePath(namefile))
    filedir <- dirname(normalizePath(namefile))
    cdo <- Sys.which("cdo")
    tempfile <- paste0(file.path(filedir, paste0("tempfile_", filename)))
    system2("rm", tempfile)
    system2(cdo, args = c(paste0(remap_method, ",", grid), namefile, tempfile))
    namefile <- tempfile
  }

  # opening file: getting variable (if namevar is given, that variable is extracted)
  printv(paste("opening file:", namefile), verbose)
  a <- nc_open(namefile)
  print(paste("Loading", namevar, "..."))

  # if no name provided load the only variable available
  if (is.null(namevar)) {
    namevar <- names(a$var)
    if (length(namevar) > 1) {
      print(namevar)
      stop("More than one var in the files, please select it with namevar=yourvar. Stopping!!!")
    }
  } else {
    if (!(namevar %in% names(a$var))) {
      print(namevar)
      print(names(a$var))
      stop("Requested a non-existeng variable. Stopping!")
    }
  }

  # load axis: updated version, looking for dimension directly stored inside the variable
  naxis <- unlist(lapply(a$var[[namevar]]$dim, function(x) x["name"]))
  for (axis in naxis) {
    assign(axis, ncvar_get(a, axis))
    printv(paste(axis, ":", length(get(axis)), "records"), verbose)
  }

  # axis definition: this are standard axis names that are used for recognition
  xlist <- c("lon", "Lon", "longitude", "Longitude", "x")
  ylist <- c("lat", "Lat", "latitude", "Latitude", "y")
  zlist <- c("lev", "Lev", "plev", "Plev", "z")

  # function to assign netcdf file dimensions
  create.dimension <- function(namedim, listdim, axisdim) {
    if (is.null(namedim)) {
      if (any(listdim %in% axisdim)) {
        out <- get(axisdim[axisdim %in% listdim])
      } else {
        warning(paste("No",listdim[1]," found"))
        out <- NA
      } 
    } else {
      out <- get(namedim) 
    }
    return(out)
  }

  # create the three dimensions
  x <- create.dimension(namelon, xlist, naxis) 
  y <- create.dimension(namelat, ylist, naxis)
  z <- create.dimension(namelev, zlist, naxis)

  # if a time dimension exists, activate the PCICt package
  if (exists("time", mode = "numeric")) {
    require(PCICt)

    # extract units and calendary type to deal with PCICt package
    units <- ncatt_get(a, "time", "units")$value
    caldata <- ncatt_get(a, "time", "calendar")$value

    # create a PCICt objecti with create.timeline()
    timeline <- create.timeline(time, units, caldata, verbose = verbose)

    if (timeflag) {
      # day frequency case
      if (min(diff(timeline)) == 1) {
        lastday_base <- paste0(max(tyears), "-", max(tmonths), "-28") # uses number.days.month, which loops to get the month change
        lastday <- as.PCICt(paste0(max(tyears), "-", max(tmonths), "-", number.days.month(lastday_base, caldata)), 
                            cal = caldata, format = "%Y-%m-%d")
        firstday <- as.PCICt(paste0(min(tyears), "-", min(tmonths), "-01"), 
                             cal = caldata, format = "%Y-%m-%d")
      } 
      # monthly data case
      if (min(diff(timeline) >= 28)) {
        firstday <- as.PCICt(paste0(min(tyears), "-", min(tmonths), "-28"), cal = caldata, format = "%Y-%m-%d")
        lastday <- as.PCICt(paste0(min(tyears), "-", min(tmonths), "-01"), cal = caldata, format = "%Y-%m-%d")
      }
      # break if the data requested is not there
      if (max(trunc(timeline,"days")) < lastday | min(trunc(timeline,"days")) > firstday) {
        print(paste("Dataset extend from", min(timeline), "up to", max(timeline)))
        print(paste("Your request is from", firstday, "up to", lastday))
        stop("This a time interval that is not present in the NetCDF. Stopping!!!")
      }
    }
  } else {
    if (timeflag) {
      warning("Unknown time axis, disabling time selection!")
      timeflag <- FALSE
    }
  }


  # time selection and variable loading
  printv("loading full field...", verbose)
  field <- ncvar_get(a, namevar)

  if (timeflag) {

    # needed time window
    printv("selecting years and months", verbose)
    time_select <- which(as.numeric(format(timeline, "%Y")) %in% tyears & as.numeric(format(timeline, "%m")) %in% tmonths)

    # Selection of data: array_indexing() replace the previous 3d selection (Jan 2019)
    time_dim <- which(dim(field) == length(time))
    field <- array_indexing(field, time_dim, time_select)
    time <- timeline[time_select]

    printv(paste("This is a", caldata, "calendar"), verbose)
    printv(paste(length(time), "records selected from", time[1], "to", time[length(time)]), verbose)

    printv(paste("Months that have been loaded are... "),verbose)
    printv(unique(format(time, "%Y-%m")),verbose)
  }

  # Level selection  (Jan 2019)
  # if a selection on levels (tlev is defined and present) select
  if (!is.null(tlev)) {
    lev_select <- match(tlev, z) 
    if (any(is.na(lev_select))) {
       warning("Requested level not found, disabling level selection!")
    } else {
      lev_dim <- which(dim(field) == length(z))
      field <- array_indexing(field, lev_dim, lev_select, drop = T)
      z <- z[lev_select]
    }
    printv(paste("Level that have been loaded are... "),verbose)
    printv(z,verbose)
  }

  # check for dimensions (presence or not of time dimension)
  dimensions <- length(dim(field))
  if (dimensions > 4) {
    stop("This file is more than 4D file. Stopping!!!")
  }


  # if dimension have been recognized modified them
  if (any(!is.na(x))) {
    printv("rotating (if needed)...",verbose)
    x <- rotation(x, rotate)
    field <- rotation(field, rotate)
    assign(naxis[naxis %in% c(xlist, namelon)], x)
  }

  if (any(!is.na(y)) & length(y) > 1) {
    if (y[2] < y[1]) {
      print("flipping...")
      y <- sort(y)
      field <- flipper(field, which(dim(field) == length(y)))
    }
    assign(naxis[naxis %in% c(ylist, namelat)], y)
  }

  if (!is.na(z[1])) {
    assign(naxis[naxis %in% c(zlist, namelon)], z)
  }

  # exporting variables to the main program (weird but may be useful!)
  if (exportlonlat) {
    assign("ics", x, envir = .GlobalEnv)
    assign("ipsilon", y, envir = .GlobalEnv)
    assign("level", z, envir = .GlobalEnv)
  }

  # close connection
  nc_close(a)

  # remove interpolated file
  if (interp2grid) {
    system2("rm", tempfile)
  }

  # showing array properties
  printv(paste(dim(field)),verbose)
  printv(paste("From", time[1], "to", time[length(time)]),verbose)

  # returning file list
  return(mget(c("field", naxis)))
}

# ncdf.opener is a simplified wrapper for ncdf.opener.universal which returns only the field, ignoring the list and no verbosity
ncdf.opener <- function(namefile, namevar = NULL, namelon = NULL, namelat = NULL, namelev = NULL, 
                        tmonths = 1:12, tyears = NULL, tlev = NULL,
                        rotate = "full", interp2grid = F, grid = "r144x73", remap_method = "remapcon2",
                        exportlonlat = TRUE, verbose = FALSE) {
  field <- ncdf.opener.universal(namefile, namevar, namelon, namelat, namelev, tmonths, tyears, tlev,
                                 rotate, interp2grid, grid, remap_method, exportlonlat, verbose)
  return(field$field)
}

##########################################################
#--------------NetCDF writing functions------------------#
##########################################################

# simple function which defines the x,y,t,z dimensions to be written in NetCDF files
# produces also an "average" climatological time and monthly time axis for EOFs
# ics and ipsilon are vectors, tempus is a PCICt object, level a numeric
ncdf.defdims <- function(ics, ipsilon, tempus, level = 50000, reftime = "1850-01-01") {

  # use reftime to create daily time axis
  tcal <- attributes(tempus)$cal
  fulltime <- as.numeric(tempus - as.PCICt(reftime, cal = tcal)) + 86400 / 2
  nametime <- paste0("secs since ", reftime, " 00:00:00")

  # define core dimensions
  x <- ncdim_def("lon", "degrees_east", ics, longname = "longitude")
  y <- ncdim_def("lat", "degrees_north", ipsilon, longname = "latitude")
  t <- ncdim_def("time", nametime, fulltime, calendar = tcal, longname = "time", unlim = T)
  z <- ncdim_def("plev", "Pa", level, longname = "pressure")

  # climatological time (median time!)
  climtime <- as.numeric(median(tempus) - as.PCICt(reftime, cal = tcal))
  tclim <- ncdim_def("time", nametime, climtime, unlim = T, calendar = tcal, longname = "time")

  # monthly time (for EOFs)
  montime <- as.numeric(tempus[which(format(tempus, "%d") == 15)] - as.PCICt(reftime, cal = tcal))
  tmon <- ncdim_def("time", nametime, montime, calendar = tcal, longname = "time", unlim = T)

  # return the 3 basic dimensions
  return(list(x = x, y = y, t = t, tmon = tmon, tclim = tclim, z = z))
}

# zero-order common NetCDF common writer
# filename is the output file
# varlist and fieldlist are two list with respectively the netcdf variables definition and
# the field that should be associated to each variable.
# Their name must be the the name of the variables
# Of course varlist and fieldlist must have the same length.
ncdf.writer <- function(filename, varlist, fieldlist) {
  print(filename)

  # get the netcdr varlist and create the file
  ncfile <- nc_create(filename, varlist)

  savelist <- names(varlist)
  # for each var, check the dims and
  for (var in savelist) {
    print(var)
    ndims <- varlist[which(savelist == var)][[1]]$ndims
    ncvar_put(ncfile, var, fieldlist[which(savelist == var)][[1]], start = rep(1, ndims), count = rep(-1, ndims))
  }
  nc_close(ncfile)
}




##########################################################
#--------------Plotting functions------------------------#
##########################################################

# function to open devices
open.plot.device <- function(figname, output_file_type, CFGSCRIPT, special = FALSE) {
  # Chose output format for figure - by JvH
  source(CFGSCRIPT)
  if (special == FALSE) {
    if (tolower(output_file_type) == "png") {
      png(filename = figname, width = png_width, height = png_height)
    } else if (tolower(output_file_type) == "pdf") {
      pdf(file = figname, width = pdf_width, height = pdf_height, onefile = T)
    } else if (tolower(output_file_type) == "eps") {
      setEPS(width = pdf_width, height = pdf_height, onefile = T, paper = "special")
      postscript(figname)
    }
  }

  # special case for TM90
  if (special == TRUE) {
    if (tolower(output_file_type) == "png") {
      png(filename = figname, width = png_width / af, height = png_height * af / 2)
    } else if (tolower(output_file_type) == "pdf") {
      pdf(file = figname, width = pdf_width / af, height = pdf_height * af / 2, onefile = T)
    } else if (tolower(output_file_type) == "eps") {
      setEPS(width = pdf_width / af, height = pdf_height * af / 2, onefile = T, paper = "special")
      postscript(figname)
    }
  }
}


# extensive filled.contour function
filled.contour3 <-
  function(x = seq(0, 1, length.out = nrow(z)),
             y = seq(0, 1, length.out = ncol(z)), z, xlim = range(x, finite = TRUE),
             ylim = range(y, finite = TRUE), zlim = range(z, finite = TRUE),
             levels = pretty(zlim, nlevels), nlevels = 20, color.palette = cm.colors,
             col = color.palette(length(levels) - 1), extend = FALSE, plot.title, plot.axes,
             key.title, key.axes, asp = NA, xaxs = "i", yaxs = "i", las = 1,
             axes = TRUE, frame.plot = axes, mar, ...) {
    # modification by Ian Taylor of the filled.contour function
    # to remove the key and facilitate overplotting with contour()
    # further modified by Carey McGilliard and Bridget Ferris
    # to allow multiple plots on one page
    # modification to allow plot outside boundaries

    if (missing(z)) {
      if (!missing(x)) {
        if (is.list(x)) {
          z <- x$z
          y <- x$y
          x <- x$x
        }
        else {
          z <- x
          x <- seq.int(0, 1, length.out = nrow(z))
        }
      }
      else {
        stop("no 'z' matrix specified")
      }
    }
    else if (is.list(x)) {
      y <- x$y
      x <- x$x
    }
    if (any(diff(x) <= 0) || any(diff(y) <= 0)) {
      stop("increasing 'x' and 'y' values expected")
    }

    if (extend) {
      z[z < min(levels)] <- min(levels)
      z[z > max(levels)] <- max(levels)
    }

    plot.new()
    plot.window(xlim, ylim, "", xaxs = xaxs, yaxs = yaxs, asp = asp)
    if (!is.matrix(z) || nrow(z) <= 1 || ncol(z) <= 1) {
      stop("no proper 'z' matrix specified")
    }
    if (!is.double(z)) {
      storage.mode(z) <- "double"
    }
    .filled.contour(as.double(x), as.double(y), z, as.double(levels),
      col = col
    )
    if (missing(plot.axes)) {
      if (axes) {
        title(main = "", xlab = "", ylab = "")
        Axis(x, side = 1, ...)
        Axis(y, side = 2, ...)
      }
    }
    else {
      plot.axes
    }
    if (frame.plot) {
      box()
    }
    if (missing(plot.title)) {
      title(...)
    } else {
      plot.title
    }
    invisible()
  }

image.scale3 <- function(z, levels, color.palette = heat.colors, colorbar.label = "image.scale", extend = T,
                         line.label = 2, line.colorbar = 0, cex.label = 1, cex.colorbar = 1, colorbar.width = 1, ...) {

  # save properties from main plotting region
  old.par <- par(no.readonly = TRUE)
  mfg.save <- par()$mfg
  old.fig <- par()$fig

  # defining plotting region with proper scaling
  # print(old.fig)
  xscal <- (old.fig[2] - old.fig[1])
  yscal <- (old.fig[4] - old.fig[3])
  lw <- colorbar.width
  lp <- line.colorbar / 100
  new.fig <- c(old.fig[2] - 0.07 * xscal * lw - lp, old.fig[2] - 0.03 * xscal - lp, old.fig[3] + 0.1 * yscal, old.fig[4] - 0.1 * yscal)

  # safety check 
  new.fig[new.fig>1] <- 1
  #print(old.fig)
  #print(new.fig)

  if (missing(levels)) {
    levels <- seq(min(z), max(z), , 12)
  }
  # fixing color palette
  col <- color.palette(length(levels) - 1)

  # starting plot
  par(mar = c(1, 1, 1, 1), fig = new.fig, new = TRUE)

  # creating polygons for legend
  poly <- vector(mode = "list", length(col))
  for (i in seq(poly)) {
    poly[[i]] <- c(levels[i], levels[i + 1], levels[i + 1], levels[i])
  }

  xlim <- c(0, 1)
  if (extend) {
    longer <- 1.5
    dl <- diff(levels)[1] * longer
    ylim <- c(min(levels) - dl, max(levels) + dl)
  } else {
    ylim <- range(levels)
  }
  plot(1, 1, t = "n", ylim = ylim, xlim = xlim, axes = FALSE, xlab = "", ylab = "", xaxs = "i", yaxs = "i", ...)
  for (i in seq(poly)) {
    polygon(c(0, 0, 1, 1), poly[[i]], col = col[i], border = NA)
  }

  if (extend) {
    polygon(c(0, 1, 1 / 2), c(levels[1], levels[1], levels[1] - dl),
      col = col[1], border = NA
    )
    polygon(c(0, 1, 1 / 2), c(levels[length(levels)], levels[length(levels)], levels[length(levels)] + dl),
      col = col[length(col)], border = NA
    )
    polygon(c(0, 0, 1 / 2, 1, 1, 1 / 2), c(
      levels[1], levels[length(levels)], levels[length(levels)] + dl, levels[length(levels)], levels[1],
      levels[1] - dl
    ), border = "black", lwd = 2)
    ylim0 <- range(levels)
    prettyspecial <- pretty(ylim0)
    prettyspecial <- prettyspecial[prettyspecial <= max(ylim0) & prettyspecial >= min(ylim0)]
    axis(4, las = 1, cex.axis = cex.colorbar, at = prettyspecial, labels = prettyspecial, ...)
  } else {
    box()
    axis(4, las = 1, cex.axis = cex.colorbar, ...)
  }

  # box, axis and leged
  mtext(colorbar.label, line = line.label, side = 4, cex = cex.label, ...)

  # resetting properties for starting a new plot (mfrow style)
  par(old.par)
  par(mfg = mfg.save, new = FALSE)
  invisible()
}

# function for interpolation and projection of a 2D field on a mapproj R projection
proj.plot <- function(lon, lat, field, lmin = NULL, proj = "azequalarea", param = NULL, orient = c(90, 0, 0), npoints = 201) {

  # default is azimuthal equal area map

  # required packages
  require(mapproj)
  require(akima)

  # it provides lower latitude limit for plots
  if (is.null(lmin)) {
    lmin <- min(lat)
  }

  # build grids
  lon.grid <- rep(lon, length(lat))
  lat.grid <- sort(rep(lat, length(lon)))

  # project grid
  proj.grid <- mapproject(lon.grid, lat.grid, projection = proj, parameters = param, orientation = orient)

  # provide limits for future plots (for polar projection)
  limiter <- mapproject(c(0, 90, 180, 270), rep(lmin, 4), proj = "", orientation = orient)
  xlims <- sort(c(limiter$x[2], limiter$x[4]))
  ylims <- sort(c(limiter$y[1], limiter$y[3]))

  # plot grid
  lon.plot <- seq(min(proj.grid$x, na.rm = T), max(proj.grid$x, na.rm = T), length.out = npoints)
  lat.plot <- seq(min(proj.grid$y, na.rm = T), max(proj.grid$y, na.rm = T), length.out = npoints)

  # interpolation (akima needed)
  good <- is.finite(field) & is.finite(proj.grid$x) & is.finite(proj.grid$y)
  projected <- interp(proj.grid$x[good], proj.grid$y[good], field[good], lon.plot, lat.plot, duplicate = "strip")
  return(projected = list(x = projected$x, y = projected$y, z = projected$z, xlim = xlims, ylim = ylims))
}

# addland function based on map which can handle projections
proj.addland <- function(lon, lat, proj = "no", orient = c(90, 0, 0), param = NULL, inter = F, color = "black") {

  # required packages
  require(maps)
  require(mapproj)

  if (proj == "no") {
    map("world", regions = ".", interior = inter, exact = F, boundary = T, add = T)
  } else {
    # get map, project and do the lines
    box()
    map("world", add = T, projection = proj, orientation = orient, parameter = param, interior = inter, exact = F, boundary = T)

    # default lines for northern hemisphere
    for (i in seq(-80, 80, 20)) {
      x0 <- lon
      y0 <- rep(i, length(lon))
      p <- mapproject(x0, y0, proj = "", orientation = orient)
      lines(p, lty = 3)
    }

    # default circles for northern hemisphere
    for (i in c(seq(-360, 360, 30))) {
      y0 <- seq(0, 90, , 90)
      x0 <- rep(i, 90)
      p <- mapproject(x0, y0, proj = "", orientation = orient)
      lines(p, lty = 3)
    }
  }
}

# rearrange arrays for use both standard plotting and proj.plot
plot.prepare <- function(ics, ipsilon, field, proj, lat_lim) {
  if (proj == "no") {
    outfile <- list(x = ics, y = ipsilon, z = field, xlim = range(ics), ylim = lat_lim, xlab = "Longitude", ylab = "Latitude", axes = T)
  } else {
    field[is.na(field)] <- 0
    p <- proj.plot(ics, ipsilon, field, lmin = lat_lim[1], proj = proj, param = NULL, orient = c(90, 0, 0), npoints = 80)
    outfile <- list(x = p$x, y = p$y, z = p$z, xlim = p$xlim, ylim = p$ylim, xlab = "", ylab = "", axes = F)
  }
  return(outfile)
}

# function that provides labels and names for Blocking Plots
field.details <- function(field) {

  # default value
  legend_distance <- 3
  lev_hist <- NULL

  # case specific
  if (field == "TM90") {
    color_field <- c("dodgerblue", "darkred")
    color_diff <- NULL
    lev_field <- c(0, 30)
    lev_diff <- NULL
    legend_unit <- "Blocked Days (%)"
    title_name <- "Instantaneous Blocking (Tibaldi & Molteni, 1990):"
  }

  if (field == "InstBlock") {
    color_field <- palette1
    color_diff <- palette2
    lev_field <- seq(0, 27, 3)
    lev_diff <- seq(-10.5, 10.5, 1)
    legend_unit <- "Blocked Days (%)"
    title_name <- "Instantaneous Blocking frequency:"
  }

  if (field == "AbsBlock") {
    color_field <- palette1
    color_diff <- palette2
    lev_field <- seq(0, 27, 3)
    lev_diff <- seq(-10.5, 10.5, 1)
    legend_unit <- "Blocked Days (%)"
    title_name <- "Instantaneous Blocking frequency (Schwierz et al., 2004):"
  }

  if (field == "ExtraBlock") {
    color_field <- palette1
    color_diff <- palette2
    lev_field <- seq(0, 27, 3)
    lev_diff <- seq(-10.5, 10.5, 1)
    legend_unit <- "Blocked Days (%)"
    title_name <- "Instantaneous Blocking frequency (GHGS2 condition):"
  }

  if (field == "BlockEvents") {
    color_field <- palette1
    color_diff <- palette2
    lev_field <- seq(0, 27, 3)
    lev_diff <- seq(-10.5, 10.5, 1)
    lev_hist <- c(0, 16)
    legend_unit <- "Blocked Days (%)"
    title_name <- "Blocking Events frequency:"
  }

  if (field == "LongBlockEvents") {
    color_field <- palette1
    color_diff <- palette2
    lev_field <- seq(0, 16, 2)
    lev_diff <- seq(-5.25, 5.25, .5)
    legend_unit <- "Blocked Days (%)"
    title_name <- "10-day Blocking Events frequency:"
  }

  if (field == "DurationEvents") {
    color_field <- palette0
    color_diff <- palette2
    lev_field <- seq(5, 11.5, .5)
    lev_diff <- seq(-2.1, 2.1, .2)
    lev_hist <- c(6, 8)
    legend_unit <- "Duration (days)"
    title_name <- "Duration of Blocking Events:"
  }

  if (field == "NumberEvents") {
    color_field <- palette0
    color_diff <- palette2
    lev_field <- seq(0, 100, 10)
    lev_diff <- seq(-42.5, 42.5, 5)
    lev_hist <- c(0, 60)
    legend_unit <- ""
    title_name <- "Number of Blocking Events:"
  }

  if (field == "Z500") {
    color_field <- palette0
    color_diff <- palette2
    lev_field <- seq(4800, 6000, 50)
    lev_diff <- seq(-310, 310, 20)
    legend_unit <- "Geopotential Height (m)"
    title_name <- "Z500:"
    legend_distance <- 4
  }

  if (field == "BI") {
    color_field <- palette0
    color_diff <- palette2
    lev_field <- seq(1, 6, 0.25)
    lev_diff <- seq(-2.1, 2.1, .2)
    legend_unit <- "BI index"
    title_name <- "Blocking Intensity (BI):"
  }

  if (field == "MGI") {
    color_field <- palette0
    color_diff <- palette2
    lev_field <- seq(0, 15, 1)
    lev_diff <- seq(-5.25, 5.25, .5)
    legend_unit <- "MGI Index"
    title_name <- "Meridional Gradient Inversion (MGI):"
  }

  if (field == "ACN" | field == "CN") {
    if (field == "ACN") {
      title_name <- "Anticyclonic Rossby wave breaking frequency:"
    }
    if (field == "CN") {
      title_name <- "Cyclonic Rossby wave breaking frequency:"
    }
    color_field <- palette1
    color_diff <- palette2
    lev_field <- seq(0, 20, 2)
    lev_diff <- seq(-5.25, 5.25, .5)
    legend_unit <- "RWB frequency (%)"
  }

  out <- list(color_field = color_field, color_diff = color_diff, lev_field = lev_field,
	      lev_diff = lev_diff, lev_hist = lev_hist, legend_unit = legend_unit,
	      legend_distance = legend_distance, title_name = title_name)
  return(out)
}

##########################################################
#------------Blocking Tracking Functions-----------------#
##########################################################

# time persistence (used for longitude filter too)
time.persistence <- function(timeseries, persistence = 5) {
  rr <- rle(timeseries)
  rr$values[which(rr$values == 1 & rr$length < persistence)] <- 0
  nn <- rep(rr$values, rr$length)
  return(nn)
}


# blocking 5 days tracking
blocking.persistence <- function(field, minduration = 5, time.array) {

  # function for persistence
  pers2 <- function(timeseries, persistence, time.array) {
    dd <- min(time.array$season):max(time.array$season)
    nn <- sapply(dd, function(x) {
      time.persistence(timeseries[which(time.array$season == x)], persistence)
    })
    xx <- c(unlist(nn))
    return(xx)
  }

  # check for etime
  if (length(time.array$month) != length(field[1, 1, ])) {
    stop("Wrong time array! Exiting...")
  }

  print("Time filtering...")
  newfield <- apply(field, c(1, 2), function(x) pers2(x, persistence = minduration, time.array))
  newfield <- aperm(newfield, c(2, 3, 1))
  print("Mean field...")
  meanfield <- apply(newfield, c(1, 2), mean, na.rm = T) * 100


  print("Events detection...")
  maxdim <- max(apply(newfield, c(1, 2), function(x) length(rle(x)$length[which(rle(x)$values == 1)])))
  events <- apply(newfield, c(1, 2), function(x) c(rle(x)$lengths[which(rle(x)$values == 1)], rep(NA, maxdim - length(rle(x)$length[which(rle(x)$values == 1)]))))
  events <- aperm(events, c(2, 3, 1))
  print("Mean Duration...")
  duration <- apply(events, c(1, 2), mean, na.rm = T)
  print("Number of Events...")
  nevents <- apply(events, c(1, 2), function(x) length(x[!is.na(x)]))

  out <- list(track = newfield, percentage = meanfield, duration = duration, events = events, nevents = nevents)
  print(quantile(meanfield))
  print(min(duration, na.rm = T))
  return(out)
}


# large scale extension with further implementation
largescale.extension.if <- function(ics, ipsilon, field) {
  print("Large Scale Extension based on fixed angle")
  fimin <- 30 # southern latitude to be analyzed
  fimax <- 75 # northern latitude to be analyzed
  yreso <- ipsilon[2] - ipsilon[1]
  xreso <- ics[2] - ics[1]
  passo <- 5 / xreso # horizontal movemenent
  vertical <- 2.5 / yreso # vertical movement
  # time=1:length(field[1,1,]) #elements of the length of the dataset
  time <- which(apply(field, 3, max) != 0) # elements length of the dataset (removing no blocked days)

  print(paste("Box dimension:", passo * 2 * xreso, "° lon x ", vertical * 2 * yreso, "° lat"))

  short <- function(ics, ipsilon, field, passo, vertical) {
    control <- field
    range <- which.min(abs(ipsilon - fimin)):which.min(abs(ipsilon - fimax)) # check range for latitude excursion
    # range=range[(1+vertical):(length(range)-vertical)] #reduce range considering border effect
    new <- rbind(field, field, field) # bind domain for cross-date line
    for (i in 1:length(ics)) {
      ii <- i + length(ics)
      if (!all(new[(ii - passo):(ii + passo), ] == 0)) { # check to speed up
        for (j in range) {
          control[i, j] <- mean(new[(ii - passo):(ii + passo), (j - vertical):(j + vertical)], na.rm = T)
        }
      }
    }
    control[control > 0] <- 1
    return(control)
  }

  tt <- length(time)
  for (t in time) {
    progression.bar(t, tt)
    field[, , t] <- short(ics, ipsilon, field[, , t], passo, vertical)
  }
  return(field)
}


# Longitude filter for minimum extension
longitude.filter <- function(ics, ipsilon, field) {
  print("Longitude filter based on fixed angle")
  out <- field
  yreso <- ipsilon[2] - ipsilon[1]
  xreso <- ics[2] - ics[1]
  startipsilon <- which.min(abs(ipsilon - 30))
  estension <- (75 - 30) / yreso
  passo <- 15 / xreso

  print(paste("Continous longitude contrain", passo * xreso, "° lon"))

  tt <- length(field[1, 1, ])
  for (t in 1:tt) {
    progression.bar(t, tt)

    new <- rbind(field[, , t], field[, , t], field[, , t])
    for (j in startipsilon:((startipsilon + estension))) {
      new[, j] <- time.persistence(new[, j], persistence = passo)
    }
    field[, , t] <- new[length(ics) + (1:length(ics)), ]
  }
  return(field)
}


##########################################################
#------------EOFs and regims functions-------------------#
##########################################################

eofs <- function(lon, lat, field, neof = 4, xlim = NULL, ylim = NULL, 
                 method = "SVD", do_standardize = F, do_regression = F, verbose = T) {
  # R tool for computing EOFs based on Singular Value Decomposition ("SVD", default)
  # or with the eigenvectors of the covariance matrix ("covariance", slower)
  # If requested, computes linear regressions and standardizes the PCs
  # If you want to use the regressions, remember to standardize the PCs
  # Take as input a 3D anomaly field.
  # Requires "personal" functions area.weight, whicher and standardize

  # area weighting, based on the root of cosine
  printv("Area Weighting...", verbose)
  ww <- area.weight(lon, lat, root = T)
  wwfield <- sweep(field, c(1, 2), ww, "*")

  # selection of the xbox and ybox if defined
  if (!is.null(xlim)) {
    lonselect <- whicher(lon, xlim[1]):whicher(lon, xlim[2])
  } else {
    lonselect <- 1:length(lon) 
  }

  if (!is.null(ylim)) {
    latselect <- whicher(lat, ylim[1]):whicher(lat, ylim[2])
  } else {
    latselect <- 1:length(lat)
  }

  # box
  box <- wwfield[lonselect, latselect, ]
  slon <- lon[lonselect]
  slat <- lat[latselect]

  # transform 3D field in a matrix
  new_box <- array(box, dim = c(dim(box)[1] * dim(box)[2], dim(box)[3]))

  # calling SVD
  if (method == "SVD") {
    printv("Calling SVD...", verbose)
    SVD <- svd(new_box, nu = neof, nv = neof)

    # extracting EOFs (loading pattern), expansions coefficient and variance explained
    pattern <- array(SVD$u, dim = c(dim(box)[1], dim(box)[2], neof))
    coefficient <- SVD$v
    variance <- (SVD$d[1:neof])^2 / sum((SVD$d)^2)
    if (do_standardize) {
      coefficient <- apply(coefficient, c(2), standardize)
    } else {
      coefficient <- sweep(coefficient, c(2), sqrt(variance), "*")
    }
  }

  # calling covariance matrix
  if (method == "covariance") {
    printv("Calling eigenvectors of the covariance matrix...", verbose)
    covma <- cov(t(new_box))
    eig <- eigen(covma)
    coef <- (t(new_box) %*% eig$vector)[, 1:neof]
    pattern <- array(eig$vectors, dim = c(dim(box)[1], dim(box)[2], dim(box)[3]))[, , 1:neof]
    variance <- eig$values[1:neof] / sum(eig$values)
    if (do_standardize) {
      coefficient <- apply(coef, c(2), standardize)
    } else {
      coefficient <- coef
    }
  }

  # linear regressions on anomalies
  regression <- NULL
  if (do_regression) {
    printv("Linear Regressions (it can takes a while)... ", verbose)
    regression <- array(NA, dim = c(length(lon), length(lat), neof))
    # for (i in 1:neof) {regression[,,i]=apply(field,c(1,2),function(x) coef(lm(x ~ coefficient[,i]))[2])}
    for (i in 1:neof) {
      regression[, , i] <- apply(field, c(1, 2), function(x) lin.fit(as.matrix(coefficient[, i], ncol = 1), x)$coefficients)
    }
  }

  # preparing output
  printv("Finalize...", verbose)
  pattern <- list(x = slon, y = slat, z = pattern)
  out <- list(pattern = pattern, coeff = coefficient, variance = variance, regression = regression)
  return(out)
}

eofs.coeff <- function(lon, lat, field, eof_object, do_standardize = F, verbose = F) {
  # Computes expansion coefficient (i.e. PCs) of a given dataset on the
  # loading pattern of EOF previously computed
  # Works only on eof_object obtained with "eofs" function

  # Area weighting, based on the root of cosine
  printv("Area Weighting...", verbose)
  ww <- area.weight(lon, lat, root = T)
  wwfield <- sweep(field, c(1, 2), ww, "*")

  # selection of the box
  xlim <- c(min(eof_object$pattern$x), max(eof_object$pattern$x))
  ylim <- c(min(eof_object$pattern$y), max(eof_object$pattern$y))
  box <- wwfield[whicher(lon, xlim[1]):whicher(lon, xlim[2]), whicher(lat, ylim[1]):whicher(lat, ylim[2]), ]

  # transform 3D field in a matrix
  new_box <- array(box, dim = c(dim(box)[1] * dim(box)[2], dim(box)[3]))
  new_pattern <- array(eof_object$pattern$z, dim = c(dim(eof_object$pattern$z)[1] * dim(eof_object$pattern$z)[2], dim(eof_object$pattern$z)[3]))

  # projects the coefficients
  coef <- (t(new_box) %*% new_pattern)

  # standardize
  if (do_standardize) {
    coefficient <- apply(coef, c(2), standardize)
  } else {
    coefficient <- coef
  }

  print("Finalize...")
  return(coefficient)
}


regimes <- function(lon, lat, field, ncluster = 4, ntime = 1000, neof = 10, xlim, ylim, alg = "Hartigan-Wong") {
  # R tool to compute cluster analysis based on k-means.
  # Requires "personal" function eofs
  # Take as input a 3D anomaly field

  # Reduce the phase space with EOFs: use SVD and do not standardize PCs
  print("Launching EOFs...")
  t0 <- proc.time()
  reducedspace <- eofs(lon, lat, field, neof = neof, xlim = xlim, ylim = ylim, method = "SVD", do_regression = F, do_standardize = F)
  t1 <- proc.time() - t0
  # print(t1)

  # extract the principal components
  PC <- reducedspace$coeff
  print(str(PC))

  # k-means computation repeat for ntime to find best solution.
  print("Computing k-means...")
  t0 <- proc.time()
  print(str(ncluster))
  regimes <- kmeans(PC, as.numeric(ncluster), nstart = ntime, iter.max = 1000, algorithm = alg)
  t1 <- proc.time() - t0
  # print(t1)

  # Extract regimes frequencyr and timeseries of occupation
  cluster <- regimes$cluster
  frequencies <- regimes$size / dim(field)[3] * 100
  print(frequencies[order(frequencies, decreasing = T)])
  # print(regimes$tot.withinss)

  print("Creating Composites...")
  compose <- aperm(apply(field, c(1, 2), by, cluster, mean), c(2, 3, 1))

  # sorting from the more frequent to the less frequent
  kk <- order(frequencies, decreasing = T)
  cluster <- cluster + 10
  for (ss in 1:ncluster) {
    cluster[cluster == (ss + 10)] <- which(kk == ss)
  }

  # prepare output
  print("Finalize...")
  out <- list(cluster = cluster, frequencies = frequencies[kk], regimes = compose[, , kk], tot.withinss = regimes$tot.withinss)
  return(out)
}


regimes2 <- function(lon, lat, field, ncluster = 4, ntime = 1000, minvar = 0.8,
                     xlim, ylim, alg = "Hartigan-Wong") {

  # R tool to compute cluster analysis based on k-means.
  # Requires "personal" function eofs (see above)
  # Take as input a 3D anomaly field

  # Reduce the phase space with EOFs: use SVD and do not standardize PCs
  print("Launching EOFs...")
  t0 <- proc.time()
  reducedspace <- eofs(lon, lat, field, neof = 20, xlim = xlim, ylim = ylim, method = "SVD", do_regression = F, do_standardize = F)
  t1 <- proc.time() - t0
  # print(t1)
  reqPC <- which(cumsum(reducedspace$variance) > minvar)[1]
  print(paste("Retaining", reqPC, "PCs to fullfil minimum explained variance required (", minvar * 100, "%)"))

  # extract the principal components
  PC <- reducedspace$coeff[, 1:reqPC]
  print(str(PC))

  # k-means computation repeat for ntime to find best solution.
  print("Computing k-means...")
  t0 <- proc.time()
  print(str(ncluster))
  regimes <- kmeans(PC, as.numeric(ncluster), nstart = ntime, iter.max = 100, algorithm = alg)
  t1 <- proc.time() - t0
  # print(t1)

  # Extract regimes frequencyr and timeseries of occupation
  cluster <- regimes$cluster
  frequencies <- regimes$size / dim(field)[3] * 100
  print(frequencies[order(frequencies, decreasing = T)])
  # print(regimes$tot.withinss)

  print("Creating Composites...")
  compose <- aperm(apply(field, c(1, 2), by, cluster, mean), c(2, 3, 1))

  # sorting from the more frequent to the less frequent
  kk <- order(frequencies, decreasing = T)
  cluster <- cluster + 10
  for (ss in 1:ncluster) {
    cluster[cluster == (ss + 10)] <- which(kk == ss)
  }

  # prepare output
  print("Finalize...")
  out <- list(cluster = cluster, frequencies = frequencies[kk], regimes = compose[, , kk], tot.withinss = regimes$tot.withinss)
  return(out)
}

##########################################################
#-------------------Time Avg functions-------------------#
##########################################################

# generalized function for time averaging based on conditon, using preallocation, vectorization and rowMeans
# use power.date.new or PCICt object to define the condition
time.mean <- function(ics, ipsilon, field, condition) {
  tmean <- array(NA, dim = c(length(ics), length(ipsilon), length(unique(condition))))
  for (t in unique(condition)) {
    tmean[, , which(t == unique(condition))] <- rowMeans(field[, , which(t == condition)], dims = 2)
  }
  return(tmean)
}


# fast function for monthly mean, using preallocation, vectorization and rowMeans
monthly.mean <- function(ics, ipsilon, field, etime) {
  condition <- paste(etime$month, etime$year)
  monthly <- array(NA, dim = c(length(ics), length(ipsilon), length(unique(condition))))
  for (t in unique(condition)) {
    monthly[, , which(t == unique(condition))] <- rowMeans(field[, , t == condition], dims = 2)
  }
  return(monthly)
}

# introduce running mean, options for loop dataset
run.mean <- function(field, n = 5, loop = F) {
  nn <- floor(n / 2)
  if (loop) {
    field = c(tail(field, n), field, head(field, n))
  }

  runfield <- field*NA
  for (t in (1 + nn):(length(field) - nn)) {
    if (!is.na(field[t])) {
      runfield[t] <- mean(field[(t - nn):(t + nn)], na.rm = T)
    } else {
      runfield[t] <- NA 
    }
  }

  if (loop) {
    runfield <- runfield[(1+n):(length(runfield)-n)] 
  }

  return(runfield)
}

# improve running mean
# use vectorization for a 5 day running mean ad-hoc function (to be generalized!)
# about 10 times faster that a standard running mean function based on for loop
run.mean5 <- function(field) {
  newfield <- rowMeans(cbind(
    c(field[3:length(field)], NA, NA),
    c(field[2:length(field)], NA), field,
    c(NA, field[1:(length(field) - 1)]),
    c(NA, NA, field[1:(length(field) - 2)])
  ),
  na.rm = T
  )
  return(newfield)
}

# this a generalization of run.mean5() using vectorization for a faster result
# result is more accurate than run.mean() since on the edges it uses all possible points
run.mean.fast <- function(field, n = 5) {

  # check for even numbers
  if (n %% 2 == 0) {
    warning("Even number, replacing with its smaller odd one")
    n <- n - 1
  }

  # prepare the loop and create n-dimensianal matrix of lagged vectors
  nn <- floor(n / 2)
  newfield <- NULL
  for (k in -nn:nn) {
    if (k < 0) {
      newfield <- cbind(newfield, c(rep(NA, abs(k)), field[1:(length(field) + k)]))
    } else {
      newfield <- cbind(newfield, c(field[(1 + k):length(field)], rep(NA, k)))
    }
  }

  # apply rowMeans to produce the final mean
  finalfield <- rowMeans(newfield, na.rm = T)
  return(finalfield)
}



# function for daily anomalies, use array predeclaration and rowMeans (40 times faster!)
daily.anom.mean <- function(ics, ipsilon, field, etime) {
  condition <- paste(etime$day, etime$month)
  daily <- array(NA, dim = c(length(ics), length(ipsilon), length(unique(condition))))
  anom <- field * NA
  for (t in unique(condition)) {
    if (sum(t == condition) == 1) {
      stop("Cannot compute a mean with a single value")
    }
    daily[, , which(t == unique(condition))] <- rowMeans(field[, , t == condition], dims = 2)
    anom[, , which(t == condition)] <- sweep(field[, , which(t == condition)], c(1, 2), daily[, , which(t == unique(condition))], "-")
  }
  return(anom)
}

# beta function for daily anomalies plus running mean (only 50% slower that standard daily avg)
daily.anom.run.mean5.old <- function(ics, ipsilon, field, etime) {
  condition <- paste(etime$day, etime$month)
  daily <- array(NA, dim = c(length(ics), length(ipsilon), length(unique(condition))))
  for (t in unique(condition)) {
    if (sum(t == condition) == 1) {
      stop("Cannot compute a mean with a single value")
    }
    daily[, , which(t == unique(condition))] <- rowMeans(field[, , t == condition], dims = 2)
  }
  rundaily <- apply(daily, c(1, 2), run.mean5)
  anom <- field * NA
  for (t in unique(condition)) {
    anom[, , which(t == condition)] <- sweep(field[, , which(t == condition)], c(1, 2), daily[, , which(t == unique(condition))], "-")
  }
  return(anom)
}

# function for daily anomalies from a 5-day running mean seasonal cycle:
# This that takes into account correctly the cross-year season as DJF
daily.anom.run.mean5 <- function(ics, ipsilon, field, etime) {

  # define condition for time selection: use numeric-compatible format
  condition <- format(etime$data,"%m%d")

  # evalute the time-ordered condition: if there is a jump, it means that there is a cross-year
  sorted <- sort(unique(condition))
  breakpoint <- which(diff(as.numeric(sorted)) > 100)
  #print(breakpoint)
  
  # if there is a cross-year, re-arrenge in order to have them as consecutive dates
  if (length(breakpoint) > 0 ) {
    sorted = c(sorted[(breakpoint + 1):length(sorted)], sorted[1:breakpoint])
  }

  #print(sorted)
  # compute the seasonal cycle
  daily <- array(NA, dim = c(length(ics), length(ipsilon), length(sorted)))
  for (t in sorted) {
    if (sum(t == condition) == 1) {
      warning("Cannot compute a mean with a single value: using the value as it is!")
      daily[, , which(t == sorted)] <- field[, , t == condition]
    } else {
      daily[, , which(t == sorted)] <- rowMeans(field[, , t == condition], dims = 2)
    }
  }

  # apply running mean on the rightly orderded seasonal cycle
  rundaily <- apply(daily, c(1, 2), run.mean5)
  
  # remove seasonal cycle
  anom <- field * NA
  for (t in sorted) {
    anom[, , which(t == condition)] <- sweep(field[, , which(t == condition)], c(1, 2), daily[, , which(t == sorted)], "-")
  }
  
  return(anom)
} 

