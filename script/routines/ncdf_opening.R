##########################################################
#--------------NetCDF loading function-------------------#
##########################################################

# universal function to open a single var 4D (x,y,z,time) ncdf files: it includes rotation, y-axis flipping, possible time selection and CDO-based interpolation
# automatically rotate matrix to place greenwich at the center (flag "rotate") and flip the latitudes in order to have increasing values
# if required (flag "interp2grid") additional interpolation with CDO can be used. "grid" can be used to specify the target grid name
# time selection based on package PCICt can be specifed with both "tmonths" and "tyears" flags
# level selection can be done with "tlev"
# exportlonlat: true to export lon lat to global environment
# force: ignore the data boundaries
# verbose: expands information
# it returns a list including its own dimensions
# last update in Oct-19
ncdf.opener.universal <- function(namefile, namevar = NULL, namelon = NULL, namelat = NULL, namelev = NULL,
                                  tmonths = NULL, tyears = NULL, tlev = NULL,
                                  rotate = "full", interp2grid = F, grid = "r144x73", remap_method = "remapbil",
                                  exportlonlat = FALSE, force = FALSE, verbose = FALSE) {

  # load package
  require(ncdf4)

  # check if timeflag (by tyears or tmomths) is activated or full file must be loaded
  if (is.null(tyears) & is.null(tmonths)) {
    timeflag <- FALSE
    printv("No years or months have been specified, loading all the data", verbose)
  } else {
    timeflag <- TRUE
    # load all the data if tmonths is not defined
    if (is.null(tmonths)) {
      tmonths <- 1:12
    }
    printv("tyears or tmonths are set!", verbose)
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
  printv(paste("Loading", namevar, "..."), verbose)
  # load axis: updated version, looking for dimension directly stored inside the variable
  naxis <- unlist(lapply(a$var[[namevar]]$dim, function(x) x["name"]))
  for (axis in naxis) {
    assign(axis, ncvar_get(a, axis))
    printv(paste(axis, ":", length(get(axis)), "records"), verbose)
  }

  # axis definition: these are standard axis names that are used for recognition
  xlist <- c("lon", "Lon", "longitude", "Longitude", "x")
  ylist <- c("lat", "Lat", "latitude", "Latitude", "y")
  zlist <- c("lev", "Lev", "plev", "Plev", "z")

  # function to assign netcdf file dimensions
  create.dimension <- function(namedim, listdim, axisdim) {
    if (is.null(namedim)) {
      if (any(listdim %in% axisdim)) {
        out <- get(axisdim[axisdim %in% listdim])
      } else {
        warning(paste("No", listdim[1], " found"))
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

  # define a time varible, time or time_counter
  if (exists("time", mode = "numeric")) {
    timevariable <- "time"
  } else if (exists("time_counter", mode = "numeric")) {
    timevariable <- "time_counter"
  }

  # if a time dimension exists, activate the PCICt package
  if (exists("timevariable")) {
    printv(print.break(), verbose)
    printv("Time axis found, trying to recognize it...", verbose)
    require(PCICt)

    # extract units and calendary type to deal with PCICt package
    units <- ncatt_get(a, timevariable, "units")$value
    caldata <- ncatt_get(a, timevariable, "calendar")$value

    # create a PCICt objecti with create.timeline()
    timeinfo <- create.timeline(get(timevariable), units, caldata, verbose = verbose)
    time <- timeline <- timeinfo$timeline

    # optional information on dataset properties: only informative
    if (length(timeline) == 1) {
      printv("Single-instant value found!", verbose)
    } else {
      delta_time <- min(diff(timeline))
      if (delta_time == 1) {
        printv("Daily data found!", verbose)
      } else if (delta_time >= 28 & delta_time <= 31) {
        printv("Monthly data found!", verbose)
      } else if (delta_time > 360) {
        printv("Yearly data found!", verbose)
      } else {
        printv(paste("Time difference of", delta_time, "days: irregular time axis?"), verbose)
      }
    }

    if (timeinfo$calendar == "unsupported") {
      print("Unsupported calendar, disabling time selection...")
      timeflag <- FALSE
    }
  }

  # time selection and variable loading
  printv(print.break(), verbose)
  printv("loading full field...", verbose)
  field <- ncvar_get(a, namevar)

  if (timeflag) {

    # needed time window: updated in May 2019 to be more versatile
    # removed check on months/daily time axis
    printv("selecting years and months", verbose)
    if (is.null(tyears)) {
      time_select <- which(as.numeric(format(timeline, "%m")) %in% tmonths)
      ff <- format(timeline, "%m")
      gg <- sprintf("%02d", tmonths)
    } else {
      time_select <- which(as.numeric(format(timeline, "%Y")) %in% tyears & as.numeric(format(timeline, "%m")) %in% tmonths)
      ff <- format(timeline, "%Y-%m")
      gg <- paste(rep(tyears, each = length(tmonths)), sprintf("%02d", tmonths), sep = "-")
    }
    if (length(time_select) == 0) {
      print(paste("Dataset extend from", min(timeline), "up to", max(timeline)))
      stop("You requested a time interval not present in the NetCDF")
    }
    if (!all(gg %in% ff)) {
      print("WARNING: Not all the data you requested has been loaded")
      print("WARNING: Following years/months are not present in the NetCDF:")
      print(gg[!gg %in% ff])
    }

    # Selection of data: array_indexing() replace the previous 3d selection (Jan 2019)
    time_dim <- which(dim(field) == length(time))
    field <- array_indexing(field, time_dim, time_select)
    time <- timeline[time_select]

    printv(paste(length(time), "records selected from", time[1], "to", time[length(time)]), verbose)

    printv(paste("Months that have been loaded are... "), verbose)
    printv(unique(format(time, "%Y-%m")), verbose)
    printv(print.break(), verbose)
  }

  # Level selection  (Jan 2019)
  # if a selection on levels (tlev is defined and present) select
  if (!is.null(tlev)) {
    lev_select <- match(tlev, z)
    if (any(is.na(lev_select))) {
      warning("Requested level not found, disabling level selection!")
    } else {
      lev_dim <- which(dim(field) == length(z))
      if (length(lev_dim) == 0) {
        printv(paste("Only one level to be loaded... "), verbose)
      } else {
        field <- array_indexing(field, lev_dim, lev_select, drop = T)
        z <- z[lev_select]
      }
    }
    printv(paste("Level that have been loaded are... "), verbose)
    printv(z, verbose)
  }

  # check for dimensions (presence or not of time dimension)
  dimensions <- length(dim(field))
  if (dimensions > 4) {
    stop("This file is more than 4D file. Stopping!!!")
  }


  # if dimension have been recognized modified them
  if (any(!is.na(x))) {
    printv("rotating (if needed)...", verbose)
    x <- rotation(x, rotate, longitude.case = TRUE)
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
    assign(naxis[naxis %in% c(zlist, namelev)], z)
  }


  # extract units
  var_units <- ncatt_get(a, namevar, "units")$value

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
  printv(paste("Final array dimension:", paste(dim(field), collapse = " ")), verbose)

  # returning file list
  return(mget(c("field", naxis, "var_units")))
}

# ncdf.opener is a simplified wrapper for ncdf.opener.universal which returns only the field, ignoring the list and no verbosity
ncdf.opener <- function(namefile, namevar = NULL, namelon = NULL, namelat = NULL, namelev = NULL,
                        tmonths = NULL, tyears = NULL, tlev = NULL,
                        rotate = "full", interp2grid = F, grid = "r144x73", remap_method = "remapcon2",
                        exportlonlat = TRUE, verbose = FALSE) {
  field <- ncdf.opener.universal(
    namefile, namevar, namelon, namelat, namelev, tmonths, tyears, tlev,
    rotate, interp2grid, grid, remap_method, exportlonlat, verbose
  )
  return(field$field)
}
