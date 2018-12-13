######################################################
#-----Blocking routines computation for MiLES--------#
#-------------P. Davini (Oct 2014)-------------------#
######################################################
miles.block.fast <- function(project, dataset, expid, ens, year1, year2, season, nomefile, FILESDIR, doforce) {

  # rm(list=ls())
  # PROGDIR="/home/paolo/MiLES"
  # FILESDIR="/work/users/paolo/miles/files"
  # source(file.path(PROGDIR,"script/basis_functions.R"))
  # dataset="ERAI"
  # year1=1979
  # year2=1990
  # season="DJF"
  # nomefile="/work/users/paolo/miles/data/ua500/ERAI/ua500_ERAI_fullfile.nc"
  # project=expid=ens=NA
  # doforce="true"

  # t0
  t0 <- proc.time()

  # setting up time domain
  years <- year1:year2
  timeseason <- season2timeseason(season)

  # define folders using file.builder function (takes care of ensembles)
  savefile1 <- file.builder(FILESDIR, "U500_Block", "BlockClim", project, dataset, expid, ens, year1, year2, season)
  savefile2 <- file.builder(FILESDIR, "U500_Block", "BlockFull", project, dataset, expid, ens, year1, year2, season)

  # check if data is already there to avoid re-run
  if (file.exists(savefile1) & file.exists(savefile2)) {
    print("Actually requested blocking data is already there!")
    print(savefile1)
    print(savefile2)
    if (as.logical(doforce)) {
      print("Running with doforce=true... re-run!")
    } else {
      print("Skipping... activate doforce=true if you want to re-run it")
      q()
    }
  }

  # new file opening
  fieldlist <- ncdf.opener.universal(nomefile, namevar = "ua", tmonths = timeseason, tyears = years, rotate = "full", exportlonlat = F)
  print(str(fieldlist))

  # assign variables
  ics <- fieldlist$lon
  ipsilon <- fieldlist$lat

  # extract calendar and time unit from the original file
  timeaxis <- fieldlist$time

  # time array to simplify time filtering
  etime <- power.date.new(timeaxis, verbose = T)
  totdays <- length(timeaxis)

  # declare variable
  U500 <- fieldlist$field
  rm(fieldlist)

  # grid resolution
  yreso <- ipsilon[2] - ipsilon[1]
  xreso <- ics[2] - ics[1]

  # reso checks: this are not needed with default 2.5 grid, but they may be relevant with
  # future envisaged power up to finer grids
  # xcritical factor is due to RWB longitudinal delta of 7.5
  # ycritical factor is due to Large Scale Extension of 2.5
  xcritical <- 2.5
  ycritical <- 2.5
  if (ycritical %% yreso != 0) {
    stop("Latitudinal resolution is not a factor of 5 deg")
  }

  if (xcritical %% xreso != 0) {
    stop("Longitudinal resolution is not a factor of 5 deg")
  }

  ##########################################################
  #--------------Tibaldi and Molteni 1990------------------#
  ##########################################################

  print("Zonald wind based tibaldi and Molteni (1990) index...")
  # TM90: parametres for blocking detection
  tm90_fi0 <- 60 # central_lat
  tm90_delta <- 20 # delta in degree from between lats
  tm90_fiN <- tm90_fi0 + tm90_delta
  tm90_fiS <- tm90_fi0 - tm90_delta # south and north lat, 80N and 40N
  tm90_central <- whicher(ipsilon, tm90_fi0)
  tm90_south <- whicher(ipsilon, tm90_fiS)
  tm90_north <- whicher(ipsilon, tm90_fiN)
  tm90_range <- seq(-5, 5, yreso) / yreso # 5 degrees to the north, 5 to the south (larger than TM90 or D'Andrea et al 1998)

  # from u500: define parameters for geostrophic approximation
  g0 <- 9.80655
  radius <- 6378137
  omega <- 7.29211 * 10^(-5)
  sinphi <- sin(ipsilon * pi / 180)
  alfa <- (2 * radius * omega / g0) * ((tm90_delta * pi / 180))
  ghgn <- ghgs <- array(NA, dim = c(length(ics), totdays, length(tm90_range)))

  # number of steps and weights for integrals
  tm90_step <- round(tm90_delta / yreso)
  tm90_ww <- c(0.5, rep(1, tm90_step - 1), 0.5)
  tm90_ww <- c(0.5, rep(1, tm90_step - 1), 0.5)


  # introduce matrix for blocking computation

  # vectorization but on different ranges
  for (erre in tm90_range) {
    ghgs[, , which(erre == tm90_range)] <- apply(sweep(
      U500[, erre + tm90_south:tm90_central, ], c(2),
      sinphi[erre + tm90_south:tm90_central] * tm90_ww, "*"
    ), c(1, 3), sum) / tm90_step
    ghgn[, , which(erre == tm90_range)] <- apply(sweep(
      U500[, erre + tm90_central:tm90_north, ], c(2),
      sinphi[erre + tm90_central:tm90_north] * tm90_ww, "*"
    ), c(1, 3), sum) / tm90_step
  }

  # adapted condition from geostropich approximatin
  tm90_check <- (ghgs < 0 & ghgn > (10 * tm90_delta / alfa)) # D'Andrea et al. 97 condition: adapted to Scaife et al. 2010
  tm90_check[tm90_check == T] <- 1
  tm90_check[tm90_check == F] <- 0
  # blocked=apply(tm90_check,c(1),max,na.rm=T)
  # frequency=apply(blocked,c(1),mean,na.rm=T)*100
  totTM90 <- apply(tm90_check, c(1, 2), max, na.rm = T)
  TM90 <- apply(totTM90, 1, mean) * 100
  print("Done!")


  ##########################################################
  #--------------Davini et al. 2012------------------------#
  ##########################################################


  # decleare main variables to be computed (considerable speed up!)
  totblocked <- totblocked2 <- U500 * 0

  # Davini et al. 2012: parameters to be set for blocking detection
  fi0 <- 30 # lowest latitude to be analyzed
  delta <- 15 # distance on which compute gradients
  step <- round(delta / yreso)
  central <- which.min(abs(ipsilon - fi0)) # lowest starting latitude
  north <- central + step # lowest north latitude
  south <- central - step # lowest sourth latitude
  maxsouth <- central - 2 * step
  fiN <- ipsilon[north]
  fiS <- ipsilon[south]
  range <- (90 - fi0 - delta) / yreso # escursion to the north for computing blocking (from 30 up to 75)

  # number of steps and weights for integrals
  ww <- c(0.5, rep(1, step - 1), 0.5)

  print("--------------------------------------------------")
  print("Zonal wind Davini et al. (2012) index and diagnostics...")
  print(c("distance for gradients:", step * yreso))
  print(paste("range of latitudes ", fi0, "-", 90 - step * yreso, " N", sep = ""))

  ##########################################################
  #--------------Istantaneous Blocking---------------------#
  ##########################################################

  #----COMPUTING BLOCKING INDICES-----
  for (t in 1:totdays) {
    progression.bar(t, totdays)

    ghgn <- ghgs <- gh2gs <- array(NA, dim = c(length(ics), length(0:range)))

    for (erre in 0:range) { # computing blocking for different latitudes

      ghgs[, which(erre == (0:range))] <- apply(sweep(
        U500[, south:central + erre, t], c(2),
        sinphi[south:central + erre] * ww, "*"
      ), c(1), sum) / step
      ghgn[, which(erre == (0:range))] <- apply(sweep(
        U500[, central:north + erre, t], c(2),
        sinphi[central:north + erre] * ww, "*"
      ), c(1), sum) / step
      gh2gs[, which(erre == (0:range))] <- apply(sweep(
        U500[, maxsouth:south + erre, t], c(2),
        sinphi[maxsouth:south + erre] * ww, "*"
      ), c(1), sum) / step
    }
    check1 <- (ghgs < 0 & ghgn > (10 * delta / alfa))
    check1[check1 == T] <- 1
    check1[check1 == F] <- 0
    totblocked[, central:(central + range), t] <- check1

    check2 <- (ghgs < 0 & ghgn > (10 * delta / alfa) & gh2gs > (5 * delta / alfa))
    check2[check2 == T] <- 1
    check2[check2 == F] <- 0
    totblocked2[, central:(central + range), t] <- check2
  }
  print(paste("Total # of days:", t))
  print("-------------------------")

  ##########################################################
  #--------------------Mean Values-------------------------#
  ##########################################################

  # compute mean values (use rowMeans that is faster when there are no NA values)
  frequency <- rowMeans(totblocked, dims = 2) * 100 # frequency of Instantaneous Blocking days
  frequency2 <- rowMeans(totblocked2, dims = 2) * 100 # frequency of Instantaneous Blocking days with GHGS2
  U500mean <- rowMeans(U500, dims = 2) # U500 mean value

  t1 <- proc.time() - t0
  print(t1)

  print("Instantaneous blocking and diagnostics done!")

  ##########################################################
  #--------------------Time filtering----------------------#
  ##########################################################

  # spatial filtering on fixed longitude distance
  spatial <- longitude.filter(ics, ipsilon, totblocked)
  # CUT=apply(spatial,c(1,2),sum,na.rm=T)/ndays*100

  # large scale extension on 10x5 box
  large <- largescale.extension.if(ics, ipsilon, spatial)
  # LARGE=apply(large,c(1,2),sum,na.rm=T)/ndays*100

  # 5-day persistence filter
  block <- blocking.persistence(large, minduration = 5, time.array = etime)

  # 10-day persistence for extreme long block
  longblock <- blocking.persistence(large, minduration = 10, time.array = etime)

  tf <- proc.time() - t1
  print(tf)


  ##########################################################
  #------------------------Save to NetCDF------------------#
  ##########################################################

  # saving output to netcdf files
  print("saving NetCDF climatologies...")

  # which fieds to plot/save
  savelist <- c("TM90", "InstBlock", "ExtraBlock", "U500", "BlockEvents", "LongBlockEvents", "DurationEvents", "NumberEvents")
  full_savelist <- c("TM90", "InstBlock", "ExtraBlock", "U500", "BlockEvents", "LongBlockEvents")

  # dimension definition (using default 1850-01-01 reftime)
  dims <- ncdf.defdims(ics, ipsilon, timeaxis)

  # pre-declare list for loading the variable definition and the fields
  nc_field <- nc_var <- sapply(savelist, function(x) NULL)
  nc_fullfield <- nc_fullvar <- sapply(full_savelist, function(x) NULL)

  # loop on vars to save
  for (var in union(savelist, full_savelist)) {

    # default dimensions for climatological and full file
    dim_ncdf <- list(dims$x, dims$y, dims$z, t = dims$tclim)
    dim_fullncdf <- list(dims$x, dims$y, dims$z, t = dims$t)

    # name of the var
    if (var == "TM90") {
      dim_ncdf <- list(dims$x, t = dims$tclim)
      dim_fullncdf <- list(dims$x, t = dims$t)
      longvar <- "Tibaldi-Molteni 1990 Instantaneous Blocking frequency"
      unit <- "%"
      field <- TM90
      full_field <- totTM90
    }
    if (var == "InstBlock") {
      longvar <- "Instantaneous Blocking frequency"
      unit <- "%"
      field <- frequency
      full_field <- totblocked
    }
    if (var == "ExtraBlock") {
      longvar <- "Instantaneous Blocking frequency (GHGS2)"
      unit <- "%"
      field <- frequency2
      full_field <- totblocked2
    }
    if (var == "U500") {
      longvar <- "Zonal wind"
      unit <- "m/s"
      field <- U500mean
      full_field <- U500
    }
    if (var == "BlockEvents") {
      longvar <- "Blocking Events frequency"
      unit <- "%"
      field <- block$percentage
      full_field <- block$track
    }
    if (var == "LongBlockEvents") {
      longvar <- "10-day Blocking Events frequency"
      unit <- "%"
      field <- longblock$percentage
      full_field <- longblock$track
    }
    if (var == "DurationEvents") {
      longvar <- "Blocking Events duration"
      unit <- "days"
      field <- block$duration
    }
    if (var == "NumberEvents") {
      longvar <- "Blocking Events number"
      unit <- ""
      field <- block$nevents
    }

    # fix eventual NaN
    field[is.nan(field)] <- NA

    # define the list for variables and fields
    if (var %in% savelist) {
      nc_var[[which(var == savelist)]] <- ncvar_def(var, unit, dim_ncdf, -999,
        longname = longvar, prec = "single",
        compression = 1
      )
      nc_field[[which(var == savelist)]] <- field
    }

    if (var %in% full_savelist) {
      nc_fullvar[[which(var == full_savelist)]] <- ncvar_def(var, unit, dim_fullncdf, -999,
        longname = longvar, prec = "single",
        compression = 1
      )
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
    miles.block.fast(project, dataset, expid, ens, year1, year2, season, z500filename, FILESDIR, doforce)
  }
}
