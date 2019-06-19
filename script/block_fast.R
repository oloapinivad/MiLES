######################################################
#-----Blocking routines computation for MiLES--------#
#-------------P. Davini (Oct 2014)-------------------#
######################################################
miles.block.fast <- function(project, dataset, expid, ens, year1, year2, season, z500filename, FILESDIR, doforce) {

  # t0
  t0 <- proc.time()

  # setting up time domain
  years <- year1:year2
  timeseason <- season2timeseason(season)

  # define folders using file.builder function (takes care of ensembles)
  savefile1 <- file.builder(FILESDIR, "Block", "BlockClim", project, dataset, expid, ens, year1, year2, season)
  savefile2 <- file.builder(FILESDIR, "Block", "BlockFull", project, dataset, expid, ens, year1, year2, season)

  # check if data is already there to avoid re-run
  if (file.exists(savefile1) & file.exists(savefile2)) {
    print("Actually requested blocking data is already there!")
    print(savefile1)
    print(savefile2)
    if (as.logical(doforce)) {
      print("Running with doforce=true... re-run!")
    } else {
      print("Skipping... activate doforce=true if you want to re-run it")
      return() #q()
    }
  }

  # new file opening
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

  # grid resolution
  yreso <- ipsilon[2] - ipsilon[1]
  xreso <- ics[2] - ics[1]

  # reso checks: this are not needed with default 2.5 grid, but they may be relevant with
  # future envisaged power up to finer grids
  # xcritical factor is due to RWB longitudinal jump of 7.5
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

  print("Tibaldi and Molteni (1990) index...")
  # TM90: parametres for blocking detection
  tm90_fi0 <- 60 # central_lat
  tm90_fiN <- tm90_fi0 + 20
  tm90_fiS <- tm90_fi0 - 20 # south and north lat, 80N and 40N
  tm90_central <- whicher(ipsilon, tm90_fi0)
  tm90_south <- whicher(ipsilon, tm90_fiS)
  tm90_north <- whicher(ipsilon, tm90_fiN)
  tm90_range <- seq(-5, 5, yreso) / yreso # 5 degrees to the north, 5 to the south (larger than TM90 or D'Andrea et al 1998)

  # TM90: beta version, the amazing power of R vectorization!
  # 6 lines to get the climatology
  tm90_ghgn <- (Z500[, tm90_north + tm90_range, ] - Z500[, tm90_central + tm90_range, ]) / (tm90_fiN - tm90_fi0)
  tm90_ghgs <- (Z500[, tm90_central + tm90_range, ] - Z500[, tm90_south + tm90_range, ]) / (tm90_fi0 - tm90_fiS)
  tm90_check <- (tm90_ghgs > 0 & tm90_ghgn < (-10)) # TM90 conditions
  tm90_check[tm90_check == T] <- 1
  tm90_check[tm90_check == F] <- 0
  totTM90 <- apply(tm90_check, c(1, 3), max, na.rm = T)
  TM90 <- apply(totTM90, 1, mean) * 100
  print("Done!")

  ##########################################################
  #--------------Schwierz et al. 2004----------------------#
  ##########################################################

  print("Schwierz et al. (2004) index...")
  # compute anomalies
  Z500anom <- daily.anom.mean(ics, ipsilon, Z500, etime)

  # threshold definition
  threshold <- quantile(Z500anom[, whicher(ipsilon,50):whicher(ipsilon,80), ], probs = 0.9)

  # defining blocking
  absblocked <- Z500anom
  absblocked[absblocked < threshold] <- 0
  absblocked[absblocked >= threshold] <- 1

  # climatology
  absfrequency <- rowMeans(absblocked, dims = 2) * 100
  print("Done!")

  ##########################################################
  #--------------Davini et al. 2012------------------------#
  ##########################################################

  # decleare main variables to be computed (considerable speed up!)
  totrwb <- totmeridional <- totBI <- Z500 * NA
  totblocked <- totblocked2 <- totblocked_test <- Z500 * 0

  # Davini et al. 2012: parameters to be set for blocking detection
  fi0 <- 30 # lowest latitude to be analyzed
  jump <- 15 # distance on which compute gradients
  step0 <- jump / yreso # number of grid points to be used
  central <- which.min(abs(ipsilon - fi0)) # lowest starting latitude
  north <- central + step0 # lowest north latitude
  south <- central - step0 # lowest sourth latitude
  maxsouth <- central - 2 * step0
  fiN <- ipsilon[north]
  fiS <- ipsilon[south]
  range <- (90 - fi0 - jump) / yreso # escursion to the north for computing blocking (from 30 up to 75)

  print("--------------------------------------------------")
  print("Davini et al. (2012) index and diagnostics...")
  print(c("distance for gradients:", step0 * diff(ics)[1]))
  print(paste("range of latitudes ", fi0, "-", 90 - step0 * diff(ics)[1], " N", sep = ""))

  ##########################################################
  #--------------Istantaneous Blocking---------------------#
  ##########################################################

  #----COMPUTING BLOCKING INDICES-----
  for (t in 1:totdays) {
    progression.bar(t, totdays)

    # multidim extension
    new_field <- rbind(Z500[, , t], Z500[, , t], Z500[, , t])

    for (delta in 0:range) { # computing blocking for different latitudes
      ghgn <- (Z500[, north + delta, t] - Z500[, central + delta, t]) / (fiN - fi0)
      ghgs <- (Z500[, central + delta, t] - Z500[, south + delta, t]) / (fi0 - fiS)
      gh2gs <- (Z500[, south + delta, t] - Z500[, maxsouth + delta, t]) / (fi0 - fiS)
      check1 <- which(ghgs > 0 & ghgn < (-10))
      check2 <- which(ghgs > 0 & ghgn < (-10) & gh2gs < (-5)) # supplementary condition
      check3 <- which(ghgs > 0 & ghgn < 0 & gh2gs < 0) # test condition

      if (length(check2) > 0) {
        totblocked2[check2, central + delta, t] <- 1
      }

      if (length(check3) > 0) {
        totblocked_test[check3, central + delta, t] <- 1
      }

      if (length(check1) > 0) {
        # 1-MATRIX FOR INSTANTANEOUS BLOCKING
        totblocked[check1, central + delta, t] <- 1

        # 2-PART ON COMPUTATION OF ROSSBY WAVEBREAKING
        r <- check1 + length(ics)
        rwb_jump <- jump / 2
        steprwb <- rwb_jump / xreso
        rwb_west <- new_field[(r - steprwb), south + delta + steprwb]
        rwb_east <- new_field[(r + steprwb), south + delta + steprwb]
        fullgh <- (rwb_west - rwb_east)

        totrwb[check1[fullgh < 0], central + delta, t] <- (-10) # gradient decreasing: cyclonic RWB
        totrwb[check1[fullgh > 0], central + delta, t] <- 10 # gradient increasing: anticyclonic RWB

        # 4-part about adapted version of blocking intensity by Wiedenmann et al. (2002)
        step <- 60 / xreso
        ii <- check1 + length(ics)
        zu <- zd <- NULL
        for (ll in ii) {
          zu <- c(zu, min(new_field[(ll - step):ll, central + delta]))
          zd <- c(zd, min(new_field[ll:(ll + step), central + delta]))
        }
        mz <- Z500[check1, central + delta, t]
        rc <- 0.5 * ((zu + mz) / 2 + (zd + mz) / 2)
        totBI[check1, central + delta, t] <- 100 * (mz / rc - 1)

        # 5 - part about meridional gradient index
        totmeridional[check1, central + delta, t] <- ghgs[check1]
      }
    }
  }

  print(paste("Total # of days:", t))
  print("-------------------------")

  ##########################################################
  #--------------------Mean Values-------------------------#
  ##########################################################

  # compute mean values (use rowMeans that is faster when there are no NA values)
  frequency <- rowMeans(totblocked, dims = 2) * 100 # frequency of Instantaneous Blocking days
  frequency2 <- rowMeans(totblocked2, dims = 2) * 100 # frequency of Instantaneous Blocking days with GHGS2
  frequency_test <- rowMeans(totblocked_test, dims = 2) * 100 # frequency of Instantaneous Blocking days test
  Z500mean <- rowMeans(Z500, dims = 2) # Z500 mean value
  BI <- apply(totBI, c(1, 2), mean, na.rm = T) # Blocking Intensity Index as Wiedenmann et al. (2002)
  MGI <- apply(totmeridional, c(1, 2), mean, na.rm = T) # Value of meridional gradient inversion

  # anticyclonic and cyclonic averages RWB
  CN <- apply(totrwb, c(1, 2), function(x) sum(x[x == (-10)], na.rm = T)) / (totdays) * (-10)
  ACN <- apply(totrwb, c(1, 2), function(x) sum(x[x == (10)], na.rm = T)) / (totdays) * (10)

  t1 <- proc.time() - t0
  print(t1)

  print("Instantaneous blocking and diagnostics done!")

  ##########################################################
  #--------------------Time filtering----------------------#
  ##########################################################

  # spatial filtering on fixed longitude distance
  spatial <- longitude.filter(ics, ipsilon, totblocked)
  # Deprecated: CUT=apply(spatial,c(1,2),sum,na.rm=T)/ndays*100

  # large scale extension on 10x5 box
  large <- largescale.extension.if(ics, ipsilon, spatial)
  # Deprecated: LARGE=apply(large,c(1,2),sum,na.rm=T)/ndays*100

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
  savelist <- c("TM90", "InstBlock", "AbsBlock", "ExtraBlock", "TestBlock", "Z500", "MGI", "BI", "CN", "ACN", "BlockEvents", "LongBlockEvents", "DurationEvents", "NumberEvents")
  full_savelist <- c("TM90", "InstBlock", "AbsBlock", "ExtraBlock", "TestBlock", "Z500", "MGI", "BI", "CN", "ACN", "BlockEvents", "LongBlockEvents")

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
    if (var == "AbsBlock") {
      longvar <- "Instantaneous Blocking frequency (Schwierz et al)"
      unit <- "%"
      field <- absfrequency
      full_field <- absblocked
    }
    if (var == "ExtraBlock") {
      longvar <- "Instantaneous Blocking frequency (GHGS2)"
      unit <- "%"
      field <- frequency2
      full_field <- totblocked2
    }
    if (var == "TestBlock") {
      longvar <- "Instantaneous Blocking frequency - TEST INDEX"
      unit <- "%"
      field <- frequency_test
      full_field <- totblocked_test
    }
    if (var == "Z500") {
      longvar <- "Geopotential Height"
      unit <- "m"
      field <- Z500mean
      full_field <- Z500
    }
    if (var == "BI") {
      longvar <- "BI index"
      unit <- ""
      field <- BI
      full_field <- totBI
    }
    if (var == "MGI") {
      longvar <- "MGI index"
      unit <- ""
      field <- MGI
      full_field <- totmeridional
    }
    if (var == "ACN") {
      longvar <- "Anticyclonic RWB frequency"
      unit <- "%"
      field <- ACN
      full_field <- totrwb / 10
      full_field[full_field == (-1)] <- NA
    }
    if (var == "CN") {
      longvar <- "Cyclonic RWB frequency"
      unit <- "%"
      field <- CN
      full_field <- totrwb / 10
      full_field[full_field == (1)] <- NA
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

    # variable definitions
    if (var %in% savelist) {
      nc_var[[which(var == savelist)]] <- ncvar_def(var, unit, dim_ncdf, -999, longname = longvar, prec = "single", compression = 1)
      nc_field[[which(var == savelist)]] <- field
    }

    # field definitions
    if (var %in% full_savelist) {
      nc_fullvar[[which(var == full_savelist)]] <- ncvar_def(var, unit, dim_fullncdf, -999, longname = longvar, prec = "single", compression = 1)
      nc_fullfield[[which(var == full_savelist)]] <- full_field
    }
  }

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
