######################################################
#-----Blocking routines computation for MiLES--------#
#-------------P. Davini (2014-2019)------------------#
######################################################
miles.block.multiple <- function(project, dataset, expid, ens, year1, year2, season, z500filename, 
                                 FILESDIR, doforce, biasccorrect, biasfile) {

  # this will track 3 different blocking indices producing 3 different files
  blocking_indices <- c("D12", "ExtraD12", "S04")

  # t0
  tstart <- proc.time()

  # setting up time domain
  years <- year1:year2
  timeseason <- season2timeseason(season)

  # brand new loop for the different tracking
  for (tracking_index in blocking_indices) {
    t0 <- proc.time()

    # name for bias correction
    if (biascorrect) {
      block_name <- "BC_Block"
    } else {
      block_name <- "Block"
    }

    # define folders using file.builder function (takes care of ensembles)
    savefile1 <- file.builder(FILESDIR, block_name, paste0(tracking_index, "_Clim"), project, dataset, expid, ens, year1, year2, season)
    savefile2 <- file.builder(FILESDIR, block_name, paste0(tracking_index, "_Full"), project, dataset, expid, ens, year1, year2, season)

    # check if data is already there to avoid re-run
    if (file.exists(savefile1) & file.exists(savefile2)) {
      print("Actually requested blocking data is already there!")
      print(savefile1)
      print(savefile2)
      if (as.logical(doforce)) {
        print("Running with doforce=true... re-run!")
      } else {
        print("Skipping... activate doforce=true if you want to re-run it")
        next
      }
    }

    # new file opening
    fieldlist <- ncdf.opener.universal(z500filename, namevar = "zg", tmonths = timeseason, tyears = years, rotate = "full", exportlonlat = F, verbose =F)
    print(str(fieldlist))

    # assign variables
    lon <- fieldlist$lon
    lat <- fieldlist$lat

    # extract calendar and time unit from the original file
    timeaxis <- fieldlist$time

    # declare variable and clean
    Z500 <- fieldlist$field
    rm(fieldlist)

    # time array to simplify time filtering
    etime <- power.date.new(timeaxis, verbose = T)
    totdays <- length(timeaxis)

    # option for bias correction
    if (biascorrect) {
      print("--------")
      print("Remove seasonal cycle and correcting it with observations...")
      # load a dataset to define the seasonal cycle
      biaslist <- ncdf.opener.universal(biasfile, namevar = "zg", tmonths = timeseason, tyears = years, rotate = "full", exportlonlat = F)
      biasetime <- power.date.new(biaslist$time, verbose = T)

      # clean records which are not included, so that they have the same length
      # this assumes that the observation has the same time span of the model data
      biaslist$field <- biaslist$field[,,biasetime$data %in% etime$data]
      biaslist$time <- biaslist$time[biasetime$data %in% etime$data]
      print("--------")
      print(str(biaslist))
      Z500 <- daily.replace.mean(lon, lat, Z500, biaslist$field, etime)
    }

    # grid resolution
    yreso <- lat[2] - lat[1]
    xreso <- lon[2] - lon[1]

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

    # TM90 and D98: parametres for blocking detection
    tm90_fi0 <- 60 # central_lat
    tm90_fiN <- tm90_fi0 + 20
    tm90_fiS <- tm90_fi0 - 20 # south and north lat, 80N and 40N
    tm90_central <- whicher(lat, tm90_fi0)
    tm90_south <- whicher(lat, tm90_fiS)
    tm90_north <- whicher(lat, tm90_fiN)
    tm90_range <- seq(-5, 5, yreso) / yreso # 5 degrees to the north, 5 to the south (larger than original TM90 or D'Andrea et al 1998 due to resolution)

    # 1D meridional gradients
    tm90_ghgn <- (Z500[, tm90_north + tm90_range, ] - Z500[, tm90_central + tm90_range, ]) / (tm90_fiN - tm90_fi0)
    tm90_ghgs <- (Z500[, tm90_central + tm90_range, ] - Z500[, tm90_south + tm90_range, ]) / (tm90_fi0 - tm90_fiS)

    print("Tibaldi and Molteni (1990) index...")
    tm90_check <- (tm90_ghgs > 0 & tm90_ghgn < (-10)) # TM90 conditions
    tm90_check[tm90_check == T] <- 1
    tm90_check[tm90_check == F] <- 0
    totTM90 <- apply(tm90_check, c(1, 3), max, na.rm = T)
    TM90 <- apply(totTM90, 1, mean) * 100

    print("D'Andrea et al. (1998) index...")
    d98_check <- (tm90_ghgs > 0 & tm90_ghgn < (-5)) # D98 conditions
    d98_check[d98_check == T] <- 1
    d98_check[d98_check == F] <- 0
    totD98 <- apply(d98_check, c(1, 3), max, na.rm = T)
    D98 <- apply(totD98, 1, mean) * 100


    ##########################################################
    #--------------Davini et al. 2012------------------------#
    ##########################################################

    # decleare main variables to be computed (considerable speed up!)
    totrwb <- totmeridional <- totBI <- Z500 * NA

    # Davini et al. 2012: parameters to be set for blocking detection
    fi0 <- 30 # lowest latitude to be analyzed
    jump <- 15 # distance on which compute gradients
    step0 <- jump / yreso # number of grid points to be used
    central <- which.min(abs(lat - fi0)) # lowest starting latitude
    north <- central + step0 # lowest north latitude
    south <- central - step0 # lowest sourth latitude
    maxsouth <- central - 2 * step0
    fiN <- lat[north]
    fiS <- lat[south]
    delta <- 0:((90 - fi0 - jump) / yreso) # escursion to the north for computing blocking (from 30 up to 75)

    print("--------------------------------------------------")
    print("Blocking indices...")

    ##########################################################
    #--------------Istantaneous Blocking---------------------#
    ##########################################################

    # Since Summer 2019 everything but the diagnostlon is vectorized: we didn't get
    # step up in speed but the code is much more compact

    print("gradients")
    # gradients for reversal index for D12 indices
    ghgn <- (Z500[, north + delta, ] - Z500[, central + delta, ]) / (fiN - fi0)
    ghgs <- (Z500[, central + delta, ] - Z500[, south + delta, ]) / (fi0 - fiS)
    gh2gs <- (Z500[, south + delta, ] - Z500[, maxsouth + delta, ]) / (fi0 - fiS)

    print("anomalies and threshold...")
    # daily anomalies and 90% for S04 absolute indices
    Z500anom <- daily.anom.mean(lon, lat, Z500, etime)
    threshold <- quantile(Z500anom[, whicher(lat, 50):whicher(lat, 80), ], probs = 0.9)
    Z500block <- Z500anom[, central + delta, ]

    # looop on the blocking indices
    for (idx in blocking_indices) {
      print(paste("Computing", idx, "instantaneous index ..."))

      if (idx == "D12") {
        # standard TM90 condition
        full_check <- which(ghgs > 0 & ghgn < (-10), arr.ind = T)
      } else if (idx == "ExtraD12") {
        # extra south condition to remove LLB
        full_check <- which(ghgs > 0 & ghgn < (-10) & gh2gs < (-5), arr.ind = T)
      } else if (idx == "S04") {
        # values exceeding 90% percentile of Z500
        full_check <- which(Z500block > threshold, arr.ind = T)
      }

      # if we detected any block (maybe remove, it is always happening on full dataset)
      if (length(full_check) > 0) {

        # create adjusted array
        full_repcheck <- cbind(full_check[, 1], full_check[, 2] + which(lat == fi0) - 1, full_check[, 3])

        # define blocking array and replace with detected instances
        totblocked <- Z500 * 0
        totblocked[full_repcheck] <- 1

        # assign to the corresponding index
        assign(paste0(idx, "_totblocked"), totblocked)

        # if the index is the one chosen for the tracking
        if (idx == tracking_index) {

          # to evaluate the diagnostic on the tracked grid point
          repcheck <- full_repcheck
          check <- full_check

          # for tracking
          totblocked_tracking <- totblocked
        }
      }
    }

    ##########################################################
    #--------------Blocking diagnostics ---------------------#
    ##########################################################

    print("Diagnostics...")

    # Creating a big abinded dataset
    print(paste("Creating abind dataset ..."))
    library("abind")
    new_field <- abind(Z500, Z500, Z500, along = 1)

    # 1- RWB orientation
    print(paste("Computing Rossby Wave Breaking..."))
    rwb_jump <- jump / 2
    steprwb <- rwb_jump / xreso

    # shift of 360 deg to estimate the longitudes
    ii <- repcheck[, 1] + length(lon)

    # computes the horizontal gradient right below the blocked point
    rwb_west <- new_field[cbind(ii - steprwb, repcheck[, 2] - steprwb, repcheck[, 3])]
    rwb_east <- new_field[cbind(ii + steprwb, repcheck[, 2] - steprwb, repcheck[, 3])]
    fullgh <- (rwb_west - rwb_east)

    # detect cyclonic or anticyclonic
    totrwb[repcheck[which(fullgh < 0), ]] <- (-10) # gradient decreasing: cyclonic RWB
    totrwb[repcheck[which(fullgh > 0), ]] <- 10 # gradient increasing: anticyclonic RWB

    # 2 -Wiedemann et al. 2002 blocking intensity: it's the slower part of the code
    print(paste("Computing Blocking Intensity..."))
    step <- 60 / xreso

    # minimum geopotential height 60 deg to the West
    zd <- sapply(1:length(ii), function(x) {
      left <- (ii[x] - step):ii[x]
      min(new_field[cbind(left, repcheck[x, 2], repcheck[x, 3])])
    })

    # minimum geopotential height 60 deg to the East
    zu <- sapply(1:length(ii), function(x) {
      right <- ii[x]:(ii[x] + step)
      min(new_field[cbind(right, repcheck[x, 2], repcheck[x, 3])])
    })

    # formula to get the intensity index
    mz <- Z500[repcheck]
    rc <- 0.5 * ((zu + mz) / 2 + (zd + mz) / 2)
    totBI[repcheck] <- 100 * (mz / rc - 1)

    # 3 - part about meridional gradient index
    print(paste("Computing Meridional Gradient Intensity..."))
    totmeridional[repcheck] <- ghgs[check]

    ##########################################################
    #--------------------Mean Values-------------------------#
    ##########################################################

    # compute mean values (use rowMeans that is faster when there are no NA values)
    for (idx in blocking_indices) {
      frequency <- rowMeans(get(paste0(idx, "_totblocked")), dims = 2) * 100 # frequency of Instantaneous Blocking days
      assign(paste0(idx, "_frequency"), frequency)
    }

    Z500mean <- rowMeans(Z500, dims = 2) # Z500 mean value
    BI <- apply(totBI, c(1, 2), mean, na.rm = T) # Blocking Intensity Index as Wiedenmann et al. (2002)
    MGI <- apply(totmeridional, c(1, 2), mean, na.rm = T) # Value of meridional gradient inversion

    # anticyclonic and cyclonic averages RWB
    CN <- apply(totrwb, c(1, 2), function(x) sum(x[x == (-10)], na.rm = T)) / (totdays) * (-10)
    ACN <- apply(totrwb, c(1, 2), function(x) sum(x[x == (10)], na.rm = T)) / (totdays) * (10)

    t1 <- proc.time() - t0
    print(t1)

    print("Instantaneous blocking and diagnostics done!")
    print(paste("Tracking", idx, "index ..."))

    ##########################################################
    #--------------------Time filtering----------------------#
    ##########################################################

    # spatial filtering on fixed longitude distance
    spatial <- longitude.filter(lon, lat, totblocked_tracking)
    # Deprecated: CUT=apply(spatial,c(1,2),sum,na.rm=T)/ndays*100

    # large scale extension on 10x5 box
    large <- largescale.extension.if(lon, lat, spatial)
    # Deprecated: LARGE=apply(large,c(1,2),sum,na.rm=T)/ndays*100

    # 5-day persistence filter
    block <- blocking.persistence(large, minduration = 5, time.array = etime)

    # 10-day persistence for extreme long block - deprecated since it worthless
    longblock <- blocking.persistence(large, minduration = 10, time.array = etime)

    t2 <- proc.time() - t1
    print(t2)


    ##########################################################
    #------------------------Save to NetCDF------------------#
    ##########################################################

    # saving output to netcdf files
    print("saving NetCDF climatologies...")

    # which fieds to plot/save
    savelist <- c("TM90", "D98", "InstBlock", "AbsBlock", "ExtraBlock", "Z500", "MGI", "BI", "CN", "ACN", "BlockEvents", "DurationEvents", "NumberEvents", "LongBlockEvents")
    full_savelist <- c("TM90", "D98", "InstBlock", "AbsBlock", "ExtraBlock", "Z500", "MGI", "BI", "CN", "ACN", "BlockEvents")

    # dimension definition (using default 1850-01-01 reftime)
    dims <- ncdf.defdims(lon, lat, timeaxis)

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
      } else if (var == "D98") {
        dim_ncdf <- list(dims$x, t = dims$tclim)
        dim_fullncdf <- list(dims$x, t = dims$t)
        longvar <- "D'Andrea et al 1998 Instantaneous Blocking frequency"
        unit <- "%"
        field <- D98
        full_field <- totD98
      } else if (var == "InstBlock") {
        longvar <- "Davini et al 2012 Instantaneous Blocking frequency"
        unit <- "%"
        field <- D12_frequency
        full_field <- D12_totblocked
      } else if (var == "AbsBlock") {
        longvar <- "Schwierz et al 2004 Instantaneous Blocking frequency"
        unit <- "%"
        field <- S04_frequency
        full_field <- S04_totblocked
      } else if (var == "ExtraBlock") {
        longvar <- "Davini et al 2012 Instantaneous Blocking frequency (GHGS2)"
        unit <- "%"
        field <- ExtraD12_frequency
        full_field <- ExtraD12_totblocked
      } else if (var == "Z500") {
        longvar <- "Geopotential Height"
        unit <- "m"
        field <- Z500mean
        full_field <- Z500
      } else if (var == "BI") {
        longvar <- "BI index"
        unit <- ""
        field <- BI
        full_field <- totBI
      } else if (var == "MGI") {
        longvar <- "MGI index"
        unit <- ""
        field <- MGI
        full_field <- totmeridional
      } else if (var == "ACN") {
        longvar <- "Anticyclonic RWB frequency"
        unit <- "%"
        field <- ACN
        full_field <- totrwb / 10
        full_field[full_field == (-1)] <- NA
      } else if (var == "CN") {
        longvar <- "Cyclonic RWB frequency"
        unit <- "%"
        field <- CN
        full_field <- totrwb / 10
        full_field[full_field == (1)] <- NA
      } else if (var == "BlockEvents") {
        longvar <- "Blocking Events frequency"
        unit <- "%"
        field <- block$percentage
        full_field <- block$track
      } else if (var == "LongBlockEvents") {
        longvar <- "10-day Blocking Events frequency"
        unit <- "%"
        field <- longblock$percentage
        full_field <- longblock$track
      } else if (var == "DurationEvents") {
        longvar <- "Blocking Events duration"
        unit <- "days"
        field <- block$duration
      } else if (var == "NumberEvents") {
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
    t3 <- proc.time() - t2
    print(t3)
  }

  tfinal <- proc.time() - tstart
  print(tfinal)
}

# blank lines
cat("\n\n\n")

# REAL EXECUTION OF THE SCRIPT
# read command line
args <- commandArgs(TRUE)

# number of required arguments from command line
name_args <- c("project", "dataset", "expid", "ens", "year1", "year2", "season", "z500filename", 
               "FILESDIR", "PROGDIR", "doforce", "biascorrect", "biasfile")

# if there arguments, check them required args and assign
if (length(args) != 0) {
  req_args <- length(name_args)
  if (length(args) != req_args) {

    #    # stop if something is wrong
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
    miles.block.multiple(project, dataset, expid, ens, year1, year2, season, z500filename, FILESDIR, doforce, biascorrect, biasfile)
  }
}
