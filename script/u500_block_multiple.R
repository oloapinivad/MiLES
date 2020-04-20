######################################################
#-Zonal wind blocking routines computation for MiLES-#
#-------------P. Davini (Dec 2018)-------------------#
######################################################
miles.u500block.multiple <- function(project, dataset, expid, ens, year1, year2, season, filename, FILESDIR, doforce) {

  # this will track 3 different blocking indices producing 3 different files
  blocking_indices <- c("D12", "ExtraD12")

  # t0
  tstart <- proc.time()

  # setting up time domain
  years <- year1:year2
  timeseason <- season2timeseason(season)

  for (tracking_index in blocking_indices) {
    t0 <- proc.time()

    # define folders using file.builder function (takes care of ensembles)
    savefile1 <- file.builder(FILESDIR, "U500_Block", paste0(tracking_index, "_Clim"), project, dataset, expid, ens, year1, year2, season)
    savefile2 <- file.builder(FILESDIR, "U500_Block", paste0(tracking_index, "_Full"), project, dataset, expid, ens, year1, year2, season)

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
    fieldlist <- ncdf.opener.universal(filename, namevar = "ua", tmonths = timeseason, tyears = years, rotate = "full", exportlonlat = F)
    print(str(fieldlist))

    # assign variables
    lon <- fieldlist$lon
    lat <- fieldlist$lat

    # extract calendar and time unit from the original file
    timeaxis <- fieldlist$time

    # time array to simplify time filtering
    etime <- power.date.new(timeaxis, verbose = T)
    totdays <- length(timeaxis)

    # declare variable
    U500 <- fieldlist$field
    rm(fieldlist)

    # grid resolution
    yreso <- lat[2] - lat[1]
    xreso <- lon[2] - lon[1]

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
    tm90_central <- whicher(lat, tm90_fi0)
    tm90_south <- whicher(lat, tm90_fiS)
    tm90_north <- whicher(lat, tm90_fiN)
    tm90_range <- seq(-5, 5, yreso) / yreso # 5 degrees to the north, 5 to the south (larger than TM90 or D'Andrea et al 1998)

    # from u500: define parameters for geostrophic approximation (constant are defined in basis_functions.R)
    sinphi <- sin(lat * pi / 180)
    alfa <- (2 * Earth.Radius * omega / g0) * ((tm90_delta * pi / 180))
    ghgn <- ghgs <- array(NA, dim = c(length(lon), totdays, length(tm90_range)))

    # number of steps and weights for integrals
    tm90_step <- round(tm90_delta / yreso)
    tm90_ww <- c(0.5, rep(1, tm90_step - 1), 0.5)

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
    print("Tibaldi and Molteni (1990) index...")
    tm90_check <- (ghgs < 0 & ghgn > (10 * tm90_delta / alfa)) # TM90 adapted to Scaife et al. 2010
    tm90_check[tm90_check == T] <- 1
    tm90_check[tm90_check == F] <- 0
    totTM90 <- apply(tm90_check, c(1, 2), max, na.rm = T)
    TM90 <- apply(totTM90, 1, mean) * 100

    print("D'Andrea et al. (1998) index...")
    d98_check <- (ghgs < 0 & ghgn > (5 * tm90_delta / alfa)) # D98 conditions adapted to Scaife et al. 2010
    d98_check[d98_check == T] <- 1
    d98_check[d98_check == F] <- 0
    totD98 <- apply(d98_check, c(1, 2), max, na.rm = T)
    D98 <- apply(totD98, 1, mean) * 100

    ##########################################################
    #--------------Davini et al. 2012------------------------#
    ##########################################################

    # Davini et al. 2012: parameters to be set for blocking detection
    fi0 <- 30 # lowest latitude to be analyzed
    delta <- 15 # distance on which compute gradients
    step <- round(delta / yreso)
    central <- which.min(abs(lat - fi0)) # lowest starting latitude
    north <- central + step # lowest north latitude
    south <- central - step # lowest sourth latitude
    maxsouth <- central - 2 * step
    fiN <- lat[north]
    fiS <- lat[south]
    range <- (90 - fi0 - delta) / yreso # escursion to the north for computing blocking (from 30 up to 75)

    # number of steps and weights for integrals
    ww <- c(0.5, rep(1, step - 1), 0.5)

    print("--------------------------------------------------")
    print("Zonal wind Davini et al. (2012) index and diagnostlon...")
    print(c("distance for gradients:", step * yreso))
    print(paste("range of latitudes ", fi0, "-", 90 - step * yreso, " N", sep = ""))

    ##########################################################
    #--------------Istantaneous Blocking---------------------#
    ##########################################################

    for (idx in blocking_indices) {
      # decleare main variables to be computed (considerable speed up!)
      totblocked <- U500 * 0


      #----COMPUTING BLOCKING INDICES-----
      for (t in 1:totdays) {
        progression.bar(t, totdays)

        ghgn <- ghgs <- gh2gs <- array(NA, dim = c(length(lon), length(0:range)))

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


          if (idx == "D12") {
            check1 <- (ghgs < 0 & ghgn > (10 * delta / alfa))
            check1[check1 == T] <- 1
            check1[check1 == F] <- 0
            totblocked[, central:(central + range), t] <- check1
          } else if (idx == "ExtraD12") {
            check1 <- (ghgs < 0 & ghgn > (10 * delta / alfa) & gh2gs > (5 * delta / alfa))
            check1[check1 == T] <- 1
            check1[check1 == F] <- 0
            totblocked[, central:(central + range), t] <- check1
          }
        }
      }
      if (tracking_index == idx) {
        totblocked_tracking <- totblocked 
      }
      assign(paste0(idx, "_totblocked"), totblocked)
      print(paste("Total # of days:", t))
      print("-------------------------")
    }

    ##########################################################
    #--------------------Mean Values-------------------------#
    ##########################################################

    # compute mean values (use rowMeans that is faster when there are no NA values)
    for (idx in blocking_indices) {
      frequency <- rowMeans(get(paste0(idx, "_totblocked")), dims = 2) * 100 # frequency of Instantaneous Blocking days
      assign(paste0(idx, "_frequency"), frequency)
    }

    U500mean <- rowMeans(U500, dims = 2) # U500 mean value

    t1 <- proc.time() - t0
    print(t1)

    print("Instantaneous blocking and diagnostlon done!")

    ##########################################################
    #--------------------Time filtering----------------------#
    ##########################################################

    # spatial filtering on fixed longitude distance
    spatial <- longitude.filter(lon, lat, totblocked_tracking)
    # CUT=apply(spatial,c(1,2),sum,na.rm=T)/ndays*100

    # large scale extension on 10x5 box
    large <- largescale.extension.if(lon, lat, spatial)
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
    savelist <- c("TM90", "D98", "InstBlock", "ExtraBlock", "U500", "BlockEvents", "DurationEvents", "NumberEvents", "LongBlockEvents")
    full_savelist <- c("TM90", "D98", "InstBlock", "ExtraBlock", "U500", "BlockEvents")

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
        longvar <- "Instantaneous Blocking frequency"
        unit <- "%"
        field <- D12_frequency
        full_field <- D12_totblocked
      } else if (var == "ExtraBlock") {
        longvar <- "Instantaneous Blocking frequency (GHGS2)"
        unit <- "%"
        field <- ExtraD12_frequency
        full_field <- ExtraD12_totblocked
      } else if (var == "U500") {
        longvar <- "Zonal wind"
        unit <- "m/s"
        field <- U500mean
        full_field <- U500
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
}

# blank lines
cat("\n\n\n")

# REAL EXECUTION OF THE SCRIPT
# read command line
args <- commandArgs(TRUE)

# number of required arguments from command line
name_args <- c("project", "dataset", "expid", "ens", "year1", "year2", "season", "u500filename", "FILESDIR", "PROGDIR", "doforce")

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
    miles.u500block.multiple(project, dataset, expid, ens, year1, year2, season, u500filename, FILESDIR, doforce)
  }
}
