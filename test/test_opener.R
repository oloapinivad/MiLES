#zero order test to check integrity of netcdf opener

rm(list=ls())
PROGDIR="/home/paolo/MiLES"
FILESDIR="/work/users/paolo/miles/files"
source(file.path(PROGDIR,"script/basis_functions.R"))
library("SDMTools")


filenames=c("/work/users/paolo/miles/data/zg500/ERAI/zg500_ERAI_fullfile.nc")
season=c("DJF")

for (filename in filenames)  {
  year1=2000
  year2=2010

  #setting up time domain
  years=year1:year2
  timeseason=season2timeseason(season)

  #load file
  fieldlist=ncdf.opener.universal(filename,namevar="zg",tmonths=timeseason,tyears=years,rotate="full",verbose=T)

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

  print("Schwierz et al. (2004) index...")
  # compute anomalies
  Z500anom <- daily.anom.mean(ics, ipsilon, Z500, etime)

  # threshold definition
  threshold <- quantile(Z500anom[, whicher(ipsilon,50):whicher(ipsilon,80), ], probs = 0.9)

  # defining blocking
  absblocked <- Z500anom
  absblocked[absblocked < threshold] <- 0
  absblocked[absblocked >= threshold] <- 1

  filtered <- absblocked * NA
  for (i in 1:2) {

    # boundary extension: overcome issues of periodic boundaries
    # extensions is adaptable, it tries to to be the cheapest as possible
    firstzero <- which(apply(absblocked[,,i], c(1) ,sum, na.rm=T) == 0)[1]
    binded <- rotation(absblocked[,,i], length(ics) - firstzero)

    # identify connected-component labeling after re-rotate
    a <- rotation(ConnCompLabel(binded), firstzero) 

    # remove components which have less than 10 points
    npointsmin <- 10
    a[a==which(table(a, exclude=0) < npointsmin)] <- 0
    a[a!=0] <- 1

    filtered[,,i] <- a

  }

  # climatology
  absfrequency <- rowMeans(absblocked, dims = 2) * 100
  print("Done!")

}


