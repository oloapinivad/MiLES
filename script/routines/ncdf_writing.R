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
