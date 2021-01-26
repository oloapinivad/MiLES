# source for routines

# load packages
library("maps")
library("ncdf4")
library("PCICt")

# extra trick to use older definition on "MILESDIR"
# this is a configuration used to call MiLES environment outside MiLES program
# plenty of the routines in MiLES are not used by MiLES itself and they are simply stored here
# if you want to use the full potential of the MiLES routines set MILESDIR before calling this script
if (!exists("PROGDIR")) { 
  ROUTDIR <- file.path(MILESDIR, "script", "routines")
  routines <- list.files(ROUTDIR)
} else {
  # list routines but the name of this file
  ROUTDIR <- file.path(PROGDIR, "script", "routines")
  routines <- c("array_manipulation.R", "constants.R", "blocking.R", "eofs_regimes.R",
                "general.R", "miles.R", "ncdf_opening.R", "ncdf_writing.R", 
                "plotting_functions.R", "spatial_averaging.R",
                "time_manipulation.R", "time_averaging.R")
}

# loading the routine files
for (routine in routines) {
  print(paste("Loading", routine))
  source(file.path(ROUTDIR,routine))
}
