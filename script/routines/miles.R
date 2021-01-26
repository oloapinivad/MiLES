# miles specific functions

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

  if (field == "D98") {
    color_field <- c("navy", "darkorange")
    color_diff <- NULL
    lev_field <- c(0, 30)
    lev_diff <- NULL
    legend_unit <- "Blocked Days (%)"
    title_name <- "Instantaneous Blocking (D'Andrea et al. 1998):"
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

  out <- list(
    color_field = color_field, color_diff = color_diff, lev_field = lev_field,
    lev_diff = lev_diff, lev_hist = lev_hist, legend_unit = legend_unit,
    legend_distance = legend_distance, title_name = title_name
  )
  return(out)
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
