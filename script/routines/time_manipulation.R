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

# leap year treu/false function
is.leapyear <- function(year) {
  return(((year %% 4 == 0) & (year %% 100 != 0)) | (year %% 400 == 0))
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
  min_interval <- min(diff(datas))
  if (min_interval == 1) {
    print("daily data!")
    breakdays <- 1
  } else if (min_interval >= 28 & min_interval <= 31) {
    print("monthly data!")
    breakdays <- 31
  } else {
    stop("Unknonw time interval, breaking!")
  }

  # create a "season" for continuous time, used by persistance tracking and season average
  startpoints <- c(0, which(diff(datas) > breakdays))
  deltapoints <- diff(c(startpoints, length(datas)))
  seas <- inverse.rle(list(lengths = deltapoints, values = seq(1, length(startpoints))))

  # identify incomplete seasons to be filtered out (weird method based on math)
  season_length <- table(seas)
  incomplete_season <- which(season_length < (mean(season_length) - sd(season_length)))

  etime <- list(
    day = as.numeric(format(datas, "%d")), month = as.numeric(format(datas, "%m")),
    year = as.numeric(format(datas, "%Y")), data = datas, season = seas,
    incomplete_season = incomplete_season
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

  printv(paste("This is a", caldata, "calendar"), verbose)

  # new method including both absolute and releative time axis (Oct 2018)
  # if "as" is present, this is an absolute time axis
  if (grepl("as", units, fixed = TRUE)) {
    kind <- "absolute"
    freq <- NULL
    timeline <- as.PCICt(as.character(time), format = "%Y%m%d", cal = caldata)
    # if a "since" is present, this is a relative time axis
  } else if (grepl("since", units, fixed = TRUE)) {
    kind <- "relative"
    # origin <- substr(gsub("[a-zA-Z ]", "", units), 1, 10)
    origin <- unlist(strsplit(units, "[a-zA-Z ]+"))[2]
    origin.pcict <- as.PCICt(origin, cal = caldata, format = "%Y-%m-%d")
    printv(paste("Origin of the time axis is:", origin.pcict), verbose)

    # distinguish between day and seconds based axis
    if (grepl("day", substr(units, 1, 3), fixed = TRUE)) {
      freq <- "(days)"
      timeline <- origin.pcict + (floor(time) * 86400)
    } else if (grepl("sec", substr(units, 1, 3), fixed = TRUE)) {
      freq <- "(seconds)"
      timeline <- origin.pcict + floor(time)
    } else if (grepl("hour", substr(units, 1, 4), fixed = TRUE)) {
      freq <- "(hours)"
      timeline <- origin.pcict + floor(time) * 3600
    } else if (grepl("month", substr(units, 1, 5), fixed = TRUE)) {
      freq <- "(months)"
      timeline <- origin.pcict + floor(time) * 86400 * 30.43
    } else if (grepl("year", substr(units, 1, 4), fixed = TRUE)) {
      freq <- "(years)"
      timeline <- origin.pcict + floor(time) * 86400 * 365
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
    kind <- "unsupported"
  }
  printv(paste("Calendar:", kind, "time axis", freq), verbose)


  return(list(timeline = timeline, calendar = kind))
}
