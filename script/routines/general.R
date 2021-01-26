# check if fast linear fit is operative (after R 3.1): 3x faster than lm.fit, 36x faster than lm
if (exists(".lm.fit")) {
  lin.fit <- .lm.fit
} else {
  lin.fit <- lm.fit
}

# check R version as numeric
R_version <- as.numeric(R.Version()$major) + as.numeric(R.Version()$minor) / 10


# normalize a time series
standardize <- function(timeseries) {
  out <- (timeseries - mean(timeseries, na.rm = T)) / sd(timeseries, na.rm = T)
  return(out)
}

# file manipulation
append.nc <- function(filename) {
  return(paste0(filename, ".nc"))
}

# fast conversion of angles into radiants (and viceversa)
rad <- function(angle, rev = F) {
  if (rev) {
    r <- angle * 180 / pi
  } else {
    r <- angle * pi / 180
  }
  return(r)
}

# verbose-only printing function
printv <- function(value, verbosity = TRUE) {
  if (verbosity) {
    print(value)
  }
}

print.break <- function(n = 50) {
  return(paste(rep("-", 50), collapse = ""))
}

# convert pressure levels to altitude, using standard atmosphere
vertical <- function(plev, reverse = F) {
  if (max(plev) > 2000) {
    slp <- p0
  } else {
    slp <- 1000
  }
  if (reverse) {
    vertical <- 100 * slp * exp(-plev / 8)
  } else {
    vertical <- -8000 * log(plev / slp) / 1000
  }
  return(vertical)
}

# progression bar
progression.bar <- function(index, total_length, each = 10) {
  if (any(index == round(seq(0, total_length, , each + 1)))) {
    progression <- paste("--->", round(index / total_length * 100), "%")
    print(progression)
  }
}
