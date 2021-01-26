##########################################################
#-------------------Filtering functions------------------#
##########################################################

fourier.filter <- function(timeseries, period, type = "nofilter") {
  # Basic timeseries filtering based on ffts
  # It works on a given time period which is assumed to work on the
  # time window of your timeseris (i.e. if you have daily data, "period"
  # should provide the period you want cut/keep)
  # type can be "bandpass", "highpass" and "lowpass"

  wavenumber <- rev(length(timeseries) / period)

  fourier <- fft(timeseries)
  magnitude <- Mod(fourier)
  phase <- Arg(fourier)

  # filter_magnitude=magnitude[1:(length(magnitude)/2+1)]
  filter_magnitude <- magnitude[1:floor((length(magnitude) / 2 + 1))]
  # print(length(filter_magnitude))
  zeroaxis <- (0:(length(filter_magnitude) - 1))

  if (type == "nofilter") {
    print("No filtering")
  }
  if (type == "bandpass") {
    filter_magnitude[zeroaxis < wavenumber[1] | zeroaxis > wavenumber[2]] <- 0
  }
  if (type == "lowpass") {
    filter_magnitude[zeroaxis > wavenumber[1]] <- 0
  }
  if (type == "highpass") {
    filter_magnitude[zeroaxis < wavenumber[1]] <- 0
  }

  # new_magnitude=c(filter_magnitude,rev(filter_magnitude[2:(length(filter_magnitude)-1)]))
  if ((length(timeseries) %% 2) == 1) {
    new_magnitude <- c(filter_magnitude, rev(filter_magnitude[2:(length(filter_magnitude))]))
  } else {
    new_magnitude <- c(filter_magnitude, rev(filter_magnitude[2:(length(filter_magnitude) - 1)]))
  }

  filter_phase <- phase
  filter_phase[which(new_magnitude == 0)] <- 0

  filter_fourier <- complex(length.out = length(new_magnitude), modulus = new_magnitude, argument = filter_phase)
  filter_y <- fft(filter_fourier, inverse = TRUE)
  filter_y <- Re(filter_y / length(filter_y))


  return(filter_y)
}

# create lanczos weights for filtering
# need to be applied with stats::filter function
# references: https://www2.atmos.umd.edu/~ekalnay/syllabi/AOSC630/METO630ClassNotes13.pdf
# original code: https://stackoverflow.com/questions/17264119/using-lanczos-low-pass-filter-in-r-program
lanczos_weights <- function(nw = 11, fc = 1 / 30, type = "lowpass") {

  # create weights
  n <- nw %/% 2
  k <- (-n):n
  w <- sin(pi * k / n) / (pi * k / n) * sin(2 * pi * fc * k) / (pi * k)

  if (type == "lowpass") {
    out <- w
    out[n + 1] <- 2 * fc
  } else if (type == "highpass") {
    out <- -w
    out[n + 1] <- 1 - 2 * fc
  }
  return(out)
}

# lanczos filtering: use lanczos weight to create a different filter.
# it is about 50% faster than fourier filtering thanks to convolution
lanczos_filter <- function(timeseries, period, window, type = "nofilter") {
  if (type == "lowpass") {
    weights <- lanczos_weights(window, 1 / period, type)
  } else if (type == "highpass") {
    weights <- lanczos_weights(window, 1 / period, type)
  } else if (type == "bandpass") {
    weights <- lanczos_weights(window, 1 / period, "lowpass") + lanczos_weights(window, 1 / period, "highpass")
  } else if (type == "nofilter") {
    print("No filtering")
  }

  filter <- stats::filter(timeseries, weights)
  return(filter)
}
