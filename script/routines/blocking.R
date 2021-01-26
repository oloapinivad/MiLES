##########################################################
#------------Blocking Tracking Functions-----------------#
##########################################################

# time persistence (used for longitude filter too)
time.persistence <- function(timeseries, persistence = 5) {
  rr <- rle(timeseries)
  rr$values[which(rr$values == 1 & rr$length < persistence)] <- 0
  nn <- rep(rr$values, rr$length)
  return(nn)
}


# blocking 5 days tracking
blocking.persistence <- function(field, minduration = 5, time.array) {

  # function for persistence
  pers2 <- function(timeseries, persistence, time.array) {
    dd <- min(time.array$season):max(time.array$season)
    nn <- sapply(dd, function(x) {
      time.persistence(timeseries[which(time.array$season == x)], persistence)
    })
    xx <- c(unlist(nn))
    return(xx)
  }

  # check for etime
  if (length(time.array$month) != length(field[1, 1, ])) {
    stop("Wrong time array! Exiting...")
  }

  print("Time filtering...")
  newfield <- apply(field, c(1, 2), function(x) pers2(x, persistence = minduration, time.array))
  newfield <- aperm(newfield, c(2, 3, 1))
  print("Mean field...")
  meanfield <- apply(newfield, c(1, 2), mean, na.rm = T) * 100


  print("Events detection...")
  maxdim <- max(apply(newfield, c(1, 2), function(x) length(rle(x)$length[which(rle(x)$values == 1)])))
  events <- apply(newfield, c(1, 2), function(x) c(rle(x)$lengths[which(rle(x)$values == 1)], rep(NA, maxdim - length(rle(x)$length[which(rle(x)$values == 1)]))))
  if (length(dim(events)) == 2) {
    print("Few events, collapsing recognized")
    events <- array(events, dim = c(dim(events), 1))
  } else {
    events <- aperm(events, c(2, 3, 1))
  }
  print("Mean Duration...")
  duration <- apply(events, c(1, 2), mean, na.rm = T)
  print("Number of Events...")
  nevents <- apply(events, c(1, 2), function(x) length(x[!is.na(x)]))

  out <- list(track = newfield, percentage = meanfield, duration = duration, events = events, nevents = nevents)
  print(quantile(meanfield))
  print(min(duration, na.rm = T))
  return(out)
}

# large scale extension with further implementation
largescale.extension.if <- function(ics, ipsilon, field) {
  print("Large Scale Extension based on fixed angle")
  fimin <- 30 # southern latitude to be analyzed
  fimax <- 75 # northern latitude to be analyzed
  yreso <- ipsilon[2] - ipsilon[1]
  xreso <- ics[2] - ics[1]
  passo <- 5 / xreso # horizontal movemenent
  vertical <- 2.5 / yreso # vertical movement
  # time=1:length(field[1,1,]) #elements of the length of the dataset
  time <- which(apply(field, 3, max) != 0) # elements length of the dataset (removing no blocked days)

  print(paste("Box dimension:", passo * 2 * xreso, "° lon x ", vertical * 2 * yreso, "° lat"))

  short <- function(ics, ipsilon, field, passo, vertical) {
    control <- field
    range <- which.min(abs(ipsilon - fimin)):which.min(abs(ipsilon - fimax)) # check range for latitude excursion
    # range=range[(1+vertical):(length(range)-vertical)] #reduce range considering border effect
    new <- rbind(field, field, field) # bind domain for cross-date line
    for (i in 1:length(ics)) {
      ii <- i + length(ics)
      if (!all(new[(ii - passo):(ii + passo), ] == 0)) { # check to speed up
        for (j in range) {
          control[i, j] <- mean(new[(ii - passo):(ii + passo), (j - vertical):(j + vertical)], na.rm = T)
        }
      }
    }
    control[control > 0] <- 1
    return(control)
  }

  tt <- length(time)
  for (t in time) {
    progression.bar(t, tt)
    field[, , t] <- short(ics, ipsilon, field[, , t], passo, vertical)
  }
  return(field)
}


# Longitude filter for minimum extension
longitude.filter <- function(ics, ipsilon, field) {
  print("Longitude filter based on fixed angle")
  out <- field
  yreso <- ipsilon[2] - ipsilon[1]
  xreso <- ics[2] - ics[1]
  startipsilon <- which.min(abs(ipsilon - 30))
  estension <- (75 - 30) / yreso
  passo <- 15 / xreso

  print(paste("Continous longitude contrain", passo * xreso, "° lon"))

  tt <- length(field[1, 1, ])
  for (t in 1:tt) {
    progression.bar(t, tt)

    new <- rbind(field[, , t], field[, , t], field[, , t])
    for (j in startipsilon:((startipsilon + estension))) {
      new[, j] <- time.persistence(new[, j], persistence = passo)
    }
    field[, , t] <- new[length(ics) + (1:length(ics)), ]
  }
  return(field)
}
