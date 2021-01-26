##########################################################
#-------------------Time Avg functions-------------------#
##########################################################

# generalized function for time averaging based on conditon, using preallocation, vectorization and rowMeans
# use power.date.new or PCICt object to define the condition
time.mean <- function(ics, ipsilon, field, condition) {
  tmean <- array(NA, dim = c(length(ics), length(ipsilon), length(unique(condition))))
  for (t in unique(condition)) {
    tmean[, , which(t == unique(condition))] <- rowMeans(field[, , which(t == condition), drop = F], dims = 2)
  }
  return(tmean)
}


# fast function for monthly mean, using preallocation, vectorization and rowMeans
monthly.mean <- function(ics, ipsilon, field, etime) {
  condition <- paste(etime$month, etime$year)
  monthly <- array(NA, dim = c(length(ics), length(ipsilon), length(unique(condition))))
  for (t in unique(condition)) {
    monthly[, , which(t == unique(condition))] <- rowMeans(field[, , t == condition, drop = F], dims = 2)
  }
  return(monthly)
}

# introduce running mean, options for loop dataset
run.mean <- function(field, n = 5, loop = F) {
  nn <- floor(n / 2)
  if (loop) {
    field <- c(tail(field, n), field, head(field, n))
  }

  runfield <- field * NA
  for (t in (1 + nn):(length(field) - nn)) {
    if (!is.na(field[t])) {
      runfield[t] <- mean(field[(t - nn):(t + nn)], na.rm = T)
    } else {
      runfield[t] <- NA
    }
  }

  if (loop) {
    runfield <- runfield[(1 + n):(length(runfield) - n)]
  }

  return(runfield)
}

# improve running mean
# use vectorization for a 5 day running mean ad-hoc function (to be generalized!)
# about 10 times faster that a standard running mean function based on for loop
run.mean5 <- function(field) {
  newfield <- rowMeans(cbind(
    c(field[3:length(field)], NA, NA),
    c(field[2:length(field)], NA), field,
    c(NA, field[1:(length(field) - 1)]),
    c(NA, NA, field[1:(length(field) - 2)])
  ),
  na.rm = T
  )
  return(newfield)
}


# this a generalization of run.mean5() using vectorization for a faster result
# result is more accurate than run.mean() since on the edges it uses all possible points
run.mean.fast <- function(field, n = 5) {

  # check for even numbers
  if (n %% 2 == 0) {
    warning("Even number, replacing with its smaller odd one")
    n <- n - 1
  }

  # prepare the loop and create n-dimensianal matrix of lagged vectors
  nn <- floor(n / 2)
  newfield <- NULL
  for (k in -nn:nn) {
    if (k < 0) {
      newfield <- cbind(newfield, c(rep(NA, abs(k)), field[1:(length(field) + k)]))
    } else {
      newfield <- cbind(newfield, c(field[(1 + k):length(field)], rep(NA, k)))
    }
  }

  # apply rowMeans to produce the final mean
  finalfield <- rowMeans(newfield, na.rm = T)
  return(finalfield)
}



# function for daily anomalies, use array predeclaration and rowMeans (40 times faster!)
daily.anom.mean <- function(ics, ipsilon, field, etime) {
  condition <- paste(etime$day, etime$month)
  daily <- array(NA, dim = c(length(ics), length(ipsilon), length(unique(condition))))
  anom <- field * NA
  for (t in unique(condition)) {
    if (sum(t == condition) == 1) {
      print("Cannot compute a mean with a single value: using climatological mean")
      # anom <- sweep(field, 1:2, apply(field, 1:2, mean), "-")
      anom[, , which(t == condition)] <- rowMeans(field, dims = 2)
    } else {
      daily[, , which(t == unique(condition))] <- rowMeans(field[, , t == condition], dims = 2)
      anom[, , which(t == condition)] <- sweep(field[, , which(t == condition)], c(1, 2), daily[, , which(t == unique(condition))], "-")
    }
  }
  return(anom)
}


# function to bias correct field1 replacing field1 climatology with field2 climatology
# this can be useful to improve blocking performance, it is used when bcblock (bias-correct-block) is used
daily.replace.mean <- function(ics, ipsilon, field1, field2, etime) {
  condition <- paste(etime$day, etime$month)
  daily <- array(NA, dim = c(length(ics), length(ipsilon), length(unique(condition))))
  anom <- field1 * NA
  for (t in unique(condition)) {
    if (sum(t == condition) == 1) {
      print("Cannot compute a mean with a single value: using climatological mean")
      # anom <- sweep(field, 1:2, apply(field, 1:2, mean), "-")
      stop()
    } else {
      daily[, , which(t == unique(condition))] <- (-rowMeans(field1[, , t == condition], dims = 2) + rowMeans(field2[, , t == condition], dims = 2))
      anom[, , which(t == condition)] <- sweep(field1[, , which(t == condition)], c(1, 2), daily[, , which(t == unique(condition))], "+")
    }
  }
  return(anom)
}


# beta function for daily anomalies plus running mean (only 50% slower that standard daily avg)
daily.anom.run.mean5.old <- function(ics, ipsilon, field, etime) {
  condition <- paste(etime$day, etime$month)
  daily <- array(NA, dim = c(length(ics), length(ipsilon), length(unique(condition))))
  for (t in unique(condition)) {
    if (sum(t == condition) == 1) {
      stop("Cannot compute a mean with a single value")
    }
    daily[, , which(t == unique(condition))] <- rowMeans(field[, , t == condition], dims = 2)
  }
  rundaily <- apply(daily, c(1, 2), run.mean5)
  anom <- field * NA
  for (t in unique(condition)) {
    anom[, , which(t == condition)] <- sweep(field[, , which(t == condition)], c(1, 2), daily[, , which(t == unique(condition))], "-")
  }
  return(anom)
}

# function for daily anomalies from a 5-day running mean seasonal cycle:
# This that takes into account correctly the cross-year season as DJF
daily.anom.run.mean5 <- function(ics, ipsilon, field, etime) {

  # define condition for time selection: use numeric-compatible format
  condition <- format(etime$data, "%m%d")

  # evalute the time-ordered condition: if there is a jump, it means that there is a cross-year
  sorted <- sort(unique(condition))
  breakpoint <- which(diff(as.numeric(sorted)) > 100)
  # print(breakpoint)

  # if there is a cross-year, re-arrenge in order to have them as consecutive dates
  if (length(breakpoint) > 0) {
    sorted <- c(sorted[(breakpoint + 1):length(sorted)], sorted[1:breakpoint])
  }

  # print(sorted)
  # compute the seasonal cycle
  daily <- array(NA, dim = c(length(ics), length(ipsilon), length(sorted)))
  for (t in sorted) {
    if (sum(t == condition) == 1) {
      warning("Cannot compute a mean with a single value: using the value as it is!")
      daily[, , which(t == sorted)] <- field[, , t == condition]
    } else {
      daily[, , which(t == sorted)] <- rowMeans(field[, , t == condition], dims = 2)
    }
  }

  # apply running mean on the rightly orderded seasonal cycle
  rundaily <- apply(daily, c(1, 2), run.mean5)

  # remove seasonal cycle
  anom <- field * NA
  for (t in sorted) {
    anom[, , which(t == condition)] <- sweep(field[, , which(t == condition)], c(1, 2), daily[, , which(t == sorted)], "-")
  }

  return(anom)
}
