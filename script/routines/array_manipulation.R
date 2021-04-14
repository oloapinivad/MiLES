#########################################################
#-----------Array manipulations functions----------------#
##########################################################

# detect ics ipsilon lat-lon
whicher <- function(axis, number) {
  out <- which.min(abs(axis - number))
  return(out)
}

# last element of a vector
last <- function(x) {
  return(x[length(x)])
}

# Function to generalize through do.call() n-dimensional array subsetting
# and array indexing. Derived from Stack Overflow issue
# https://stackoverflow.com/questions/14500707/select-along-one-of-n-dimensions-in-array
array_indexing <- function(field, dim, value, drop = FALSE) {

  # Create list representing arguments supplied to [
  # bquote() creates an object corresponding to a missing argument
  indices <- rep(list(bquote()), length(dim(field)))
  indices[[dim]] <- value

  # do.call on the indices
  out <- do.call("[", c(list(field), indices, list(drop = drop)))

  return(out)
}

# n-dimensional generalized evolution of flipper()
# it reverts the ipsilon a required dimension: by default uses the second dimension
# used by ncdf.opener.universal()
flipper <- function(field, dim = 2) {
  ydim <- dim
  ll <- dim(field)[ydim]
  field <- array_indexing(field, ydim, ll:1)
  return(field)
}

# define rotate along longitudes based on array_indexing (Jan 2019)
# universal to every dimension, 50% faster than previous rotation() functon
# can be used also for ad hoc rotation
# used by ncdf.opener.universal()
rotation <- function(line, rotate, longitude.case = FALSE) {

  # dimension of the first dimension (the one to be rotated)
  if (is.null(dim(line))) {
    dims <- 1
    ll <- length(line)
  } else {
    ll <- dim(line)[1]
    dims <- length(dim(line))
  }

  # default options
  if (is.character(rotate)) {
    if (rotate == "full" | rotate == "reverse") {
      # 180 degrees rotation of longitude
      move1 <- 1 / 2 * ll
    } else if (rotate == "+half") {
      # 90 degree rotation (useful for TM90)
      move1 <- 1 / 4 * ll
    } else if (rotate == "-half") {
      # 270 degree rotation (useful for TM90)
      move1 <- 3 / 4 * ll
    } else if (rotate == "no") {
      return(line)
    }

    # numeric options
  } else {
    if (rotate == 0 | rotate >= ll) {
      print("Nothing to do")
      return(line)
    } else {
      move1 <- rotate + 1 # always add one to avoid crash
    }
  }

  # move2 as the difference of the number of points
  move2 <- ll - move1

  # special case for flipping longitudes
  if (longitude.case) {
    newline <- c(line[(move2 + 2):ll] - 360, line[1:(move2 + 1)])
    if (all(newline <= 0.5)) {
      newline <- newline + 360
    }
  } else {
    # create new elements order
    elements <- c((move2 + 2):ll, 1:(move2 + 1))

    # run the selection using array_indexing
    newline <- array_indexing(line, 1, elements)
  }

  return(newline)
}
