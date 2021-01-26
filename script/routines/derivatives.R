# Second order centered derivative, under testing
# generalized function for any array from 1D up to 3D
# it works on cartesian coordinates (i.e. meters)
# to be extended for cases with lon-lat and for pressure
der.2nd <- function(field, x, y = NA, z = NA, along = "x") {
  fder <- function(ff, axis, lat) {
    out <- c(NA, diff(ff, 2) / diff(axis, 2), NA)
    return(out)
  }

  # check dimension of the field
  checkdim <- c(all(!is.na(x)), all(!is.na(y)), all(!is.na(z)))
  dims <- c("x", "y", "z")[checkdim]
  d <- length(dims)

  # 1d case
  if (d == 1) {
    outfield <- fder(field, get(along))
    # 2d-3d cases
  } else {

    # set the direction on which call apply
    applydir <- which(dims != along)

    # set array rotation (due to apply transposed output)
    rot <- rep(NA, d)
    rot[which(c("x", "y", "z") == along)] <- 1
    rot[which(is.na(rot))] <- c(2, 3)[1:d - 1]

    # apply derivative
    outfield <- apply(field, applydir, function(ff) {
      fder(ff, get(along))
    })

    # array permutation
    outfield <- aperm(outfield, rot)
  }

  return(outfield)
}

# function to compute the lon/lat derivative of a 2d field (2nd order)
# wind flag is used to activate/deactivate the derivation of the unity vector
# lon-lat should be in radiants
der.2nd.spherical <- function(lon, lat, field, along = "lat", wind = FALSE) {
  axis <- get(along)

  # safety check on lon/lat radiants
  if (max(lon) > 2 * pi) {
    print("Longitude is not in radiant, converting...")
    lon <- rad(lon)
  }
  if (max(lat) > pi) {
    print("Latitude is not in radiant, converting...")
    lat <- rad(lat)
  }

  # factor for derivation of unity vectors
  if (wind == TRUE) {
    cos_factor <- cos(lat)
  } else {
    cos_factor <- rep(1, length(lat))
  }

  out <- field * NA
  # print(paste("Derivative along", along))

  if (along == "lat") {
    incr <- Earth.Radius * diff(axis, 2) * cos_factor[2:(length(lat) - 1)]
    for (x in 1:length(lon)) {
      out[x, ] <- c(NA, diff(field[x, ] * cos_factor, 2) / incr, NA)
    }
  } else if (along == "lon") {
    for (y in 1:length(lat)) {
      incr <- Earth.Radius * diff(axis, 2) * cos(lat[y])
      loopincr <- c(incr[length(incr)], incr, incr[1])
      out[, y] <- diff(c(field[length(lon), y], field[, y], field[1, y]), 2) / loopincr
    }
  }
  return(out)
}


# define the divergence of a vector field
divergence <- function(field) {

}

# define the gradient of a scalar field
gradient <- function(field) {

}


# 2nd order derivative for lon/lat fields along the vertical direction
# need to provide the two levels
der.2nd.dz <- function(lev1, lev2, field1, field2) {
  dz <- (vertical(lev2) - vertical(lev1)) * 1000
  df <- (field2 - field1) / dz
  return(df)
}

# 2nd order derivative for lon/lat fields along the vertical direction
# need to provide the two levels
der.2nd.dp <- function(lev1, lev2, field1, field2) {
  df <- (field2 - field1) / (lev2 - lev1)
  return(df)
}
