# spatial operations, weights, averages and correlations...

# function to compute area cell (i.e. weights)
# starting from lon-lat bounds
area_weights <- function(lon, lat) {
  # bounds (not tested for regular grid)
  lon_len <- diff(na.omit(stats::filter(c(last(lon) - 360, lon, lon[1] + 360), c(0.5, 0.5), sides = 2))) *
    2 * pi * Earth.Radius / 360
  lat_len <- diff(c(-90, na.omit(stats::filter(lat, c(0.5, 0.5), sides = 2)), 90)) *
    2 * pi * Earth.Radius / 360

  # 2d grid
  grid <- grid.extend(lon_len, lat_len)

  # area
  area <- sweep(grid$x, 2, cos(rad(lat)), "*") * grid$y

  return(area)
}

# produce a 2d matrix of area weight
# WARNING, REDUNDANCY
area.weight <- function(ics, ipsilon, root = T) {
  field <- array(NA, dim = c(length(ics), length(ipsilon)))
  if (root == T) {
    for (j in 1:length(ipsilon)) {
      field[, j] <- sqrt(cos(pi / 180 * ipsilon[j]))
    }
  }

  if (root == F) {
    for (j in 1:length(ipsilon)) {
      field[, j] <- cos(pi / 180 * ipsilon[j])
    }
  }

  return(field)
}


# function to compute pressure level weights based on
# pressure level bounds - used for zonal mean
zonal_weights <- function(lat, plev) {
  lat_len <- diff(c(-90, na.omit(filter(lat, c(0.5, 0.5), sides = 2)), 90)) *
    2 * pi * Earth.Radius / 360
  plev_depth <- -diff(c(100000, na.omit(filter(c(plev), c(0.5, 0.5))), 0))

  # 2d grid
  grid <- grid.extend(lat_len, plev_depth)

  # area (do not take it as real area, useful only for weights)
  area <- grid$x * grid$y

  return(area)
}

# weighted sum (normalize weights)
weighted.sum <- function(x, w) {
  w <- w / sum(w)
  return(sum(x * w))
}

# weighted correlation
weighted.cor <- function(x, y, w) {
  w.mean.x <- sum(w * x) / sum(w)
  w.mean.y <- sum(w * y) / sum(w)

  w.cov.xy <- sum(w * (x - w.mean.x) * (y - w.mean.y)) / sum(w)
  w.var.y <- sum(w * (y - w.mean.y) * (y - w.mean.y)) / sum(w)
  w.var.x <- sum(w * (x - w.mean.x) * (x - w.mean.x)) / sum(w)

  corr <- w.cov.xy / sqrt(w.var.x * w.var.y)
  return(corr)
}

# weighted standard deviations
weighted.sd <- function(x, w) {
  w.mean <- sum(w * x) / sum(w)
  v1 <- sum(w)
  v2 <- sum(w^2)
  var <- v1 / (v1^2 - v2) * sum(w * (x - w.mean)^2)
  sdd <- sqrt(var)
  return(sdd)
}


# function to expand the grid for latitude and longitude
grid.extend <- function(longitude, latitude) {
  latexp <- array(rep(latitude, each = length(longitude)),
    dim = c(length(longitude), length(latitude))
  )
  lonexp <- array(rep(longitude, length(latitude)),
    dim = c(length(longitude), length(latitude))
  )
  return(list(x = lonexp, y = latexp))
}
