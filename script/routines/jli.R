jli <- function(lon, lat, field, jetstyle = "atlantic", bw = 0.5,
                npoints = 128, lon_sel = NULL, lat_sel = NULL) {

  # bandwitdht
  band <- diff(lon)[1] * bw

  # select a type of JLI
  if (jetstyle == "atlantic") {
    lon_bounds <- c(-60, 0)
    lat_bounds <- c(20, 70)
  } else if (jetstyle == "global") {
    lon_bounds <- c(-180, 180)
    lat_bounds <- c(0, 70)
  } else {
    lon_bounds <- lon_sel
    lat_bounds <- lat_sel
  }

  # boundary selection
  ilon_bounds <- whicher(lon, lon_bounds[1]):whicher(lon, lon_bounds[2])
  ilat_bounds <- whicher(lat, lat_bounds[1]):whicher(lat, lat_bounds[2])
  short_lat <- lat[ilat_bounds]

  # subselection and cleaning from NA

  subfield <- field[ilon_bounds, ilat_bounds, which(!is.na(field[1, 1, ]))]
  print(paste("Removing", length(which(is.na(field[1, 1, ]))), "NA records"))
  print(str(subfield))

  # 1D JLI
  zon_mean <- colMeans(subfield, dims = 1)
  jli_lat <- short_lat[apply(zon_mean, 2, which.max)]
  jli_speed <- apply(zon_mean, 2, max)

  # 2D JLI
  jli2d_speed <- apply(subfield, c(1, 3), max)
  jli2d_lat <- apply(subfield, c(1, 3), function(x) short_lat[which.max(x)])

  # densities
  density_2d <- t(apply(jli2d_lat, 1, function(x) {
    density(x, n = npoints, bw = band, from = lat_bounds[1], to = lat_bounds[2])$y
  }))
  density_1d <- density(jli_lat, n = npoints, bw = band, from = lat_bounds[1], to = lat_bounds[2])



  out <- list(
    lon = lon[ilon_bounds], lat = short_lat, jli_lat = jli_lat,
    jli_speed = jli_speed, jli_lat_2d = jli2d_lat, jli_speed_2d = jli2d_speed,
    density_1d = density_1d, density_2d = density_2d, density_lat = density_1d$x,
    npoints = npoints
  )

  return(out)
}
