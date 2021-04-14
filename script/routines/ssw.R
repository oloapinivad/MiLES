#nothing more than the PDF of the U10 at 60N, it has to be expanded

ssw <- function(lon, lat, field) {

  zon_mean <- colMeans(field[,whicher(lat,60),])
  density_1d <- density(zon_mean, n = 128, bw = 5, from = -40, to = 80)
  out <- list(speed = zon_mean, density_1d = density_1d)
  return(out)
} 


