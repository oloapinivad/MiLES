##########################################################
#-------Dynamics and thermodynamics functions------------#
##########################################################


# function to compute the beta parameter
beta.param <- function(longitude, latitude) {
  beta1d <- 2 * omega * cos(rad(latitude)) / Earth.Radius
  beta2d <- array(rep(beta1d, each = length(longitude)),
    dim = c(length(longitude), length(latitude))
  )
  return(list(beta1d = beta1d, beta2d = beta2d))
}

# function to compute the f parameter
coriolis.param <- function(longitude, latitude) {
  f1d <- 2 * omega * sin(rad(latitude))
  f2d <- array(rep(f1d, each = length(longitude)),
    dim = c(length(longitude), length(latitude))
  )
  return(list(f1d = f1d, f2d = f2d))
}

# function to compute the stationary wave number aka the refraction index
# based on Willis et al 2019, spherical coordinates
# verified (but r*cos(lat) has to be added, unclear)
refraction.index.willis <- function(lon, lat, ua) {
  der_ua <- der.2nd.spherical(rad(lon), rad(lat), ua, along = "lat", wind = T)
  der2_ua <- der.2nd.spherical(rad(lon), rad(lat), der_ua, along = "lat")
  beta <- beta.param(lon, lat)$beta2d
  refraction <- sqrt((beta - der2_ua) / ua) * Earth.Radius * cos(grid.extend(rad(lon), rad(lat))$y)
  return(refraction)
}


# function to compute the stationary wave number aka the refraction index
# based on the Mercator coordinates by Hoskins and Ambrizzi 1993 (verified)
refraction.index <- function(lon, lat, ua) {
  grid <- grid.extend(rad(lon), rad(lat))
  wind <- ua / (cos(grid$y) * Earth.Radius)
  der_wind <- der.2nd(wind * cos(grid$y)^2, rad(lon), rad(lat), along = "y") / cos(grid$y)
  der_wind2 <- der.2nd(der_wind, rad(lon), rad(lat), along = "y") / cos(grid$y)
  beta_m <- (2 * omega - der_wind2) * cos(grid$y)^2 / Earth.Radius
  refraction <- sqrt(Earth.Radius * beta_m / wind)
  refraction[ua < 0] <- NA
  return(list(ks = refraction, ks_int = round(refraction), bm = beta_m))
}

# relative vorticity 
# computed the derivative in spherical coordinates
# verified with ERA5 output
relative.vorticity <- function(lon, lat, ua, va) {
  dvdx <- der.2nd.spherical(rad(lon), rad(lat), va, along = "lon", wind = T)
  dudy <- der.2nd.spherical(rad(lon), rad(lat), ua, along = "lat", wind = T)
  rv <- dvdx - dudy
  return(rv)
}

# planetary vorticity, taken from Coriolis parameter
planetary.vorticity <- function(lon, lat) {
  return(coriolis.param(lon, lat)$f2d)
}

# divergent wind from velocity potential
divergent.wind <- function(lon, lat, velopot) {
  x <- der.2nd.spherical(rad(lon), rad(lat), velopot, along = "lon")
  y <- der.2nd.spherical(rad(lon), rad(lat), velopot, along = "lat")
  div <- list(x = x, y = y)
  return(div)
}

# rossby wave source (from Scaife et al.) 
# to be verified
rws <- function(lon, lat, div_wind, absvorticity) {
  ft <- absvorticity * (
    der.2nd.spherical(rad(lon), rad(lat), div_wind$x, along = "lon", wind = T) +
      der.2nd.spherical(rad(lon), rad(lat), div_wind$y, along = "lat", wind = T)
  )
  st <- div_wind$x * der.2nd.spherical(rad(lon), rad(lat), absvorticity, along = "lon") +
    div_wind$y * der.2nd.spherical(rad(lon), rad(lat), absvorticity, along = "lat")
  r <- -(ft + st)
  return(r)
}

# Deformation vector used in Davini et 2017, from Cai and Mak 1990
deformation <- function(lon, lat, ua, va) {
  x <- der.2nd.spherical(rad(lon), rad(lat), ua, along = "lon", wind = T) -
    der.2nd.spherical(rad(lon), rad(lat), va, along = "lat", wind = T)
  y <- der.2nd.spherical(rad(lon), rad(lat), va, along = "lon", wind = T) +
    der.2nd.spherical(rad(lon), rad(lat), ua, along = "lat", wind = T)
  d <- list(x = x, y = y)
  return(d)
}

# barotropic energy conversion: scalar product between modified hoskins
# E-vector (has to be computed before) and deformation vector
# verified
barotropic.energy.conversion <- function(hoskins, deform) {
  bec <- hoskins$x * deform$x + hoskins$y * deform$y
  return(bec)
}

# Momentum convergence as a divergence of the E-vector
# factor of two induced by the estimation of Ex used in my CDO code
# to be verified
momentum.convergence <- function(lon, lat, ex, ey) {
  m <- der.2nd.spherical(rad(lon), rad(lat), ex / 2, along = "lon", wind = T) +
    der.2nd.spherical(rad(lon), rad(lat), ey, along = "lat", wind = T)
  return(m)
}

# Exner function
exner <- function(level) {
  return((p0 / level)^(Rd / cp))
}

# Potential Temperature
potential.temperature <- function(level, ta) {
  pt <- ta * exner(level)
  return(pt)
}

# Brunt Vaisal frequency
brunt.vaisala <- function(lon, lat, level, ta, level_down, ta_down, level_up, ta_up) {
  dthetadz <- der.2nd.dz(
    level_down, level_up, potential.temperature(level_down, ta_down),
    potential.temperature(level_up, ta_up)
  )
  n <- sqrt(g0 / potential.temperature(level, ta) * dthetadz)
  return(n)
}

# Eady Growth Rate
eady.growth <- function(lon, lat, level_down, ua_down, level_up, ua_up, vaisala) {
  dudz <- der.2nd.dz(level_down, level_up, ua_down, ua_up)
  sigma <- 0.3068 * coriolis.param(lon, lat)$f2d * abs(dudz) / vaisala
  return(sigma)
}

# Static stability
static.stability <- function(level, dthetadp) {
  s <- -dthetadp / (Rd / p0 * (p0 / level)^(cv / cp))
  return(s)
}

# Baroclinic energy conversion
baroclinic.energy.conversion <- function(lon, lat, static, heatflux, theta) {
  bec <- -1 / static * heatflux * der.2nd.spherical(rad(lon), rad(lat), theta, along = "lat")
  return(bec)
}

# Potential vorticity on isobaric level
# at first order only the component on z axis is enough
# however differences with ERA5 PV arise (perhaps due to computation on model levels?)
potential.vorticity <- function(lon, lat, ua, va, dthetadp) {
  absolute_vorticity <- relative.vorticity(lon, lat, ua, va) + planetary.vorticity(lon, lat)
  tt <- absolute_vorticity * dthetadp
  PV <- -g0 * tt
  return(PV)
}


# potential vorticity on pressure level (complete)
# a little bit more precise, but still differences found with ERA5
potential.vorticity.complete <- function(lon, lat, ua, va, theta, duadp, dvadp, dthetadp) {
  absolute_vorticity <- relative.vorticity(lon, lat, ua, va) + planetary.vorticity(lon, lat)
  ft <- dvadp * der.2nd.spherical(rad(lon), rad(lat), theta, along = "lon")
  st <- duadp * der.2nd.spherical(rad(lon), rad(lat), theta, along = "lat")
  tt <- absolute_vorticity * dthetadp
  PV <- -g0 * (-ft + st + tt)
  return(PV)
}
