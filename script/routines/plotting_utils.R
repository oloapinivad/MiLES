# plotting utils
lettering <- paste0("(", letters, ")")
cex.letter <- 2
map_projection <- "azequalarea"
lat_lim <- c(30, 90)
plotpar <- list(
  cex.main = 2.5, cex.axis = 1.5, cex.lab = 1.5,
  mar = c(5, 2, 5, 9), oma = c(1, 1, 3, 3)
)

# Rcolorbrewer palette
library("RColorBrewer")
palette.rdbl <- colorRampPalette(rev(brewer.pal(11, "RdBu")))
palette.bupu <- colorRampPalette((brewer.pal(9, "BuPu")))
palette.puor <- colorRampPalette(rev(brewer.pal(11, "PuOr")))
palette.spct <- colorRampPalette(rev(brewer.pal(9, "Spectral")))

# imagescale3 color bar details
imgscl_colorbar <- 1.4
imgscl_label <- 1.5
imgscl_line <- 3

# dark colors
darken <- function(color, factor = 1.5) {
  col <- col2rgb(color)
  col <- col / factor
  col <- rgb(t(col), maxColorValue = 255)
  col
}

# function to creata a rectangle of lon/lat
create.box <- function(lons, lats, map_projection = "no") {
  rect <- c(NA, NA)
  rect <- rbind(rect, cbind(rep(lons[1], length(lats[1]:lats[2])), lats[1]:lats[2]))
  rect <- rbind(rect, cbind(lons[1]:lons[2], rep(lats[2], length(lons[1]:lons[2]))))
  rect <- rbind(rect, cbind(rep(lons[2], length(lats[1]:lats[2])), lats[2]:lats[1]))
  rect <- rbind(rect, cbind(lons[2]:lons[1], rep(lats[1], length(lons[2]:lons[1]))))
  rect <- rect[2:length(rect[, 1]), ]
  if (map_projection == "no") {
    rr <- list(x = rect[, 1], y = rect[, 2])
  } else {
    rr <- mapproject(rect[, 1], rect[, 2], projection = map_projection, orientation = c(90, 0, 0))
  }
  return(rr)
}
