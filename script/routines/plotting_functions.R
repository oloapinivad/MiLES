##########################################################
#--------------Plotting functions------------------------#
##########################################################

# function to open devices
open.plot.device <- function(figname, output_file_type, CFGSCRIPT, special = FALSE) {
  # Chose output format for figure - by JvH
  source(CFGSCRIPT)
  if (special == FALSE) {
    if (tolower(output_file_type) == "png") {
      png(filename = figname, width = png_width, height = png_height)
    } else if (tolower(output_file_type) == "pdf") {
      pdf(file = figname, width = pdf_width, height = pdf_height, onefile = T)
    } else if (tolower(output_file_type) == "eps") {
      setEPS(width = pdf_width, height = pdf_height, onefile = T, paper = "special")
      postscript(figname)
    }
  }

  # special case for TM90
  if (special == TRUE) {
    if (tolower(output_file_type) == "png") {
      png(filename = figname, width = png_width / af, height = png_height * af / 2)
    } else if (tolower(output_file_type) == "pdf") {
      pdf(file = figname, width = pdf_width / af, height = pdf_height * af / 2, onefile = T)
    } else if (tolower(output_file_type) == "eps") {
      setEPS(width = pdf_width / af, height = pdf_height * af / 2, onefile = T, paper = "special")
      postscript(figname)
    }
  }
}

# extensive filled.contour function
filled.contour3 <-
  function(x = seq(0, 1, length.out = nrow(z)),
           y = seq(0, 1, length.out = ncol(z)), z, xlim = range(x, finite = TRUE),
           ylim = range(y, finite = TRUE), zlim = range(z, finite = TRUE),
           levels = pretty(zlim, nlevels), nlevels = 20, color.palette = cm.colors,
           col = color.palette(length(levels) - 1), extend = FALSE, plot.title, plot.axes,
           key.title, key.axes, asp = NA, xaxs = "i", yaxs = "i", las = 1,
           image.scale = FALSE,
           axes = TRUE, frame.plot = axes, mar, ...) {
    # modification by Ian Taylor of the filled.contour function
    # to remove the key and facilitate overplotting with contour()
    # further modified by Carey McGilliard and Bridget Ferris
    # to allow multiple plots on one page
    # modification to allow plot outside boundaries

    if (missing(z)) {
      if (!missing(x)) {
        if (is.list(x)) {
          z <- x$z
          y <- x$y
          x <- x$x
        }
        else {
          z <- x
          x <- seq.int(0, 1, length.out = nrow(z))
        }
      }
      else {
        stop("no 'z' matrix specified")
      }
    }
    else if (is.list(x)) {
      y <- x$y
      x <- x$x
    }
    if (any(diff(x) <= 0) || any(diff(y) <= 0)) {
      stop("increasing 'x' and 'y' values expected")
    }

    if (extend) {
      z[z < min(levels)] <- min(levels)
      z[z > max(levels)] <- max(levels)
    }

    plot.new()
    plot.window(xlim, ylim, "", xaxs = xaxs, yaxs = yaxs, asp = asp)
    if (!is.matrix(z) || nrow(z) <= 1 || ncol(z) <= 1) {
      stop("no proper 'z' matrix specified")
    }
    if (!is.double(z)) {
      storage.mode(z) <- "double"
    }
    .filled.contour(as.double(x), as.double(y), z, as.double(levels),
      col = col
    )
    if (image.scale) {
      image.scale3(z, levels = as.double(levels), colorbar.label = "", color.palette = color.palette)
    }
    if (missing(plot.axes)) {
      if (axes) {
        title(main = "", xlab = "", ylab = "")
        Axis(x, side = 1, ...)
        Axis(y, side = 2, ...)
      }
    }
    else {
      plot.axes
    }
    if (frame.plot) {
      box()
    }
    if (missing(plot.title)) {
      title(...)
    } else {
      plot.title
    }
    invisible()
  }

image.scale3 <- function(z, levels, color.palette = heat.colors, colorbar.label = "image.scale", extend = T,
                         line.label = 3, line.colorbar = 1.5, cex.label = 1.5, cex.colorbar = 1.7, colorbar.width = 1.5, ...) {

  # save properties from main plotting region
  old.par <- par(no.readonly = TRUE)
  mfg.save <- par()$mfg
  old.fig <- par()$fig

  # defining plotting region with proper scaling
  # print(old.fig)
  xscal <- (old.fig[2] - old.fig[1])
  yscal <- (old.fig[4] - old.fig[3])
  lw <- colorbar.width
  lp <- line.colorbar / 100
  new.fig <- c(old.fig[2] - 0.07 * xscal * lw - lp, old.fig[2] - 0.03 * xscal - lp, old.fig[3] + 0.1 * yscal, old.fig[4] - 0.1 * yscal)

  # safety check
  new.fig[new.fig > 1] <- 1
  # print(old.fig)
  # print(new.fig)

  if (missing(levels)) {
    levels <- seq(min(z), max(z), , 12)
  }
  # fixing color palette
  col <- color.palette(length(levels) - 1)

  # starting plot
  par(mar = c(1, 1, 1, 1), fig = new.fig, new = TRUE)

  # creating polygons for legend
  poly <- vector(mode = "list", length(col))
  for (i in seq(poly)) {
    poly[[i]] <- c(levels[i], levels[i + 1], levels[i + 1], levels[i])
  }

  xlim <- c(0, 1)
  if (extend) {
    longer <- 1.5
    dl <- diff(levels)[1] * longer
    ylim <- c(min(levels) - dl, max(levels) + dl)
  } else {
    ylim <- range(levels)
  }
  plot(1, 1, t = "n", ylim = ylim, xlim = xlim, axes = FALSE, xlab = "", ylab = "", xaxs = "i", yaxs = "i", ...)
  for (i in seq(poly)) {
    polygon(c(0, 0, 1, 1), poly[[i]], col = col[i], border = NA)
  }

  if (extend) {
    polygon(c(0, 1, 1 / 2), c(levels[1], levels[1], levels[1] - dl),
      col = col[1], border = NA
    )
    polygon(c(0, 1, 1 / 2), c(levels[length(levels)], levels[length(levels)], levels[length(levels)] + dl),
      col = col[length(col)], border = NA
    )
    polygon(c(0, 0, 1 / 2, 1, 1, 1 / 2), c(
      levels[1], levels[length(levels)], levels[length(levels)] + dl, levels[length(levels)], levels[1],
      levels[1] - dl
    ), border = "black", lwd = 2)
    ylim0 <- range(levels)
    prettyspecial <- pretty(ylim0)
    prettyspecial <- prettyspecial[prettyspecial <= max(ylim0) & prettyspecial >= min(ylim0)]
    axis(4, las = 1, cex.axis = cex.colorbar, at = prettyspecial, labels = prettyspecial, ...)
  } else {
    box()
    axis(4, las = 1, cex.axis = cex.colorbar, ...)
  }

  # box, axis and leged
  mtext(colorbar.label, line = line.label, side = 4, cex = cex.label, ...)

  # resetting properties for starting a new plot (mfrow style)
  par(old.par)
  par(mfg = mfg.save, new = FALSE)
  invisible()
}

# function for interpolation and projection of a 2D field on a mapproj R projection
proj.plot <- function(lon, lat, field, lmin = NULL, proj = "azequalarea", param = NULL, orient = c(90, 0, 0), npoints = 201) {

  # default is azimuthal equal area map

  # required packages
  require(mapproj)
  require(akima)

  # it provides lower latitude limit for plots
  if (is.null(lmin)) {
    lmin <- min(lat)
  }

  # build grids
  lon.grid <- rep(lon, length(lat))
  lat.grid <- sort(rep(lat, length(lon)))

  # project grid
  proj.grid <- mapproject(lon.grid, lat.grid, projection = proj, parameters = param, orientation = orient)

  # provide limits for future plots (for polar projection)
  limiter <- mapproject(c(0, 90, 180, 270), rep(lmin, 4), proj = "", orientation = orient)
  xlims <- sort(c(limiter$x[2], limiter$x[4]))
  ylims <- sort(c(limiter$y[1], limiter$y[3]))

  # plot grid
  lon.plot <- seq(min(proj.grid$x, na.rm = T), max(proj.grid$x, na.rm = T), length.out = npoints)
  lat.plot <- seq(min(proj.grid$y, na.rm = T), max(proj.grid$y, na.rm = T), length.out = npoints)

  # interpolation (akima needed)
  good <- is.finite(field) & is.finite(proj.grid$x) & is.finite(proj.grid$y)
  projected <- interp(proj.grid$x[good], proj.grid$y[good], field[good], lon.plot, lat.plot, duplicate = "strip")
  return(projected = list(x = projected$x, y = projected$y, z = projected$z, xlim = xlims, ylim = ylims))
}

# addland function based on map which can handle projections
proj.addland <- function(lon, lat, proj = "no", orient = c(90, 0, 0), param = NULL, inter = F, color = "black") {

  # required packages
  require(maps)
  require(mapproj)

  if (proj == "no") {
    map("world", regions = ".", interior = inter, exact = F, boundary = T, add = T)
  } else {
    # get map, project and do the lines
    box()
    map("world", add = T, projection = proj, orientation = orient, parameter = param, interior = inter, exact = F, boundary = T)

    # default lines for northern hemisphere
    for (i in seq(-80, 80, 20)) {
      x0 <- lon
      y0 <- rep(i, length(lon))
      p <- mapproject(x0, y0, proj = "", orientation = orient)
      lines(p, lty = 3)
    }

    # default circles for northern hemisphere
    for (i in c(seq(-360, 360, 30))) {
      y0 <- seq(0, 90, , 90)
      x0 <- rep(i, 90)
      p <- mapproject(x0, y0, proj = "", orientation = orient)
      lines(p, lty = 3)
    }
  }
}


# rearrange arrays for use both standard plotting and proj.plot
plot.prepare <- function(ics, ipsilon, field, proj, lat_lim) {
  if (proj == "no") {
    outfile <- list(x = ics, y = ipsilon, z = field, xlim = range(ics), ylim = lat_lim, asp = NULL, xlab = "Longitude", ylab = "Latitude", axes = T)
  } else {
    field[is.na(field)] <- 0
    p <- proj.plot(ics, ipsilon, field, lmin = lat_lim[1], proj = proj, param = NULL, orient = c(90, 0, 0), npoints = 80)
    outfile <- list(x = p$x, y = p$y, z = p$z, xlim = p$xlim, ylim = p$ylim, xlab = "", ylab = "", axes = F, asp = 1)
  }
  return(outfile)
}
