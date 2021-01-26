##########################################################
#------------EOFs and regims functions-------------------#
##########################################################

eofs <- function(lon, lat, field, neof = 4, xlim = NULL, ylim = NULL, weight = T,
                 method = "SVD", do_standardize = F, do_regression = F, verbose = T, true_coeff = F) {
  # R tool for computing EOFs based on Singular Value Decomposition ("SVD", default)
  # or with the eigenvectors of the covariance matrix ("covariance", slower)
  # If requested, computes linear regressions and standardizes the PCs
  # If you want to use the regressions, remember to standardize the PCs
  # Take as input a 3D anomaly field.
  # Requires "personal" functions area.weight, whicher and standardize

  # area weighting, based on the root of cosine
  if (weight) {
    printv("Area Weighting...", verbose)
    ww <- area.weight(lon, lat, root = T)
    wwfield <- sweep(field, c(1, 2), ww, "*")
  } else {
    wwfield <- field
  }

  # selection of the xbox and ybox if defined
  if (!is.null(xlim)) {
    lonselect <- whicher(lon, xlim[1]):whicher(lon, xlim[2])
  } else {
    lonselect <- 1:length(lon)
  }

  if (!is.null(ylim)) {
    latselect <- whicher(lat, ylim[1]):whicher(lat, ylim[2])
  } else {
    latselect <- 1:length(lat)
  }

  # box
  box <- wwfield[lonselect, latselect, , drop = F]
  slon <- lon[lonselect]
  slat <- lat[latselect]

  # transform 3D field in a matrix
  new_box <- array(box, dim = c(dim(box)[1] * dim(box)[2], dim(box)[3]))

  # calling SVD
  if (method == "SVD") {
    printv("Calling SVD...", verbose)
    SVD <- svd(new_box, nu = neof, nv = neof)

    # extracting EOFs (loading pattern), expansions coefficient and variance explained
    pattern <- array(SVD$u, dim = c(dim(box)[1], dim(box)[2], neof))
    coefficient <- SVD$v
    true_variance <- SVD$d[1:neof]^2
    variance <- (SVD$d[1:neof])^2 / sum((SVD$d)^2)
    if (do_standardize) {
      coefficient <- apply(coefficient, c(2), standardize)
    } else {
      coefficient <- sweep(coefficient, c(2), sqrt(variance), "*")
      # if (true_coeff) {
      #  coefficient <- sweep(coefficient, c(2), SVD$d[1:neof], "*")
      # }
    }
  }

  # calling covariance matrix
  if (method == "covariance") {
    printv("Calling eigenvectors of the covariance matrix...", verbose)
    covma <- cov(t(new_box))
    eig <- eigen(covma)
    coef <- (t(new_box) %*% eig$vector)[, 1:neof]
    pattern <- array(eig$vectors, dim = c(dim(box)[1], dim(box)[2], dim(box)[3]))[, , 1:neof]
    variance <- eig$values[1:neof] / sum(eig$values)
    if (do_standardize) {
      coefficient <- apply(coef, c(2), standardize)
    } else {
      coefficient <- coef
    }
  }

  # linear regressions on anomalies
  regression <- NULL
  if (do_regression) {
    printv("Linear Regressions (it can takes a while)... ", verbose)
    regression <- array(NA, dim = c(length(lon), length(lat), neof))
    # for (i in 1:neof) {regression[,,i]=apply(field,c(1,2),function(x) coef(lm(x ~ coefficient[,i]))[2])}
    for (i in 1:neof) {
      regression[, , i] <- apply(field, c(1, 2), function(x) lin.fit(as.matrix(coefficient[, i], ncol = 1), x)$coefficients)
    }
  }

  # preparing output
  printv("Finalize...", verbose)
  pattern <- list(x = slon, y = slat, z = pattern)
  out <- list(pattern = pattern, coeff = coefficient, variance = variance, regression = regression, true_variance = true_variance)
  return(out)
}

eofs.coeff <- function(lon, lat, field, eof_object, do_standardize = F, verbose = F) {
  # Computes expansion coefficient (i.e. PCs) of a given dataset on the
  # loading pattern of EOF previously computed
  # Works only on eof_object obtained with "eofs" function

  # Area weighting, based on the root of cosine
  printv("Area Weighting...", verbose)
  ww <- area.weight(lon, lat, root = T)
  wwfield <- sweep(field, c(1, 2), ww, "*")

  # selection of the box
  xlim <- c(min(eof_object$pattern$x), max(eof_object$pattern$x))
  ylim <- c(min(eof_object$pattern$y), max(eof_object$pattern$y))
  box <- wwfield[whicher(lon, xlim[1]):whicher(lon, xlim[2]), whicher(lat, ylim[1]):whicher(lat, ylim[2]), ]

  # transform 3D field in a matrix
  new_box <- array(box, dim = c(dim(box)[1] * dim(box)[2], dim(box)[3]))
  new_pattern <- array(eof_object$pattern$z, dim = c(dim(eof_object$pattern$z)[1] * dim(eof_object$pattern$z)[2], dim(eof_object$pattern$z)[3]))

  # projects the coefficients
  coef <- (t(new_box) %*% new_pattern)

  # standardize
  if (do_standardize) {
    coefficient <- apply(coef, c(2), standardize)
  } else {
    coefficient <- coef
  }

  print("Finalize...")
  return(coefficient)
}


regimes <- function(lon, lat, field, ncluster = 4, ntime = 1000, neof = 10, xlim, ylim, alg = "Hartigan-Wong") {
  # R tool to compute cluster analysis based on k-means.
  # Requires "personal" function eofs
  # Take as input a 3D anomaly field

  # Reduce the phase space with EOFs: use SVD and do not standardize PCs
  print("Launching EOFs...")
  t0 <- proc.time()
  reducedspace <- eofs(lon, lat, field, neof = neof, xlim = xlim, ylim = ylim, method = "SVD", do_regression = F, do_standardize = F)
  t1 <- proc.time() - t0
  # print(t1)

  # extract the principal components
  PC <- reducedspace$coeff
  print(str(PC))

  # k-means computation repeat for ntime to find best solution.
  print("Computing k-means...")
  t0 <- proc.time()
  print(str(ncluster))
  regimes <- kmeans(PC, as.numeric(ncluster), nstart = ntime, iter.max = 1000, algorithm = alg)
  t1 <- proc.time() - t0
  # print(t1)

  # Extract regimes frequencyr and timeseries of occupation
  cluster <- regimes$cluster
  frequencies <- regimes$size / dim(field)[3] * 100
  print(frequencies[order(frequencies, decreasing = T)])
  # print(regimes$tot.withinss)

  print("Creating Composites...")
  compose <- aperm(apply(field, c(1, 2), by, cluster, mean), c(2, 3, 1))

  # sorting from the more frequent to the less frequent
  kk <- order(frequencies, decreasing = T)
  cluster <- cluster + 10
  for (ss in 1:ncluster) {
    cluster[cluster == (ss + 10)] <- which(kk == ss)
  }

  # prepare output
  print("Finalize...")
  out <- list(cluster = cluster, frequencies = frequencies[kk], regimes = compose[, , kk], tot.withinss = regimes$tot.withinss)
  return(out)
}

regimes2 <- function(lon, lat, field, ncluster = 4, ntime = 1000, minvar = 0.8,
                     xlim, ylim, alg = "Hartigan-Wong") {

  # R tool to compute cluster analysis based on k-means.
  # Requires "personal" function eofs (see above)
  # Take as input a 3D anomaly field

  # Reduce the phase space with EOFs: use SVD and do not standardize PCs
  print("Launching EOFs...")
  t0 <- proc.time()
  reducedspace <- eofs(lon, lat, field, neof = 20, xlim = xlim, ylim = ylim, method = "SVD", do_regression = F, do_standardize = F)
  t1 <- proc.time() - t0
  # print(t1)
  reqPC <- which(cumsum(reducedspace$variance) > minvar)[1]
  print(paste("Retaining", reqPC, "PCs to fullfil minimum explained variance required (", minvar * 100, "%)"))

  # extract the principal components
  PC <- reducedspace$coeff[, 1:reqPC]
  print(str(PC))

  # k-means computation repeat for ntime to find best solution.
  print("Computing k-means...")
  t0 <- proc.time()
  print(str(ncluster))
  regimes <- kmeans(PC, as.numeric(ncluster), nstart = ntime, iter.max = 100, algorithm = alg)
  t1 <- proc.time() - t0
  # print(t1)

  # Extract regimes frequencyr and timeseries of occupation
  cluster <- regimes$cluster
  frequencies <- regimes$size / dim(field)[3] * 100
  print(frequencies[order(frequencies, decreasing = T)])
  # print(regimes$tot.withinss)

  print("Creating Composites...")
  compose <- aperm(apply(field, c(1, 2), by, cluster, mean), c(2, 3, 1))

  # sorting from the more frequent to the less frequent
  kk <- order(frequencies, decreasing = T)
  cluster <- cluster + 10
  for (ss in 1:ncluster) {
    cluster[cluster == (ss + 10)] <- which(kk == ss)
  }

  # prepare output
  print("Finalize...")
  out <- list(cluster = cluster, frequencies = frequencies[kk], regimes = compose[, , kk], tot.withinss = regimes$tot.withinss)
  return(out)
}
