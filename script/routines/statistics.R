# t student
t.corr <- function(r, sample) {
  t <- r * (sqrt(sample - 2)) / sqrt(1 - r^2)
  xd <- seq(0.001, 0.999, 0.001)
  distr <- qt(xd, df = sample - 2)
  significance <- xd[which.min(abs(distr - t))]
  out <- list(corr = r, tvalue = t, signif = significance)
  return(out)
}

# extract significant for correlation
t.sign.corr <- function(significance, sample, method = "two-tailed", verbose = F) {
  if (method == "two-tailed") {
    sign <- significance + (1 - significance) / 2
  }
  if (method == "one-tailed") {
    sign <- significance
  }

  t <- qt(sign, df = sample - 2)
  r <- t / (sqrt(sample - 2 + t^2))
  if (verbose) {
    print(paste("Signifance level is", significance * 100, "% with", method, "distribution"))
    print(paste("Corresponding t-Student quantile is", t))
    print(paste("Critical Pearson's correlation coefficient is", r))
  }

  return(r)
}


# function to extract the first significant year (when all the following are significant)
year_mannkendall <- function(onepoint) {
  ff <- NA
  # remove empty values and set a startyear
  onepoint2 <- onepoint[!is.na(onepoint)]
  startyear <- 1990
  if (length(onepoint2) > 0 & !all(onepoint2 == 0)) {
    lastyear <- max(as.numeric(names(onepoint2)))
    # apply mannkendall
    # dput(onepoint2)
    test <- sapply(startyear:lastyear, function(x) {
      serie <- onepoint2[which(names(onepoint2) %in% (1951:x))]
      if (!all(serie == 0)) {
        return(MannKendall(serie)$sl)
      } else {
        return(NULL)
      }
    })
    names(test) <- startyear:lastyear
    check <- which(unlist(test) < 0.05)
    # print(check)
    # if there is a sequence, and the last value is from the last year, go for it
    if (length(check) > 0) {
      if (names(check[length(check)]) == as.numeric(names(onepoint2[length(onepoint2)]))) {
        if (length(check) == 1) {
          ff <- lastyear
        } else {
          # all the values after the first guess are significant
          if (rle(rev(diff(check)))$values[1] == 1) {
            lastsequence <- rle(rev(diff(check)))$lengths[1] - 1
            ff <- lastyear - lastsequence
          }
        }
      }
    }
  }
  return(ff)
}

# kendall
round.pi <- function(pivalue) {
  if (round(pivalue, 2) == 1) {
    value <- 0.99
  } else if (round(pivalue, 2) == 0.99 | round(pivalue, 2) == 0.98) {
    value <- round(pivalue, 2)
  } else {
    value <- 0.05 * round(round(pivalue, 2) / 0.05)
  }
  return(value)
}
