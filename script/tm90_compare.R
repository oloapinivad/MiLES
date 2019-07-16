######################################################
#------Blocking extra project comparer MiLES----------#
#-------------P. Davini (Jul 19)---------------------#
######################################################

# DECLARING THE FUNCTION: EXECUTION IS AT THE BOTTOM OF THE SCRIPT
#miles.tm90.compare <- function(projects, dataset_obs, expid, ens, year1, year2,
#                                FIGDIR, FILESDIR, REFDIR, CFGSCRIPT) {

# folders that should come from outside
PROGDIR <- "/home/paolo/MiLES"
CFGSCRIPT <- "/home/paolo/MiLES/config/R_config.R"
FILESDIR="/work/users/paolo/miles/files"
FIGDIR="/work/users/paolo/miles/figures"

# sourcing
source(file.path(PROGDIR, "script/basis_functions.R"))
source(CFGSCRIPT)

# This is simple script which build the TM90 ensemble mean for a series of
# different project over a common time window and season
# Only one ensemble member is pick for each model (the lowest number available)
# Figures will be different from Davini and D'Andrea 2016 due to slightly 
# different blocking definition
year1 <- 1961
year2 <- 2000
season <- "DJFM"
expid <- "historical"
projects <- c("OBS","CMIP3","CMIP5","CMIP6")
datasets_obs <- c("ERAEXT", "NCEP", "JRA55")

# assign colors
file_idx="D12"
index <- "TM90"
multicol <- c("gray10", "darkorange", "cyan3", "darkorchid3")
names(multicol) <- projects


#----no need to change below this line-----#

# block or U500_block function to select
block_select <- function(project) {
  return(ifelse(project=="CMIP3","U500_Block","Block"))
}

FIGDIR <- file.path(FIGDIR, "TM90_compare")
dir.create(FIGDIR)

# loop on projects
freq_mean <- array(NA, dim = c(144, length(projects)), dimnames <- list(NULL, projects))
for (project in projects) {

  # find datasets
  if (project == "OBS") {
    ens <- exp <- NA
    datasets <- datasets_obs
  } else {
    datasets <- list.files(file.path(FILESDIR,block_select(project), project))
    exp <- expid
    # ad hoc tuning
  }
  print(datasets)
   
  # loop on datasets: if exists open it
  freq_full <- array(NA, dim = c(144, length(datasets)), dimnames <- list(NULL, datasets))
  for (dataset in datasets) {

    # ens selection: pick the first one in the folder
    if (project != "OBS") {
      ens <- list.files(file.path(FILESDIR,block_select(project), project, dataset, exp))[1]
    }

    # build file path
    filename <- file.builder(file.path(FILESDIR), block_select(project), paste0(file_idx,"_Clim"), project, dataset, exp, ens, year1, year2, season)
    if (file.exists(filename)) {
      freq <- ncdf.opener.universal(filename, namevar = index, rotate = "-half")
      freq_full[,which(dataset==datasets)] <- freq$field
    }
  }

  # remove missing datasets and compute ensemble mean
  remove_empty <- which(is.na(freq_full[1,]))
  if (length(remove_empty)!=0) {
    freq_full <- freq_full[,-remove_empty]
  }
  assign(paste0(project,"_freq_full"),freq_full)
  freq_mean[,which(project == projects)] <- rowMeans(freq_full)

}

# remove empty projects
remove_empty <- which(is.na(freq_mean[1,]))
if (length(remove_empty)!=0) {
  freq_mean <- freq_mean[,-remove_empty]
}
real_projects <-  names(freq_mean[1,])


lon <- freq$lon

# function to represent the error ribbons (as std) in an ensemble
error.polygon<-function(lon,ensemble,mode="sd",col="grey",extend=F) {
  if (mode == "sd") {
    sd <- apply(ensemble,1,sd)
    mm <- apply(ensemble,1,mean)
    upper <- mm+sd
    lower <- mm-sd
  } else if (mode == "qnt") {
    lower <- apply(ensemble,1,quantile)["10%",]
    upper <- apply(ensemble,1,quantile)["90%",]
  }
  if (extend) {
    lon<-extend.limits(lon, lon.case=T)
    lower <- extend.limits(lower)
    upper <- extend.limits(upper)
  }
  polygon(c(lon,rev(lon)),c(lower,rev(upper)),col=col, border=col)
}

# function to extend limit for a correct period plotting
extend.limits<-function(field,lon.case=F) {
  return(c(ifelse(rep(lon.case,6), tail(field)-360, tail(field)),
           field,
           ifelse(rep(lon.case,6), head(field)+360, head(field)))
  )

}

#graphical details
fp <- field.details(index)
legend_names <- NULL
lwdline <- 6
xtick <- seq(-60,240,30)
xlabels <- c("60W","30W","0","30E","60E","90E","120E","150E","180E","150W","120W")
ytick <- seq(0,35,5)
ymax <- c(0,35)

# figure plotting: keep very basic with plot
name <- file.path(FIGDIR,paste0(index,"_compare_",season,"_",year1,"-",year2,".pdf"))
legend_names <- NULL
pdf(file = name, width = 14, height = 10, onefile = T, bg = "white", family = "Helvetica")
par(cex.main=2.2, cex.lab=1.8, cex.axis=1.7, mar = c(4, 5, 4, 4), oma = c(1, 1, 2, 1), mgp=c(2.8, 1, 0), xaxs="i", yaxs="i")
plot(lon, lon, type = "n", ylim = ymax, main = paste(fp$title_name,season,paste0(year1,"-",year2)), xlab = "Longitude (deg)", ylab = fp$legend_unit, xlim=c(-90,270), xaxt="n", yaxt="n")
axis(side=1, at=xtick, labels = xlabels)
axis(side=2, at=ytick, labels = ytick)
grid(nx=length(xtick)+1, ny=length(ytick)-1)
for (project in real_projects) {
  points(extend.limits(lon, lon.case=T), extend.limits(freq_mean[,project]), type = "l", lwd = lwdline, lty = 1, col = multicol[project])
  freq_full <- get(paste0(project,"_freq_full"))
  legend_names <- c(legend_names, paste0(project, ": ", length(freq_full[1,])," datasets"))
  error.polygon(lon, freq_full, col = adjustcolor(multicol[project], alpha=0.1), extend=T)
}
legend("topleft", legend = legend_names, lwd = lwdline, lty = 1, col = multicol[real_projects], bg = "white", cex = 1.8)
dev.off()

