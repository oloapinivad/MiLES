################################################
#------- MidLatitude Evaluation System --------#
#------------------MiLES v0.8------------------#
#---------Jul 2019, P. Davini, CNR-ISAC--------#
#----------------------------------------------#
################################################

# This is the R wrapper to handle all the MiLES computation
# Based on yaml namelists

library("yaml")
namelist <- yaml.load_file("~/MiLES/namelist/era_test.yml")
config <- yaml.load_file(paste0("~/MiLES/config/config_",namelist$machine,".yml"))

# converting variables
config$PROGDIR <- system(paste("echo",config$PROGDIR),intern=TRUE)
config$OUTPUTDIR <- system(paste("echo",config$OUTPUTDIR),intern=TRUE)

# source functions
script_functions <- c("basis_functions", "block_fast", "block_figures")
for (script_function in script_functions) {
  source(file.path(config$PROGDIR,"script", paste0(script_function, ".R")))
}

# extract the datasets from the namelist
dsets <- namelist$datasets

# set the reference dataset
if (namelist$std_clim) {
  rf <- list(dataset = "ERAI_clim", year1 = 1979, year2 = 2017)
} else {
  rf <- namelist$reference
  dsets[[length(dsets) + 1]] <- rf
}

# temporary hardcoded options
dirout <- config$OUTPUTDIR
refdir <- file.path(config$PROGDIR,"clim")
cfgscript <- file.path(config$PROGDIR,"config","R_config.R")

stop()

# computation block: loop on seasons and datasets
for (season in namelist$seasons) {

  for (i in 1:length(dsets)) {
    ds <- dsets[[i]]
    filein <- file.path("/work/users/paolo/miles/data", namelist$findvar, ds$dataset, 
                        paste(namelist$findvar, ds$dataset, "fullfile.nc", sep = "_"))
    if (any(namelist$options == "block")) {
       miles.block.fast(project = ds$project, dataset = ds$dataset, expid = ds$expid, ens = ds$ens,
                        year1 = ds$year1, year2 = ds$year2, season = season,
                        z500filename = filein, FILESDIR = dirout, doforce = F)
    }
  }
}


# figures block
if (any(namelist$options == "figures")) {

  for (season in namelist$seasons) {

    for (i in 1:length(dsets)) {
      ds <- dsets[[i]]
      if (any(namelist$options == "block")) {
        miles.block.figures(project = ds$project, dataset = ds$dataset, expid = ds$expid, 
                            ens = ds$ens, year1 = ds$year1, year2 = ds$year2, 
                            project_ref = rf$project, dataset_ref = rf$dataset,
                            expid_ref = rf$expid, ens_ref = rf$expid, 
                            year1_ref = rf$year1, year2_ref = rf$year2,
                            varname = "zg", season = season, FIGDIR = dirout, 
                            FILESDIR = dirout, REFDIR = refdir,
                            CFGSCRIPT = cfgscript) 
      }
    }
  }
}


