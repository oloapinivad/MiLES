#config name: create your own config file for your machine.
machine=wilma

# part of the code you want to run. Possible options are:
#"eofs": eofs parts
#"block": blocking
#"regimes": regimes
#"meandering": meandering index
#"figures": figures for all selected parts
#"forcedata": recompute the fullfile
#"forceanl" : recompute the analysis of each part
options="block forceanl"

# dataset-experiment-ensemble identificators: important for the folder structure and the naming.
# if you have more than one run or experiments of the same model use this variable to distinguish them
# Beware that $ens_list is a list of ensemble on which you can run a loop
# set also years 
dataset_exp=ERAI
#project_exp=CMIP5
#expid_exp=T255
#ens_exp=lab2
year1_exp=1979
year2_exp=2017

# std_clim flag: this is used to choose which climatology compare with results
# or with a user specified one: standard climatology is ERAINTERIM 1979-2017
# if std_clim=true ERAINTERIM 1979-2017 is used
# if std_clim=false a MiLES-generated different climatology can be specified
std_clim=true

# only valid if std_clim=false
dataset_ref=NCEP
#project_ref=CMIP5
#expid_ref=historical
#ens_ref=r1
year1_ref=1986
year2_ref=2005

# please specify one or more of the 4 standard seasons using 3 characters. 
# std_clim is supported for these 4 seasons only. 
# To analyse the whole year use "ALL"
# Beta: now you can define your own season putting together 3-character string for each consecutive month you 
# want to include, for example "Jan_Feb_Mar".
#seasons="Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec"
#seasons="Nov_Dec_Jan_Feb_Mar_Apr"
seasons="DJF"


# select which teleconnection pattern EOFs you want to compute
# "NAO": the 4 first  EOFs of North Atlantic, i.e. North Atlantic Oscillation as EOF1
# "AO" : the 4 first EOFs of Northern Hemispiere, i.e. Arctic Oscillation as EOF1 
# "PNA": the 4 first EOFs of North Pacific, i.e. Pacific North American Pattern as EOF1 (beta)
# "lon1_lon2_lat1_lat2" : custom regions for EOFs: beware that std_clim will be set to false!
#teles="NAO AO PNA"
teles="PNA NAO"
#tele="-50_20_10_80"

#variables to analyze: don't change unless you know what are you doing!
varname=zg
level=500
