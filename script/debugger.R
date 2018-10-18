rm(list=ls())
PROGDIR="/home/paolo/MiLES"
FILESDIR="/work/users/paolo/miles/files"
source(file.path(PROGDIR,"script/basis_functions.R"))
source(file.path(PROGDIR,"script/meandering_functions.R"))
dataset="ERAI"
expid="no"
ens="no"
year1=1979
year2=2008
season="DJF"
z500filename="/work/users/paolo/miles/Z500/ERAI/Z500_ERAI_fullfile.nc"

#setting up time domain
years=year1:year2
timeseason=season2timeseason(season)

#load file
fieldlist=ncdf.opener.universal(z500filename,namevar="zg",tmonths=timeseason,tyears=years,rotate="full")
print(str(fieldlist))

#extract calendar and time unit from the original file
tcal=attributes(fieldlist$time)$cal
tunit=attributes(fieldlist$time)$units

#time array to simplify time filtering
etime=power.date.new(fieldlist$time)
totdays=length(fieldlist$time)

#declare variable
Z500=fieldlist$field

# list of isohypses on which evaluate the MI
isolvls=seq(4900,6200,5)

# reference latitude (60N following Di Capua et al., 2016)
ref_lat=60

#Running the real code
MI_lat=MI_value=1:totdays*NA

