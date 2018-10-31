#zero order test to check integrity of netcdf opener

rm(list=ls())
PROGDIR="/home/paolo/MiLES"
FILESDIR="/work/users/paolo/miles/files"
source(file.path(PROGDIR,"script/basis_functions.R"))


filenames=c("/work/users/paolo/miles/Z500/ERAI/Z500_ERAI_fullfile.nc")
seasons=c("DJF","JJA")

for (filename in filenames)  {
	for (season in seasons)  {

		year1=2000
		year2=2001

	#setting up time domain
	years=year1:year2
	timeseason=season2timeseason(season)

	#load file
	fieldlist=ncdf.opener.universal(filename,namevar="zg",tmonths=timeseason,tyears=years,rotate="full",verbose=false)


	}
}


