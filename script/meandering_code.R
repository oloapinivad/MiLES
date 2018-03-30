######################################################
#-----Blocking routines computation for MiLES--------#
#-------------P. Davini (Oct 2014)-------------------#
######################################################
#miles.meandering<-function(exp,ens,year1,year2,season,z500filename,FILESDIR,doforce)
#{

PROGDIR="/home/paolo/MiLES"
source(paste0(PROGDIR,"/script/basis_functions.R"))
exp="ERAI"
ens="NO"
year1=1979
year2=2017
season="DJF"
z500filename="/home/paolo/work/miles/Z500/ERAI/Z500_ERAI_fullfile.nc"
FILESDIR="/work/scratch/users/paolo"



#t0
t0<-proc.time()

#setting up time domain
years=year1:year2
timeseason=season2timeseason(season)

#define folders using file.builder function (takes care of ensembles)
savefile1=file.builder(FILESDIR,"MI","MIClim",exp,ens,year1,year2,season)
savefile2=file.builder(FILESDIR,"MI","MIFull",exp,ens,year1,year2,season)

#check if data is already there to avoid re-run
#if (file.exists(savefile1)) {
#	print("Actually requested Meandering Index data is already there!")
#	if (doforce=="true") {
#	 	print("Running with doforce=true... re-run!")
#	} else	{
#	print("Skipping... activate doforce=true if you want to re-run it"); q() 
#	}
#}

#new file opening
nomefile=z500filename
fieldlist=ncdf.opener.time(nomefile,"zg",tmonths=timeseason,tyears=years,rotate="full")
print(str(fieldlist))

#extract calendar and time unit from the original file
tcal=attributes(fieldlist$time)$cal
tunit=attributes(fieldlist$time)$units

#time array to simplify time filtering
etime=power.date.new(fieldlist$time)
totdays=length(fieldlist$time)

#declare variable
Z500=fieldlist$field

########################################################################################
# MEANDERING INDEX Originally from https://github.com/giorgiadicapua/MeanderingIndex  #
########################################################################################

# list of isohypses on which evaluate the MI
isolvls=seq(4800,6200,5)

# reference latitude (60N following Di Capua et al., 2016)
ref_lat=60

# The function calculates the meandering index for a single isohypse (contour of geopotential height) 
# the curviness value is calculated at the chosen reference latitude (60° N suggested for the Polar Front region)
# the script needs a daily field of geopotential height as input, latitudes and longitudes vectors, and the isolevel 
# and returns the value of the meandering index and its latitude

#Original function
MI_lat_ref <- function (longitudes,latitudes,hgt_field,isolvl,ref_lat=ref_lat,verbose=F) {
        if (verbose == TRUE) {
                print("MI_lat_ref")
                print("latitude of reference")
                print(ref_lat)
                print("isolvl")
                print(isolvl)
                print("latitudes")
                print(latitudes)
                print("longitudes")
                print(longitudes)
        }

        #init values
        curv <- 0
        curv_lat <- 0

        # calculate position of the isopleth
        isopleth = contourLines(longitudes,latitudes,hgt_field,nlevels=1,levels=isolvl)

        # longest_iso returns 0 if the isopleth does not exist, 1 if only one isopleth is found or the number of the longest found isopleth 
        i_longest_iso <- longest_iso(isopleth)
        #print(i_longest_iso)

        if (i_longest_iso == 0) {
                # isopleth does not exist
                if (verbose==TRUE) {print("isopleth == 0")}
                curv <- 0; curv_lat <- 0
        } else {
                # isopleth does exist
                longest_isop <- isopleth[[i_longest_iso]]
                # check whether the last and the first point of the isopleth[[i_longest_iso]]y are equal
                #print(abs(longest_isop$y[1]-longest_isop$y[length(longest_isop$y)]))

                if (abs(longest_isop$y[1]-longest_isop$y[length(longest_isop$y)])>3) {
                        if (verbose==TRUE) {print("check whether the last and the first point of the isopleth[[i_longest_iso]]y are equal")}
                        if (verbose==TRUE) {print(paste("isop[1]",longest_isop$y[1],"isop[last]",longest_isop$y[length(longest_isop$y)] ))}

                        # how to solve the problem
                        #lon_translated <-list(lon1=c(61:240,1:60),lon2=c(121:240,1:120),lon3=c(181:240,1:180))
                        # PD: introduce resolution independet more versatile rotation
                        lon0=1:length(longitudes)
                        deltas=seq(length(longitudes)/4,,length(longitudes)/4,3) #this can be further improved
                        lon_translated=lapply(deltas,function(x) {c(tail(lon0,length(lon0)-x),head(lon0,x))}) #use lapply to recursively create list

                        for (ii in 1:length(deltas)) {
                                if (verbose==TRUE) {print("inside the for")}
                                isopleth = contourLines(longitudes,latitudes,hgt_field[lon_translated[[ii]],],nlevels=1,levels=isolvl)
                                i_longest_iso <- longest_iso(isopleth)
                                longest_isop <- isopleth[[i_longest_iso]]

                                if (abs(longest_isop$y[1]-longest_isop$y[length(longest_isop$y)])<3) {
                                        output <- curviness.ref.lat.function(longitudes, latitudes,hgt_field,longest_isop, isolvl, ref_lat)
                                        curv <- output[[1]]
                                        curv_lat <- output[[2]]
                                        if ( curv!=0 & verbose==TRUE) { #PD: add option to avoid plots
                                                contour(longitudes, latitudes,hgt_field[lon_translated[[ii]],], nlevel=10,levels=seq(5000,5900,100),main=curv)
                                                lines(longest_isop$x,longest_isop$y, col=189,lwd=3)
                                                }
                                        break
                                } else {
                                        curv <- 0
                                        curv_lat <- 0
                                }
                        } # end loop over different latitudes

                } else {
                        # calc isopleth curviness at the reference mean latitude and store the original mean latitude 
                        # Isopleths that hit the latitude-borders are disgarded
                        output <- curviness.ref.lat.function(longitudes, latitudes,hgt_field,longest_isop, isolvl, ref_lat)
                        curv <- output[[1]]
                        curv_lat <- output[[2]]
                } # end if which check whether the isopleth is cut at the border of the field
                }
        output_list <- list(curv=curv, curv_lat=curv_lat)
        return(output_list)
}

#Updated function, cleanup by PD (slower!)
MI_lat_ref2 <- function (longitudes,latitudes,hgt_field,isolvl,ref_lat=ref_lat,verbose=F) {
	if (verbose == TRUE) {
    		#print(paste("latitude of reference",ref_lat,"deg"))
        	print(paste("isohypse level",isolvl))
	}

  	#init values
  	curv <- 0
  	curv_lat <- 0

	# introduce resolution independet more versatile rotation	
	# automatic start from no rotation, then try 90 deg every time
	lon0=1:length(longitudes)
        deltas=seq(0,,length(longitudes)/4,4) #this can be further improved
	lon_translated=lapply(deltas,function(x) {c(tail(lon0,length(lon0)-x),head(lon0,x))}) #use lapply to recursively create list
	for (ii in 1:length(deltas)) {

    		# calculate position of the isopleth
		isopleth = contourLines(longitudes,latitudes,hgt_field[lon_translated[[ii]],],nlevels=1,levels=isolvl) #PD is it correct to use longitudes here?
      
      		# longest_iso returns 0 if the isopleth does not exist, 1 if only one isopleth is found or the number of the longest found isopleth 
      		i_longest_iso <- longest_iso(isopleth)
     
		# if isopleth does not exist 
      		if (i_longest_iso == 0) {
	        	if (verbose==TRUE) {print("isopleth == 0")}
			break #exit the loop

		#  # isopleth does exist
		} else {

			longest_isop <- isopleth[[i_longest_iso]]

			#if the longest_isop is shorter than longitudes it cannot be global, exit!
			if (length(longest_isop$x) < length(longitudes) ) { break } 
		
			# check whether the last and the first point of the isopleth[[i_longest_iso]]y are equal (circular)	
			if (verbose==TRUE) {print("Check if the isopleth is circular...")}
                        if (verbose==TRUE) {print(paste("isop[1]",longest_isop$y[1],"isop[last]",longest_isop$y[length(longest_isop$y)] ))}			

			#if isohypse are closed compute the curviness
			if (abs(longest_isop$y[1]-longest_isop$y[length(longest_isop$y)]) < 3) { #PD: why 3? is this resolution dependent?)
				
				#compute curviness
                        	output <- curviness.ref.lat.function(longitudes, latitudes,hgt_field,longest_isop, isolvl, ref_lat)
                                curv <- output[[1]]
                                curv_lat <- output[[2]]
				
				# if you find something plot the line!
                                if ( curv!=0 & verbose==TRUE) { #add option to avoid plots
                                	contour(longitudes, latitudes,hgt_field[lon_translated[[ii]],], nlevel=10,levels=seq(5000,5900,100),main=curv)
                                        lines(longest_isop$x,longest_isop$y, col=189,lwd=3)
                                }
				
				if (verbose==TRUE) {print("Isohypse found!")}	
				# you have your isohypse, break the cycle!
                                break
			}
		if (verbose==TRUE) {print("going to next rotation...")}
                }
	}
        output_list <- list(curv=curv, curv_lat=curv_lat)
        return(output_list)
}



# script for calculate curviness and curviness at a reference lat

curviness.ref.lat.function <- function(lon, lat,hgt_field, isop_old, isolvl, ref_lat, verbose=FALSE) {
	#print(isop_old)
	#print(paste(min(isop_old$x),lon[1],max(isop_old$x),lon[length(lon)]))
	if ( 
	# isopleth should not move outside latitudinal range
	((min(isop_old$y) > lat[2]) & (max(isop_old$y) < lat[length(lat)-1]))
	&
	# isopleth should circle whole globe  
	((min(isop_old$x) == lon[1]) & (max(isop_old$x) == lon[length(lon)]))
	) {
		# first, the mean latitude at which the original isopleth is found must be calculated because we need to know the difference between the reference
		# latitude at which we want to calculate the isopleth
		if (verbose==TRUE) {print("calculating curviness at the latitude of reference")}
		lat_or <- mean(isop_old$y)
        	curv_lat <- lat_or
        	# the new translated laditude must be defined (nothing happens to the longitudes)
        	# pay attention: the reference latitude must be negative in the Southern hemisphere
        	if (lat_or > ref_lat) {
			lat_diff <- abs(lat_or-ref_lat)
	      		lat_new <- lat-lat_diff
	        } else if (lat_or < ref_lat) {
			lat_diff <- abs(lat_or-ref_lat)
	        	lat_new <- lat+lat_diff
		}
	    	# calculate isopleth curviness at the reference latitude
	    	# the isopleth must be recalculated
	    	isopleth_new = contourLines(lon,lat_new,hgt_field,nlevels=1,levels=isolvl)
	    	# the isopleth could not to exist: control
	    	i_longest_iso <- longest_iso(isopleth_new)
	        isop_new <- isopleth_new[[i_longest_iso]]
	        if (i_longest_iso == 0) {
			curv <- 0
		}
		if ( 
		# isopleth should not move outside latitudinal range
		((min(isop_new$y) > lat_new[2]) & (max(isop_new$y) < lat_new[length(lat_new)-1]))
		&
		# isopleth should circle whole globe  
		((min(isop_new$x) == lon[1]) & (max(isop_new$x) == lon[length(lon)]))
		) {
			curv = calculate.isopleth.curviness(isop_new)
		} else #PD: weird behaviour here, add to add it
		{
			curv <- 0
			curv_lat <- 0
		}
		
	} else {
	      if (verbose==TRUE) {print("isopleth is outside latitudinal range or doesn't circle the whole globe ")}
	      curv <- 0
	      curv_lat <- 0 
		          
	}
  
  	output <- list(curv=curv, curv_lat=curv_lat)
  	return(output)
}# end function


# calculate the isopleth curviness function
#
calculate.isopleth.curviness <- function ( isop ) {
	#print("calculate.isopleth.curviness")
	# determine x / y coordinates in [m] for this isopleth
	x_coor   = isop$x
  	y_coor   = isop$y
  	# check correct order (i.e. from west to east)
    	if (x_coor[1] > x_coor[length(x_coor)]) {
	        # invert coordinate vectors
	        x_coor = rev(x_coor)
      		y_coor = rev(y_coor)
        }
    	lat_mean = mean(y_coor)
    	# shift to meter-based grid
    	x_coor   = lon2meter( x_coor, lat_mean)
        y_coor   = lat2meter( y_coor, lat_mean)
        x_coor = duplicate.dateline.x (x_coor, lat_mean)
        y_coor = duplicate.dateline.y (y_coor)
	#plot(x_coor,y_coor)  
	# calculate distance in meters
	distance = isoline.distance(x_coor, y_coor)
	cos_lat       = cos( lat_mean*pi/180. )
	circumference = 2*pi*Earth.Radius * cos_lat
	return( distance / circumference ) # definition of curviness
}

# finde the longest isopleth

longest_iso <- function(isopleth) {
	
	if(length(isopleth)==0) { 
		i_longest_iso = 0
    	} else {
	      	# get index of longest isopleth
	      	if(length(isopleth)>1){
			#print(length(isopleth))
		        distance = 0
        		for(iso in 1:length(isopleth)){
		        	if(length(isopleth[[iso]]$x) > distance){
					distance = length(isopleth[[iso]]$x)
	        			i_longest_iso = iso
		        	} 
	  		}
	   	 } else {
	         #print(length(isopleth))
	         i_longest_iso = 1
		}
    	}
  	return(i_longest_iso)
}

isoline.distance <- function ( x, y ) {
	#print("isoline.distance")
	  
	# calculates the total eucledian distance along a line given by coordinate vectors x and y
	# x and y coordinates should come in units of meters
	# check input:
	if( length(x) != length(y) ) {
		stop("Error: isoline.distance: coordinate vectors have unequal size!")
	}
  	# loop over coordinates, start at i = 2
  	total_distance = 0
  	for(i in 2:length(x)){
		# pythagoras
	        total_distance = total_distance + sqrt( (x[i]-x[i-1])^2 + (y[i]-y[i-1])^2 )
	}
    # done. return total distance
    return(total_distance)
}

lon2meter <- function ( lon, lat ) {
	#print("lon2meter")
	# returns distance in meters for input vector consisting of longitudinal coordinates "lon" at constant latitude "lat" (both in degrees)
	# check input:
	if( length(lat) != 1 ) { 
		stop("Error: lon2meter: single latitude expected!")
  	}
  	if( min(lon) != lon[1] ) { 
	      stop("Error: lon2meter: method expects first coordinate to be the most westerly coordinate!")
    	}
    	cos_lat      = cos( lat*pi/180. )
    	x0           = lon[1] * pi * Earth.Radius / 180. * cos_lat # most westerly point
      	# most westerly point is set to x = 0 meter
      	# thus, most easterly point should have coordinate (close to) circumpherence of Earth at this latitude
      	return ( (lon-lon[1]) * pi * Earth.Radius / 180. * cos_lat);
}

lat2meter <- function ( lat, lat_zero ) {
	# print("lat2meter")
  
	# returns distance in meters of input vector "lat" to latitude "lat_zero" (both in degrees)
	# lat_zero is taken as zero-axis for new meter-based coordinate system
	# check input:
	if( length(lat_zero) != 1 ) { 
		stop("Error: lat2meter: single lat_zero latitude expected!")
  	}
    	return ( (lat-lat_zero) * pi * Earth.Radius / 180. );
}

duplicate.dateline.x <- function ( x, lat_mean ) {
	#print("duplicate.dateline.x")
	  
	# a vector "x_" is returned with size length(x)+1
	# the last coordinate of "x" is duplicated and stored in x_[1]
	# check input:
	if( length(lat_mean) != 1 ) { 
		stop("Error: duplicate.dateline.x: single lat_mean latitude expected!")
	}
  	eq_circ = 2*pi*6367500
  	x_ = array(0,dim=c(length(x)+1))
  	x_[1:length(x)] = x
  	x_[length(x)+1] = cos(lat_mean*pi/180)*eq_circ
  	return(x_)
}

duplicate.dateline.y <- function ( x ) {
	#print("duplicate.dateline.y")
	  
	# a vector "x_" is returned with size length(x)+1
	# the first coordinate of "x" is duplicated and stored in last coordinate of x_
	x_ = array(0,dim=c(length(x)+1))
	x_[1:length(x)] = x
    	x_[length(x)+1] = x[1]
    	return(x_)
} 

#debugs
#prova=MI_lat_ref(ics,ipsilon,Z500[,,1],5850,60,verbose=T)
#isopleth = contourLines(ics,ipsilon,Z500[,,1],nlevels=1,levels=5850)
#contour(ics,ipsilon,Z500[,,1],nlevels=1,levels=5850)
#MI_list=sapply(isolvls,function(x) {MI_lat_ref(ics,ipsilon,Z500[,,1],x,60,verbose=F)})
#MI_list=apply(Z500,c(3),function(y) {sapply(isolvls,function(x) {MI_lat_ref(ics,ipsilon,y,x,60,verbose=F)})})
#library("microbenchmark")
#microbenchmark(sapply(isolvls,function(x) {MI_lat_ref2(ics,ipsilon,Z500[,,t],x,ref_lat,verbose=F)}),sapply(isolvls,function(x) {MI_lat_ref(ics,ipsilon,Z500[,,t],x,ref_lat,verbose=F)}))

t1=proc.time()-t0
print(t1)

#Running the real code
MI_lat=MI_value=1:totdays*NA
for (t in 1:totdays) {
	if (any(t==round(seq(0,totdays,,21))))	{
		print(paste("--->",round(t/totdays*100),"%"))
	}
	#sapply on the MI function
	MI_list=sapply(isolvls,function(x) {MI_lat_ref2(ics,ipsilon,Z500[,,t],x,ref_lat,verbose=F)})
	MI_value[t]=max(unlist(MI_list[1,]))
	MI_lat[t]=unlist(MI_list[2,])[which.max(unlist(MI_list[1,]))]
}




tf=proc.time()-t1
print(tf)


##########################################################
#------------------------Save to NetCDF------------------#
##########################################################

#saving output to netcdf files
#print("saving NetCDF climatologies...")

#which fieds to plot/save
fieldlist=c("MI","MI_lat")
full_fieldlist=c("MI","MI_lat")

# dimensions definition
fulltime=as.numeric(etime$data)-as.numeric(etime$data)[1]
TIME=paste(tunit," since ",year1,"-",timeseason[1],"-01 00:00:00",sep="")
LEVEL=50000
x <- ncdim_def( "Lon", "degrees_east", 0, longname="Longitude")
t1 <- ncdim_def( "Time", TIME, 1, unlim=T, calendar=tcal, longname="Time")
t2 <- ncdim_def( "Time", TIME, fulltime,unlim=T, calendar=tcal, longname="Time")


for (var in fieldlist)
{
        #name of the var
	if (var=="MI") 
		{longvar="Meandering Index"; unit=""; field=mean(MI_value); full_field=MI_value}
	if (var=="MI_lat")   
		{longvar="Meandering Index Latitude"; unit="deg"; field=mean(MI_lat); full_field=MI_lat}

        #variable definitions
	var_ncdf=ncvar_def(var,unit,list(x,t=t1),-999,longname=longvar,prec="single",compression=1)
	full_var_ncdf=ncvar_def(var,unit,list(x,t=t2),-999,longname=longvar,prec="single",compression=1)
        
	assign(paste0("var",var),var_ncdf)
	assign(paste0("full_var",var),full_var_ncdf)
        assign(paste0("field",var),field)
	assign(paste0("full_field",var),full_field)
}

#Climatologies Netcdf file creation
print(savefile1)
namelist1=paste0("var",fieldlist)
nclist1 <- mget(namelist1)
ncfile1 <- nc_create(savefile1,nclist1)
for (var in fieldlist) {
        # put variables into the ncdf file
	ncvar_put(ncfile1, fieldlist[which(var==fieldlist)], get(paste0("field",var)), start = c(1, 1),  count = c(-1,-1))
}
nc_close(ncfile1)

#Fullfield Netcdf file creation
print(savefile2)
namelist2=paste0("full_var",full_fieldlist)
nclist2 <- mget(namelist2)
ncfile2 <- nc_create(savefile2,nclist2)
for (var in full_fieldlist) {
        # put variables into the ncdf file
        ncvar_put(ncfile2, full_fieldlist[which(var==full_fieldlist)], get(paste0("full_field",var)), start = c(1, 1),  count = c(-1,-1))

}
nc_close(ncfile2)

#}

#blank lines
cat("\n\n\n")

# REAL EXECUTION OF THE SCRIPT 
# read command line
#args <- commandArgs(TRUE)

# number of required arguments from command line
#name_args=c("exp","ens","year1","year2","season","z500filename","FILESDIR","PROGDIR","doforce")
#req_args=length(name_args)

# print error message if uncorrect number of command 
#if (length(args)!=0) {
#    if (length(args)!=req_args) {
#        print(paste("Not enough or too many arguments received: please specify the following",req_args,"arguments:"))
#        print(name_args)
#    } else {
# when the number of arguments is ok run the function()
#        for (k in 1:req_args) {assign(name_args[k],args[k])}
#        source(paste0(PROGDIR,"/script/basis_functions.R"))
#        miles.meandering(exp,ens,year1,year2,season,z500filename,FILESDIR,doforce)
#    }
#}
