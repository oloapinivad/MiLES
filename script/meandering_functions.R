######################################################
#-----Meandering routines computation for MiLES------#
#--------G. Di Capua & P. Davini (Apr 2018)----------#
######################################################

#############################################################################
# --------MEANDERING INDEX Originally from Di Capua et al. (2016) ----------#
# ---------------- functions updated by P. Davini --------------------------#
#############################################################################


# Updated from the the original MI_lat_ref. It includes:
# 1. Versatile longitudine rotation independent on resolution and/or grid
# 2. curv and curv_lat error values for debugging: 
#    NA (no assignation, i.e. something gone bad) ,
#    0 (no length found),
#    -1 (not circulating around the globel),
#    -2  (outside latitudinal range) 
#    -3 (multiple longitudinal rotation failing)  
# 3. First isopleth detection and following longitudinal rotation are integrated in a single loop.
MI.fast <- function (longitudes,latitudes,hgt_field,isolvl,ref_lat=ref_lat,verbose=F) {

	printv<-function(value) {if (verbose) {print(value)} }
        printv(paste("latitude of reference",ref_lat,"deg"))
        printv(paste("isohypse level",isolvl))

        #init values
        curv=NA; curv_lat=NA

        # PD: introduce resolution independet more versatil2e rotation       
        # PD: automatic start from no rotation, then try 90 deg every time
        lon0=1:length(longitudes)
        deltas=seq(0,,length(longitudes)/4,4) #4 steps rotation, this can be further improved
        for (ii in 1:length(deltas)) {

                # avoid translation if it is not necessary, save time
                if (ii==1) {
			lon_translated=lon0
			longitudes2=longitudes
		} else {
			lon_translated=c(tail(lon0,length(lon0)-deltas[ii]),head(lon0,deltas[ii]))
			longitudes2=c(tail(longitudes,length(lon0)-deltas[ii]),head(longitudes+360,deltas[ii]))
		}

                # calculate position of the isoplet
                isopleth = contourLines(longitudes2,latitudes,hgt_field[lon_translated,],nlevels=1,levels=isolvl) #PD is it correct to use longitudes here?

                # longest_iso returns 0 if the isopleth does not exist, 1 if only one isopleth is found or the number of the longest found isopleth 
                i_longest_iso <- longest.iso.fast(isopleth)

                # if isopleth does not exist 
                if (i_longest_iso == 0) {
                        printv("isopleth == 0")
                        curv=0; curv_lat=0
                        break #exit the loop

                #  isopleth does exist
                } else {
                        longest_isop <- isopleth[[i_longest_iso]]

                        # if you find something plot the line!
                        if (verbose==TRUE) { #add option to avoid plots
                                contour(longitudes2, latitudes,hgt_field[lon_translated,], nlevel=10,levels=seq(4800,6200,100),main=curv)
                                lines(longest_isop$x,longest_isop$y, col=189,lwd=3)
                                }

                        #if the longest_isop is shorter than longitudes it cannot be global, exit!
                        if ((min(longest_isop$x) != longitudes2[1]) | (max(longest_isop$x) != longitudes2[length(longitudes2)])) {
                                curv=-1; curv_lat=-1;
                                printv("isopleth doesn't circle the whole globe ")
                                break
                        }

                        #if isop is outside of the reasonable latitude range, exit!
                        if ((min(longest_isop$y) < ipsilon[2]) | (max(longest_isop$y) > ipsilon[length(ipsilon)-1])) {
                                curv=-2; curv_lat=-2;
                                printv("isopleth is outside latitudinal range")
                                break
                        }

                        # check whether the last and the first point of the isopleth[[i_longest_iso]]y are equal (circular)     
                        printv("Check if the isopleth is closed around dateline...")
                        printv(paste("isop[1]",longest_isop$y[1],"isop[last]",longest_isop$y[length(longest_isop$y)] ))

                        #if isohypse are closed compute the curviness
                        if (abs(longest_isop$y[1]-longest_isop$y[length(longest_isop$y)]) < diff(latitudes)[1]) { #PD: one grid point

                                #compute curviness
				printv("Closed isohypse found!")
                                output <- ref.lat.fast(longitudes2,latitudes,hgt_field[lon_translated,],longest_isop, ref_lat,verbose=verbose)
                                curv <- output[[1]]
                                curv_lat <- output[[2]]

                                # you have your isohypse, break the cycle!
                                break
                        }
                printv("going to next rotation...")
                }
        }
        if (ii==4 && is.na(curv)==TRUE) {
               curv=-3; curv_lat=-3
               printv("No closed hysoipse found...")
        }
        output_list <- list(curv=curv, curv_lat=curv_lat)
        return(output_list)
}

# Updated from curviness.ref.lat.function2.g: It includes:
# 1. Abrupt simplification: a simple tranlation of latitudes is done without recomputing the longes isoline
ref.lat.fast <- function(lon, lat,hgt_field, isop_old,ref_lat, verbose=FALSE) {
	printv<-function(value) {if (verbose) {print(value)} }
        printv(paste(round(min(isop_old$y),2),lat[2],round(max(isop_old$y),2),lat[length(lat)-1]))
        printv(paste(min(isop_old$x),lon[1],max(isop_old$x),lon[length(lon)]))
        printv("calculating curviness at the latitude of reference")
        # first, the mean latitude at which the original isopleth is found must be calculated because we need to know the difference between the reference
        # latitude at which we want to calculate the isopleth
        lat_or <- mean(isop_old$y)
        curv_lat <- lat_or

        # extra options for different reference latitudes!
        #if (is.null(ref_lat)) {ref_lat=lat_or}
        #if (is.na(ref_lat)) {ref_lat=weighted.mean(isop_old$y,cos(isop_old$y*pi/180))}

        # the new translated laditude must be defined (nothing happens to the longitudes)
        # pay attention: the reference latitude must be negative in the Southern hemisphere     
        lat_new=lat+(ref_lat-lat_or)
        isop_new=isop_old

        #PD: we don't need to recompute the isophlets, just update the latitudes with the new reference
        isop_new$y=isop_new$y+(ref_lat-lat_or)

        #compute curviness
        curv = curviness.fast(lon,lat,isop_new)

        output <- list(curv=curv, curv_lat=curv_lat)
        return(output)
}

# Updated from calculate.isopleth.curvines. It includes:
# 1. Removed dependence on lon2meter and lat2meter
# 2. Removed dependence on dateline x/y functions
curviness.fast <- function ( lon,lat,isop ) {
        
	# determine x / y coordinates in [m] for this isopleth
        lon   = isop$x
        lat   = isop$y
        
	# check correct order (i.e. from west to east)
        if (lon[1] > lon[length(lon)]) {
                # invert coordinate vectors
                lon = rev(lon)
                lat = rev(lat)
        }
        lat_mean = mean(lat)
        cos_lat = cos( lat_mean*pi/180. )

        # shift to meter-based grid
        # PD: remove dependenc on lon2meter and lat2metere functions
        x_coor =  (lon-lon[1]) * pi * Earth.Radius / 180. * cos_lat
        y_coor =  (lat-lat_mean) * pi * Earth.Radius / 180.

        # PD: extend dimensions, remove dependence on function
        x_coor = c(x_coor, cos_lat*2*pi*Earth.Radius)
        y_coor = c(y_coor, y_coor[1])

        distance = isoline.distance.fast(x_coor, y_coor)
        circumference = 2*pi*Earth.Radius * cos_lat
        return( distance / circumference ) # definition of curviness
}

# Updated from longest.iso (slower but cleaner using lapply)
longest.iso.fast <- function(isopleth) {

        if (length(isopleth)<2 ) {
                i_longest_iso = length(isopleth)
        } else {
                i_longest_iso = which.max(unlist(lapply(isopleth,function(x) length(x$x))))
                #i_longest_iso = which.max(unlist(lapply(isopleth,function(x) length(x[[2]]))))
                }
        return(i_longest_iso)
}

# Updated from isoline.distance: use vectorization to improve speed
isoline.distance.fast <- function ( x, y ) {

        # calculates the total eucledian distance along a line given by coordinate vectors x and y
        # x and y coordinates should come in units of meters
        # check input:
        if( length(x) != length(y) ) {
                stop("Error: isoline.distance: coordinate vectors have unequal size!")
        }
        #total_distance=sum(sapply(2:length(x), function(i) {sqrt( (x[i]-x[i-1])^2 + (y[i]-y[i-1])^2 )}))
        total_distance=sum(sqrt((x[-1]-x[-length(x)])^2+(y[-1]-y[-length(y)])^2))

        # done. return total distance
	#print(total_distance)
        return(total_distance)
}

