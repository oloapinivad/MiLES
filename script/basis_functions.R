#basis functions
R_LIBLOC=.libPaths()[1]

##########################################################
#------------------------Packages------------------------#
##########################################################

#loadin packages
library("maps",lib.loc=R_LIBLOC)
library("ncdf4",lib.loc=R_LIBLOC)

#check if fast linear fit is operative (after R 3.1): 3x faster than lm.fit, 36x faster than lm
if (exists(".lm.fit")) {lin.fit=.lm.fit} else {lin.fit=lm.fit}

##########################################################
#-----------------Basic functions------------------------#
##########################################################

#normalize a time series
standardize<-function(timeseries)
{
	out=(timeseries-mean(timeseries,na.rm=T))/sd(timeseries,na.rm=T)
	return(out)
}


#detect ics ipsilon lat-lon
whicher<-function(axis,number)
{
	out=which.min(abs(axis-number))
	return(out)
}

#produce a 2d matrix of area weight
area.weight<-function(ics,ipsilon,root=T)
{
field=array(NA,dim=c(length(ics),length(ipsilon)))
if (root==T)
{
        for (j in 1:length(ipsilon))
        {field[,j]=sqrt(cos(pi/180*ipsilon[j]))}
}

if (root==F)
{
        for (j in 1:length(ipsilon))
        {field[,j]=cos(pi/180*ipsilon[j])}
}
return(field)

}

##########################################################
#--------------Time Based functions----------------------#
##########################################################

# to convert season charname to months number
season2timeseason<-function(season)
{
	if (season=="ALL")  {timeseason=1:12}
	if (season=="JJA")  {timeseason=6:8}
	if (season=="DJF")  {timeseason=c(1,2,12)}
	if (season=="MAM")  {timeseason=3:5}
	if (season=="SON")  {timeseason=9:11}
	if (!exists("timeseason")) {stop("wrong season selected!")}
	return(timeseason)
}

#leap year treu/false function
is.leapyear=function(year)
{
	return(((year %% 4 == 0) & (year %% 100 != 0)) | (year %% 400 == 0))
}

power.date<-function(season,ANNO1,ANNO2)
{
#evalute the number of days that will analyze in order
#to create arrays of the needed dimensions

	#create continous calendar
	p1<-as.Date(paste0(ANNO1,"-01-01"))
	p2<-as.Date(paste0(ANNO2,"-12-31"))
	datas=seq(p1,p2,by="day")

	#select only days correspondeing to the needed season
	timeseason=season2timeseason(season)
	month=as.numeric(format(datas,"%m"))
	whichdays=which(month %in% timeseason)

	#create a "season" for continuous time, used by persistance tracking
	seas=whichdays*1; ss=1
	for (i in 1:(length(whichdays)-1)) 
		{
		if (diff(whichdays)[i]>1)  {ss=ss+1}
		seas[i+1]=ss
		}

	#produce a final timeseries of dates
	datas=datas[whichdays]
	dataline=list(day=as.numeric(format(datas,"%d")),month=as.numeric(format(datas,"%m")),year=as.numeric(format(datas,"%Y")),season=seas,data=datas)
	print("Time Array Built")
	print(paste("Length:",length(seas),"days for",season,"season"))
	print(paste("From",datas[1],"to",datas[length(seas)]))

	return(dataline)
}

power.date.no.leap<-function(season,ANNO1,ANNO2)
{
	#apply to power.date object to clean out elements for leap years
	e=power.date(season,ANNO1,ANNO2)
	leap.days=which(e$month==2 & e$day==29)
	dataline.leap=list(day=e$day[-leap.days],month=e$month[-leap.days],year=e$year[-leap.days],season=e$season[-leap.days],data=e$data[-leap.days])
	print("FIXED FOR NO LEAP CALENDAR: Time Array Built")
	print(paste("Length:",length(dataline.leap$season),"days for",season,"season"))
	print(paste("From",dataline.leap$data[1],"to",dataline.leap$data[length(dataline.leap$season)]))
	return(dataline.leap)
}

power.date.30day<-function(season,ANNO1,ANNO2)
{
	#apply to power.date object to clean out elements for leap years
	nmonths=length(season2timeseason(season))
	nyears=as.numeric(ANNO2)-as.numeric(ANNO1)+1
	dd=rep(seq(1,30),nmonths*nyears)
	mm=rep(rep(season2timeseason(season),each=30),nyears)
	#create a "season" for continuous time, used by persistance tracking
	seas=mm*0+1; ss=1
	for (i in 1:(length(mm)-1))
	        {
	        if (diff(mm)[i]>1)  {ss=ss+1}
	        seas[i+1]=ss
	        }
	dataline.30day=list(day=dd,month=mm,season=seas)
	print("SIMPLIFIED CALENDAR FOR 30-day CALENDAR: Time Array Built")
	print(paste("Length:",length(dataline.30day$season),"days for",season,"season"))
	return(dataline.30day)
}

##########################################################
#--------------NetCDF loading function-------------------#
##########################################################

#function to open ncdf files: old but reliable version
ncdf.opener<-function(namefile,namevar,namelon,namelat,rotate=T)
{

#opening file
a=nc_open(namefile)
b=ncvar_get(a,namevar)

#check for 2d or 3d dimensions (presence or not of time dimension)
daily=b
dimensions=length(dim(daily))

#traslating arrays to center on Greenwich
if (dimensions>1)
{
	ics=ncvar_get(a,namelon); ipsilon=ncvar_get(a,namelat); nc_close(a)

	if (rotate==TRUE)
        	{
        	x=vettore=ics
        	x[(length(ics)/2):length(ics)]=vettore[1:(length(ics)/2+1)]
        	x[1:(length(ics)/2-1)]=vettore[(length(ics)/2+2):length(ics)]-360
        	ics=x
        	}
}

if (dimensions==2)
{
        if (ipsilon[1]>ipsilon[2])
                {
                ipsilon=rev(ipsilon)
                b=b[,length(ipsilon):1]
                }
        daily=b


        if (rotate==TRUE)
                {
                matrice=b
                daily[(length(ics)/2):length(ics),]=matrice[1:(length(ics)/2+1),]
                daily[1:(length(ics)/2-1),]=matrice[(length(ics)/2+2):length(ics),]
                }
}


if (dimensions==3)
{
        if (ipsilon[1]>ipsilon[2])
                {
                ipsilon=rev(ipsilon)
                b=b[,length(ipsilon):1,]
                }
        daily=b


        if (rotate==TRUE)
                {
                for (i in 1:length(b[1,1,]))
                        {
                        matrice=b[,,i]
                        daily[(length(ics)/2):length(ics),,i]=matrice[1:(length(ics)/2+1),]
                        daily[1:(length(ics)/2-1),,i]=matrice[(length(ics)/2+2):length(ics),]
                        }
                }
}

#showing array properties
#print(str(daily))

#exporting variables to the main program

if (dimensions>1)
{
	assign("ipsilon",ipsilon, envir = .GlobalEnv)
	assign("ics",ics, envir = .GlobalEnv)
}
daily
}

##########################################################
#--------------Plotting functions------------------------#
##########################################################


#extensive filled.contour function 
filled.contour3 <-
  function (x = seq(0, 1, length.out = nrow(z)),
            y = seq(0, 1, length.out = ncol(z)), z, xlim = range(x, finite = TRUE), 
            ylim = range(y, finite = TRUE), zlim = range(z, finite = TRUE), 
            levels = pretty(zlim, nlevels), nlevels = 20, color.palette = cm.colors, 
            col = color.palette(length(levels) - 1), plot.title, plot.axes, 
            key.title, key.axes, asp = NA, xaxs = "i", yaxs = "i", las = 1, 
            axes = TRUE, frame.plot = axes,mar, ...) 
{
  # modification by Ian Taylor of the filled.contour function
  # to remove the key and facilitate overplotting with contour()
  # further modified by Carey McGilliard and Bridget Ferris
  # to allow multiple plots on one page

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
    else stop("no 'z' matrix specified")
  }
  else if (is.list(x)) {
    y <- x$y
    x <- x$x
  }
  if (any(diff(x) <= 0) || any(diff(y) <= 0)) 
    stop("increasing 'x' and 'y' values expected")
 # mar.orig <- (par.orig <- par(c("mar", "las", "mfrow")))$mar
 # on.exit(par(par.orig))
 # w <- (3 + mar.orig[2]) * par("csi") * 2.54
 # par(las = las)
 # mar <- mar.orig
 plot.new()
 # par(mar=mar)
  plot.window(xlim, ylim, "", xaxs = xaxs, yaxs = yaxs, asp = asp)
  if (!is.matrix(z) || nrow(z) <= 1 || ncol(z) <= 1) 
    stop("no proper 'z' matrix specified")
  if (!is.double(z)) 
    storage.mode(z) <- "double"
  .filled.contour(as.double(x), as.double(y), z, as.double(levels), 
                          col = col)
  if (missing(plot.axes)) {
    if (axes) {
      title(main = "", xlab = "", ylab = "")
      Axis(x, side = 1,...)
      Axis(y, side = 2,...)
    }
  }
  else plot.axes
  if (frame.plot) 
    box()
  if (missing(plot.title)) 
    title(...)
  else plot.title
  invisible()
}

image.scale3 <- function(z,levels,color.palette=heat.colors,colorbar.label="image.scale",
line.label=2,line.colorbar=0,cex.label=1,cex.colorbar=1,colorbar.width=1,...){

 #save properties from main plotting region
 old.par <- par(no.readonly = TRUE)
 mfg.save <- par()$mfg
 old.fig=par()$fig

 #defining plotting region with proper scaling
 #print(old.fig)
 xscal=(old.fig[2]-old.fig[1]); yscal=(old.fig[4]-old.fig[3]); lw=colorbar.width; lp=line.colorbar/100
 new.fig=c(old.fig[2]-0.07*xscal*lw-lp,old.fig[2]-0.03*xscal-lp,old.fig[3]+0.1*yscal,old.fig[4]-0.1*yscal)
 #print(new.fig)

 if (missing(levels)) { levels=seq(min(z),max(z),,12)}
 #fixing color palette
 col=color.palette(length(levels)-1)

 #starting plot
 par(mar=c(1,1,1,1),fig=new.fig,new=TRUE)

 #creating polygons for legend 
 poly <- vector(mode="list", length(col))
 for(i in seq(poly))
  {poly[[i]] <- c(levels[i], levels[i+1], levels[i+1], levels[i])}
 ylim<-range(levels); xlim<-c(0,1)
 plot(1,1,t="n",ylim=ylim, xlim=xlim, axes=FALSE, xlab="", ylab="", xaxs="i", yaxs="i", ...)
 for(i in seq(poly))
   {polygon(c(0,0,1,1), poly[[i]], col=col[i], border=NA)}

 #box, axis and leged
 box()
 axis(4,las=1,cex.axis=cex.colorbar,...)
 mtext(colorbar.label,line=line.label,side=4,cex=cex.label,...)

 #resetting properties for starting a new plot (mfrow style)
 par(old.par)
 par(mfg = mfg.save, new = FALSE)
 invisible()

}

##########################################################
#------------Blocking Tracking Functions-----------------#
##########################################################

#time persistence (used for longitude filter too)
time.persistence<-function(timeseries,persistence=5)
{
rr=rle(timeseries)
rr$values[which(rr$values==1 & rr$length<persistence)]=0
nn=rep(rr$values,rr$length)
return(nn)
}


#blocking 5 days tracking
blocking.persistence<-function(field,persistence=5,time.array)
{

#function for persistence
pers<-function(timeseries,persistence,time.array)
{
        xx=NULL
        for (s in min(time.array$season):max(time.array$season))
                {
                yy=timeseries[which(time.array$season==s)]
                nn=time.persistence(yy,persistence)
                xx=append(xx,nn)
                }
        return(xx)
}

#function for persistence
pers2<-function(timeseries,persistence,time.array)
{     
dd=min(time.array$season):max(time.array$season) 
nn=sapply(dd, function(x) {time.persistence(timeseries[which(time.array$season==x)],persistence)})
xx=c(unlist(nn))
return(xx)
}

# check for etime
if (length(time.array$month)!=length(field[1,1,])) { stop("Wrong time array! Exiting...") }

print("Time filtering...")
newfield=apply(field,c(1,2),function(x) pers2(x,persistence=5,time.array))
newfield=aperm(newfield,c(2,3,1))
print("Mean field...")
meanfield=apply(newfield,c(1,2),mean,na.rm=T)*100


print("Events detection...")
maxdim=max(apply(newfield,c(1,2),function(x) length(rle(x)$length[which(rle(x)$values==1)])))
events=apply(newfield,c(1,2),function(x) c(rle(x)$lengths[which(rle(x)$values==1)],rep(NA,maxdim-length(rle(x)$length[which(rle(x)$values==1)]))))
events=aperm(events,c(2,3,1))
print("Mean Duration...")
duration=apply(events,c(1,2),mean,na.rm=T)
print("Number of Events...")
nevents=apply(events,c(1,2),function(x) length(x[!is.na(x)]))

out=list(track=newfield,percentage=meanfield,duration=duration,events=events,nevents=nevents)
return(out)
}


#large scale extension with further implementation
largescale.extension.if<-function(ics,ipsilon,field)
{
print("Large Scale Extension based on fixed angle")
fimin=30 #southern latitude to be analyzed
fimax=75 #northern latitude to be analyzed
deltaics=diff(ics)[1]
deltaips=diff(ipsilon)[1]
passo=round(5/deltaics)  #horizontal movemenent
vertical=round(2.5/deltaips) #vertical movement
#time=1:length(field[1,1,]) #elements of the length of the dataset
time=which(apply(field,3,max)!=0) #elements length of the dataset (removing no blocked days)

print(paste("Box dimension:",passo*2*deltaics,"° lon x ",vertical*2*deltaips,"° lat"))

short<-function(ics,ipsilon,field,passo,vertical) {
	control=field
	range=which.min(abs(ipsilon-fimin)):which.min(abs(ipsilon-fimax)) #check range for latitude excursion
	#range=range[(1+vertical):(length(range)-vertical)] #reduce range considering border effect
	
	new=rbind(field,field,field) #bind domain for cross-date line
	for (i in 1:length(ics))
		{
		ii=i+length(ics)
		if (!all(new[(ii-passo):(ii+passo),]==0)) #check to speed up
			{
			for (j in range)
				{
				control[i,j]=mean(new[(ii-passo):(ii+passo),(j-vertical):(j+vertical)],na.rm=T)
				}
			}
		}
	control[control>0]=1
	return(control)
}


for (t in time)
{
	if (any(t==round(seq(0,length(field[1,1,]),,11))))
        	{print(paste("--->",round(t/length(field[1,1,])*100),"%"))}
			{field[,,t]=short(ics,ipsilon,field[,,t],passo,vertical)}
}
return(field)
}

#large scale extension
largescale.extension2<-function(ics,ipsilon,field)
{
print("Large Scale Extension based on fixed angle")
deltaics=(ics[20]-ics[19])
deltaips=(ipsilon[20]-ipsilon[19])
passo=round(5/(ics[20]-ics[19]))
vertical=round(2.5/(ipsilon[20]-ipsilon[19]))

print(paste("Box dimension:",passo*2*deltaics,"° lon x ",vertical*2*deltaips,"° lat"))

short<-function(ics,ipsilon,field,passo,vertical)
{
out=field
startipsilon=which.min(abs(ipsilon-30))
estension=round((75-30)/(ipsilon[20]-ipsilon[19]))
new=rbind(field,field,field)
for (i in 1:length(ics))
{ii=i+length(ics)
for (j in startipsilon:(startipsilon+estension))
{
control=mean(new[(ii-passo):(ii+passo),(j-vertical):(j+vertical)],na.rm=T)
if (control>0)
{out[i,j]=1}
}
}
return(out)
}

for (t in 1:length(field[1,1,]))
{
if (any(t==round(seq(0,length(field[1,1,]),,11))))
        {print(paste("--->",round(t/length(field[1,1,])*100),"%"))}
if (all(!is.na(field[,,t])))
{field[,,t]=short(ics,ipsilon,field[,,t],passo,vertical)}
}
return(field)
}

#Longitude filter for minimum extension
longitude.filter<-function(ics,ipsilon,field)
{
print("Longitude filter based on fixed angle")
out=field
deltaics=(ics[20]-ics[19])
startipsilon=which.min(abs(ipsilon-30))
estension=round((75-30)/(ipsilon[20]-ipsilon[19]))
passo=round(15/(ics[20]-ics[19]))

print(paste("Continous longitude contrain",passo*deltaics,"° lon"))

for (t in 1:length(field[1,1,]))
{
if (any(t==round(seq(0,length(field[1,1,]),,11))))
        {print(paste("--->",round(t/length(field[1,1,])*100),"%"))}

new=rbind(field[,,t],field[,,t],field[,,t])
for (j in startipsilon:((startipsilon+estension)))
{
new[,j]=time.persistence(new[,j],persistence=passo)
}
field[,,t]=new[length(ics)+(1:length(ics)),]
}
return(field)
}


##########################################################
#------------EOFs and regims functions-------------------#
##########################################################

eofs<-function(lon,lat,field,neof=4,xlim,ylim,method="SVD",do_standardize=F,do_regression=F)
{
# R tool for computing EOFs based on Singular Value Decomposition ("SVD", default)
# or with the eigenvectors of the covariance matrix ("covariance", slower) 
# If requested, computes linear regressions and standardizes the PCs
# If you want to use the regressions, remember to standardize the PCs
# Take as input a 3D anomaly field.
# Requires "personal" functions area.weight, whicher and standardize

#area weighting, based on the root of cosine
print("Area Weighting...")
ww=area.weight(lon,lat,root=T)
wwfield=sweep(field,c(1,2),ww,"*")

#selection of the box
box=wwfield[whicher(lon,xlim[1]):whicher(lon,xlim[2]),whicher(lat,ylim[1]):whicher(lat,ylim[2]),]
slon=lon[whicher(lon,xlim[1]):whicher(lon,xlim[2])]
slat=lat[whicher(lat,ylim[1]):whicher(lat,ylim[2])]

#transform 3D field in a matrix
new_box=array(box,dim=c(dim(box)[1]*dim(box)[2],dim(box)[3]))

#calling SVD
if (method=="SVD")
{       
        print("Calling SVD...")
        SVD=svd(new_box,nu=neof,nv=neof)
        
        #extracting EOFs (loading pattern), expansions coefficient and variance explained
        pattern=array(SVD$u,dim=c(dim(box)[1],dim(box)[2],neof))
        coefficient=SVD$v
        variance=(SVD$d[1:neof])^2/sum((SVD$d)^2)
        if (do_standardize)
                { coefficient=apply(coefficient,c(2),standardize) }
                else
                { coefficient=sweep(coefficient,c(2),sqrt(variance),"*") }
}

#calling covariance matrix
if (method=="covariance")
{       
        print("Calling eigenvectors of the covariance matrix...")
        covma=cov(t(new_box))
        eig=eigen(covma)
        coef=(t(new_box)%*%eig$vector)[,1:neof]
        pattern=array(eig$vectors,dim=c(dim(box)[1],dim(box)[2],dim(box)[3]))[,,1:neof]
        variance=eig$values[1:neof]/sum(eig$values)
           if (do_standardize)
                { coefficient=apply(coef,c(2),standardize) }
                else
                { coefficient=coef }
}

#linear regressions on anomalies
regression=NULL
if (do_regression)
{       
        print("Linear Regressions (it can takes a while)... ")
        regression=array(NA,dim=c(length(lon),length(lat),neof))
        for (i in 1:neof) {regression[,,i]=apply(field,c(1,2),function(x) coef(lm(x ~ coefficient[,i]))[2])}
}

#preparing output
print("Finalize...")
pattern=list(x=slon,y=slat,z=pattern)
out=list(pattern=pattern,coeff=coefficient,variance=variance,regression=regression)
return(out)
}


regimes<-function(lon,lat,field,ncluster=4,ntime=1000,neof=10,xlim,ylim)
{
# R tool to compute cluster analysis based on k-means.
# Requires "personal" function eofs
# Take as input a 3D anomaly field

#Reduce the phase space with EOFs: use SVD and do not standardize PCs
print("Launching EOFs...")
reducedspace=eofs(lon,lat,field,neof=neof,xlim=xlim,ylim=ylim,method="SVD",do_regression=F,do_standardize=F)

#extract the principal components
PC=reducedspace$coeff

#k-means computation repeat for ntime to find best solution. 
print("Computing k-means...")
regimes=kmeans(PC,ncluster,nstart=ntime,iter.max=100)

#Extract regimes frequencyr and timeseries of occupation
cluster=regimes$cluster
frequencies=regimes$size/dim(field)[3]*100

#Create composites...
print("Creating Composites...")
compose=aperm(apply(field,c(1,2),by,cluster,mean),c(2,3,1))

#prepare output
print("Finalize...")
out=list(cluster=cluster,frequencies=frequencies,regimes=compose)
return(out)
}
