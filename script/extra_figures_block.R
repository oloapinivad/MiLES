
source("/home/ms/it/ccpd/MiLES/config/config.R")
source("/home/ms/it/ccpd/MiLES/script/basis_functions.R")
BASEDIR="/home/ms/it/ccpd/scratch/miles"
DATADIR=file.path(BASEDIR,"files")
FIGDIR=file.path(BASEDIR,"figures","ExtraMultiModel")
dir.create(FIGDIR,recursive=T)


filebuilder<-function(DATADIR,dataset,ens,year1,year2,season) {
	if (ens=="NO") {
		filedir=file.path(DATADIR,dataset,"Block",paste0(year1,"_",year2),season)
		filename=paste0("BlockClim_",dataset,"_",year1,"_",year2,"_",season,".nc")	
	} else {
		filedir=file.path(DATADIR,dataset,ens,"Block",paste0(year1,"_",year2),season)
		filename=paste0("BlockClim_",dataset,"_",ens,"_",year1,"_",year2,"_",season,".nc")
	}
	return(paste0(filedir,"/",filename))
}

area.weight2<-function(ics,ipsilon) {
		field=array(NA,dim=c(length(ics),length(ipsilon)))
        re=6371220
        for (j in 2:length(ipsilon)) {
			dx=(ics[2]-ics[1])*re*pi/180*cos(pi/180*ipsilon[j])
    		dy=(ipsilon[j]-ipsilon[j-1])*re*pi/180
   			field[,j]=dy*dx
	}
    field[,1]=field[,length(ipsilon)]
	return(field)

}

field.details<-function(field) {
   
	if (field=="BlockEvents") {
        color_field=palette1; color_diff=palette2
        lev_field=seq(0,27,3); lev_diff=seq(-10.5,10.5,1); lev_hist=c(0,16)
        legend_unit="Blocked Days (%)"; title_name="Blocking Events frequency:";
    }

    if (field=="DurationEvents") {
        color_field=palette0; color_diff=palette2
        lev_field=seq(5,11.5,.5); lev_diff=seq(-2.1,2.1,.2); lev_hist=c(0,8) 
        legend_unit="Duration (days)"; title_name="Duration of Blocking Events:";
    }

    if (field=="NumberEvents") {
        color_field=palette0; color_diff=palette2
        lev_field=seq(0,100,10); lev_diff=seq(-42.5,42.5,5); lev_hist=c(0,60)
        legend_unit=""; title_name="Number of Blocking Events:";
    }

out=list(color_field=color_field,color_diff=color_diff,lev_field=lev_field,lev_diff=lev_diff,lev_hist=lev_hist,legend_unit=legend_unit,title_name=title_name)
return(out)
}



ensfinder<-function(dataset) {
	if (dataset=="ERAI") {
	    enslist=c("NO")
	}

	if (dataset=="S4" | dataset=="S5") {
	    enslist=c("00","01","02","03","04","05","06","07","08","09")
	}

return(enslist)
}
	
year1=1982
year2=2016
season="DJF"
dataset_ref="ERAI"
datasets=c(dataset_ref,"S4","S5")
SECTORS=c("Euro","Azores","Greenland","FullPacific")
variables=c("BlockEvents","DurationEvents")
KOL=c("black","darkgreen","blue","darkorange","red","violet","grey50","black")

for (variable in variables) {

	fp=field.details(variable)

for (dataset in datasets) {

	enslist=ensfinder(dataset)
	freq_full=array(NA,dim=c(144,37,length(enslist)))
	for (ens in enslist) {

		print(paste(dataset,ens,variable))
		freq=ncdf.opener(filebuilder(DATADIR,dataset,ens,year1,year2,season),namevar=variable,namelon="Lon",namelat="Lat",rotate="no")
		freq_full[,,which(ens==enslist)]=freq
	}
	assign(paste(variable,dataset,sep="_"),freq_full)
	
}


print("2D polar plots...")

nn=1.5
cex.letter=4
PROJ="azequalarea"
ORIENT=c(90,0,0)
lettering=c("(a)","(b)","(c)","(d)","(e)","(f)","(g)","(h)","(i)","(l)")

name=paste(FIGDIR,"/",variable,"_MultiModelComparison_",year1,"_",year2,"_",season,".pdf",sep="")
pdf(file=name,width=5*length(datasets)*2,height=10*2,onefile=T,bg="white",family='Helvetica')
par(mfrow=c(2,length(datasets)),mar=c(6,6,6,6),oma=c(1,1,5,8),mgp=c(4,1,0),cex.main=6)

for (dataset in datasets) {
	field=apply(get(paste(variable,dataset,sep="_")),c(1,2),mean,na.rm=T)
	polar=proj.plot(ics,ipsilon,field,proj=PROJ,orient=ORIENT,lmin=30,npoints=100)
	filled.contour3(polar$x,polar$y,polar$z,main=paste(toupper(dataset)),xlab="",ylab="",levels=fp$lev_field,color.palette=palette0,axes=F,xlim=polar$xlim,ylim=polar$ylim)
	mtext(lettering[which(dataset==datasets)],line=1.5,adj=0,cex=cex.letter)
	if (dataset==datasets[1]) {mtext(fp$title_name,side=3,line=.5,outer=TRUE,cex=5,font=2)}
	box();
	proj.addland(proj=PROJ,orient=ORIENT)
}
image.scale3(field,color.palette=palette0,levels=fp$lev_field,colorbar.width=1.4,colorbar.label=fp$legend_unit,line.colorbar=-0.5,line.label=7,cex.colorbar=3.5,cex.label=3)

for (dataset in datasets) {
	field=apply(get(paste(variable,dataset,sep="_")),c(1,2),mean,na.rm=T)
	field_ref=apply(get(paste(variable,dataset_ref,sep="_")),c(1,2),mean,na.rm=T)
	polar=proj.plot(ics,ipsilon,field-field_ref,proj=PROJ,orient=ORIENT,lmin=30,npoints=100)
	filled.contour3(polar$x,polar$y,polar$z,main=paste(toupper(dataset),"-",dataset_ref),xlab="",ylab="",levels=fp$lev_diff,color.palette=palette2,axes=F,xlim=polar$xlim,ylim=polar$ylim)
	polar=proj.plot(ics,ipsilon,field_ref,proj=PROJ,orient=ORIENT,lmin=30,npoints=100)
	contour(polar$x,polar$y,polar$z,levels=fp$lev_field[-1],add=T,lwd=2,drawlabels=F)
	mtext(lettering[which(dataset==datasets)+5],line=1.5,adj=0,cex=cex.letter)
	box(); #addland.stereo.lambert(color="gray40",lat.min=l.min)
	proj.addland(proj=PROJ,orient=ORIENT)

	if (dataset==datasets[1]) {
		for (SECTOR in SECTORS) {
			if (SECTOR=="FullPacific") {SECTOR="FullPacific2"}
			v=sector.details(SECTOR)
			rect=c(NA,NA)
			rect=rbind(rect,cbind(rep(v$lons[1],length(v$lats[1]:v$lats[2])),v$lats[1]:v$lats[2]))
			rect=rbind(rect,cbind(v$lons[1]:v$lons[2],rep(v$lats[2],length(v$lons[1]:v$lons[2]))))
			rect=rbind(rect,cbind(rep(v$lons[2],length(v$lats[1]:v$lats[2])),v$lats[2]:v$lats[1]))
			rect=rbind(rect,cbind(v$lons[2]:v$lons[1],rep(v$lats[1],length(v$lons[2]:v$lons[1]))))
			rect=rect[2:length(rect[,1]),]
			rr=mapproject(rect[,1],rect[,2],projection=PROJ,orientation=ORIENT)
			lines(rr$x,rr$y,lwd=6,lty=1,col="red")
		}
	}
}

image.scale3(field,color.palette=palette2,levels=fp$lev_diff,colorbar.width=1.4,colorbar.label=fp$legend_unit,line.colorbar=-0.5,line.label=7,cex.colorbar=3.5,cex.label=3)

dev.off()

print("Histogram regional sectors...")

regional=regional2=array(NA,dim=c(length(SECTORS),length(datasets)))
for (SECTOR in SECTORS) {
	v=sector.details(SECTOR)
	for (dataset in datasets) {
        field=get(paste(variable,dataset,sep="_"))
        field_mean=apply(field[v$lonssel,v$latssel,,drop=F],3,mean,na.rm=T)
        FREQ=mean(field_mean,na.rm=T); SDFREQ=sd(field_mean,na.rm=T)
        regional[which(SECTOR==SECTORS),which(dataset==datasets)]=FREQ
        regional2[which(SECTOR==SECTORS),which(dataset==datasets)]=SDFREQ

    }
}

name=paste(FIGDIR,"/",variable,"_RegionalBlocking_",year1,"_",year2,"_",season,".pdf",sep="")
pdf(file=name,width=20,height=20,onefile=T)
par(mfrow=c(2,2),mar=c(8,6,6,6),oma=c(6,1,5,1),cex.axis=2.5,cex.main=4.5,cex.lab=2.5,mgp=c(4,1,0))
spacing=c(0.5,0.2,0.2)
for (k in 1:length(SECTORS)) {
	v=sector.details(SECTORS[k])
	barplot(regional[k,],names.arg=datasets,col=KOL,main=paste(v$name),density=100,space=spacing,las=3,ylim=fp$lev_hist,ylab=fp$legend_unit)
	if (k==1) {mtext(paste0(fp$title_name),side=3,line=.5,outer=TRUE,cex=5,font=2)}
	abline(h=regional[k,1],lwd=2,lty=3)
	segdist=seq(0.5,length(datasets)-0.5,1)+cumsum(spacing)
	segments(segdist,regional[k,]-regional2[k,],segdist,regional[k,]+regional2[k,],lwd=3)
	epsilon = .1
	segments(segdist-epsilon,regional[k,]-regional2[k,],segdist+epsilon,regional[k,]-regional2[k,],lwd=3)
	segments(segdist-epsilon,regional[k,]+regional2[k,],segdist+epsilon,regional[k,]+regional2[k,],lwd=3)
	mtext(lettering[k],line=1.5,adj=0.02,cex=cex.letter)
}
dev.off()

}

print("Taylor diagrams...")

variable="BlockEvents"
fp=field.details(variable)

selectname="nonzero"; normalize=T
field_ref=get(paste(variable,dataset_ref,sep="_"))
WW0=area.weight(ics,ipsilon,root=F)
if (selectname=="nonzero") {nonzero=which(field_ref!=0)}
if (selectname=="fulldomain") {nonzero=1:length(field_ref)}
if (normalize) {
		SIGMA0=1
		} else {
        SIGMA0=weighted.sd(field_ERA[nonzero],WW0[nonzero])
		}
radius=round(SIGMA0*2); delta=radius/8
xminplot=0; xmaxplot=radius
yminplot=0; ymaxplot=radius

name=paste(FIGDIR,"/",variable,"_TaylorDiagram_",year1,"_",year2,"_",season,".pdf",sep="")
pdf(file=name,width=24,height=24,onefile=T)
par(mar=c(8,8,8,4),oma=c(1,1,1,1),cex.lab=3.5,cex.main=4,cex.axis=3.5,mgp=c(4,1,0))

plot(SIGMA0,0,main=paste("Taylor diagram for",year1,"-",year2,season,"\n",fp$title_name,"climatology"),xlim=c(xminplot,xmaxplot),ylim=c(yminplot,xmaxplot),axes=F,xlab="Standard Deviation",ylab="StandardDeviation",col="black",cex=5,pch=16)
abline(v=c(xminplot),h=c(yminplot),lwd=3)
axis(1,at=seq(xminplot,xmaxplot,delta),labels=seq(xminplot,xmaxplot,delta),tick=F)
axis(2,at=seq(yminplot,ymaxplot,delta),labels=seq(yminplot,ymaxplot,delta),tick=F)

TOTESSE=NULL

PCHS=c(17,rep(16,length(datasets)-1))
for (dataset in datasets)
        {
		field_full=get(paste(variable,dataset,sep="_"))
		for (ens in 1:length(ensfinder(dataset))) {
			field=field_full[,,ens,drop=F]
        	ERRE=weighted.cor(field[nonzero],field_ref[nonzero],WW0[nonzero])
        	if (normalize) {
        		SIGMAF=weighted.sd(field[nonzero],WW0[nonzero])/weighted.sd(field_ref[nonzero],WW0[nonzero])
        	    } else {
        	    SIGMAF=weighted.sd(field[nonzero],WW0[nonzero])
       		}
        	errezero=1
        	ESSE=(4*(1+ERRE)^4)/((SIGMAF+1/SIGMAF)^2*(1+errezero)^4)
        	TOTESSE=c(TOTESSE,ESSE)
        	points(SIGMAF*cos(acos(ERRE)),SIGMAF*sin(acos(ERRE)),cex=3,pch=PCHS[which(dataset==datasets)],bg=KOL[which(dataset==datasets)],col=KOL[which(dataset==datasets)],lwd=3)
	}
}

cexlines=4
aspdisp=radius/40

#details correlation coefficient
ll=c(0,radius)
for (rr in c(.2,.4,.6,.8,.9,.95,.99)) {
	points(ll*cos(acos(rr)),ll*sin(acos(rr)),type="l",lty=3,col="dark blue",lwd=3)
	text(xmaxplot*cos(acos(rr)),ymaxplot*sin(acos(rr)),rr,col="dark blue",cex=2.5,pos=4)
}

#external circle
rrext=seq(0,1,.001);  ll=radius
points(ll*cos(acos(rrext)),ll*sin(acos(rrext)),type="l",lty=1,lwd=cexlines)

# rms circles
for (ll in seq(radius/8,radius,radius/8)) {
	rr=seq(-1,1,0.001)
	sx=rr[which.min(abs(SIGMA0+ll*cos(acos(rr))))]
	#dx=(5-ll^2)/2
	dx=(SIGMA0^2+radius^2-ll^2)/(2*SIGMA0)
	dx=rr[which.min(abs(SIGMA0+ll*cos(acos(rr))-dx))]
	rrnew=seq(sx,dx,0.001)
	points(SIGMA0+ll*cos(acos(rrnew)),ll*sin(acos(rrnew)),type="l",lty=3,col="dark green",lwd=cexlines)
	if (ll<(radius*3/4))
		text(SIGMA0,ll,ll,col="dark green",cex=2.5,pos=3)
}

#standard deviation circles
rr=seq(0,1,0.001)
for (ll in seq(0,radius,delta)) {
	points(ll*cos(acos(rr)),ll*sin(acos(rr)),type="l",lty=3,lwd=cexlines)
}

#legend
legend(xminplot+0.05,ymaxplot-aspdisp,toupper(datasets),pch=PCHS,col=KOL,bg="white",cex=4)

dev.off()

name=paste(FIGDIR,"/","Duration_vs_events_",year1,"_",year2,"_",season,".pdf",sep="")
pdf(file=name,width=20,height=20,onefile=T)
par(mar=c(6,6,6,6),oma=c(1,1,5,1),mfrow=c(2,2),cex.lab=3,cex.main=3.5,cex.axis=3,mgp=c(4,1,0))

for (SECTOR in SECTORS) {
	v=sector.details(SECTOR)
	plot(1,1,main=paste(v$name),xlim=c(6.2,8),ylim=c(0,60),type="n",xlab="Avg Duration (days)",ylab="Avg Number of events")
	grid()
	PCH=16

	if (SECTOR==SECTOR[1]) {mtext(paste0("Blocking Events duration vs. Number of events"),side=3,line=.5,outer=TRUE,cex=4,font=2)}

	for (dataset in datasets)
       {
	   evfield=get(paste("NumberEvents",dataset,sep="_"))
	   durfield=get(paste("DurationEvents",dataset,sep="_"))
       for (ens in 1:length(ensfinder(dataset)))
                {
                EV=mean(evfield[v$lonssel,v$latssel,ens],na.rm=T)
                DUR=mean(durfield[v$lonssel,v$latssel,ens],na.rm=T)
                points(DUR,EV,cex=3,pch=PCH,col=KOL[which(dataset==datasets)],bg=KOL[which(dataset==datasets)])
                }
		EV=mean(evfield[v$lonssel,v$latssel,],na.rm=T)
        DUR=mean(durfield[v$lonssel,v$latssel,],na.rm=T)
		points(DUR,EV,cex=6,pch=PCH+2,col=KOL[which(dataset==datasets)],bg="black",lwd=3)
        text(DUR+0.02,EV,toupper(paste(dataset)),pos=4,cex=3)
	}

}
dev.off()

