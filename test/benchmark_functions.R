source("/home/paolo/MiLES/script/basis_functions.R")
filename="/work/users/paolo/miles/data/zg500/SPHINX/EC-Earth/T799/hab4/zg500_EC-Earth_T799_hab4_fullfile.nc"

fieldlist=field1=ncdf.opener.universal(filename,namevar="zg",tyears=1980:1982,tmonths=c(1,2,12))

#extract calendar and time unit from the original file
tcal=attributes(fieldlist$time)$cal
tunit=attributes(fieldlist$time)$units

reftime="1970-01-01"
fulltime=as.numeric(fieldlist$time-as.PCICt(reftime,cal=tcal))
#fulltime=as.numeric(etime$data)-as.numeric(etime$data)[1]
TIME=paste0(tunit," since ",reftime,"-01 00:00:00")

#test=daily.anom.mean(field1$lon,field1$lat,field1$field,power.date.new(field$time))

#for (f in seq(1,1)) {
#	field=get(paste0("field",f))
#	print(all.equal(power.date.new(field$time),power.date.new2(field$time)))
#}

#library(rbenchmark)
#print(benchmark(power.date.new(field$time),power.date.old(field$time),replications=10))

#etime=power.date.old(field1$time)
#condition=paste(etime$day,etime$month)
#daily=array(NA,dim=c(length(field1$lon),length(field1$lat),length(unique(condition))))
#for (t in unique(condition)) {
#	 if (sum(t==condition)==1) {
#                       stop("Cannot compute a mean with a single value")
#                }
#        daily[,,which(t==unique(condition))]=rowMeans(field1$field[,,t==condition],dims=2)
#}


rundaily=apply(daily,c(1,2),run.mean5)

field=sample(1:100,2000,replace=T)
rundaily=apply(daily,c(1,2),run.mean5)

n=5
r1=run.mean(field,n=n)
r2=run.mean5(field)
r3=run.mean.general(field,n=n)

print(benchmark(run.mean(field,n=n),run.mean5(field),run.mean.general(field,n=n),replications=200))

#extract calendar and time unit from the original file
tcal=attributes(field1$time)$cal
tunit=attributes(field1$time)$units
etime=power.date.new(field1$time)
fulltime=as.numeric(etime$data)-as.numeric(etime$data)[1]
fulltime2=as.PC
