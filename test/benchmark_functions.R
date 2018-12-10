source("/home/paolo/MiLES/script/basis_functions.R")
filename="/work/users/paolo/miles/data/zg500/SPHINX/EC-Earth/T799/hab4/zg500_EC-Earth_T799_hab4_fullfile.nc"

field1=ncdf.opener.universal(filename,namevar="zg",tyears=1980:1982,tmonths=c(1,2,12))



for (f in seq(1,4)) {
	field=get(paste0("field",f))
	print(all.equal(power.date.new(field$time),power.date.new2(field$time)))
}

library(rbenchmark)
print(benchmark(power.date.new(field$time),power.date.new2(field$time),replications=10))


field=sample(1:100,2000,replace=T)

n=5
r1=run.mean(field,n=n)
r2=run.mean5(field)
r3=run.mean.general(field,n=n)

print(benchmark(run.mean(field,n=n),run.mean5(field),run.mean.general(field,n=n),replications=200))

