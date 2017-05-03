#R_LIBLOC=.libPaths()[1]
#dir.create(R_LIBLOC,recursive=T)

#preference is for web-based installation
web=1

if (web==1)
{
# web based installation
packages=c("spam","maps","fields","ncdf4")
new.packages <- packages[!(packages %in% installed.packages()[,"Package"])]

if (length(new.packages)==0) {print("All packages are there, no need to install anything")}
if (length(new.packages)!=0) {print(paste("Installing",length(new.packages),"R packages... follow on-screen instruction"))}

for (pack in new.packages)
{
print(paste("Installing... ",pack))
install.packages(pack,repos="http://cran.irsn.fr",type="source")
}
}


if (web==0)
{
#manual installation
packages=c("spam_1.0-1.tar.gz","maps_2.3-9.tar.gz","fields_7.1.tar.gz","ncdf4_1.13.tar.gz")
PROGDIR<-Sys.getenv(c("PROGDIR"))
R_PACKDIR=paste0(PROGDIR,"/R_packages")
for (pack in new.packages)
{
packageinstall=paste(R_PACKDIR,"/",pack,sep="")
install.packages(packageinstall,repos=NULL,type="source")
}
}

