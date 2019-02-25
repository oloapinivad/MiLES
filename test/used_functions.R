# 
rm(list=ls())
library("mvbutils")

listfile=list.files("script",pattern="*.R")
for (l in listfile) {
	print(l)
	if (l=="extra_figures_block.R") next
	source(paste0("script/",l))
}

foodweb()

