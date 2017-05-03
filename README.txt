MiLES v0.2
Mid-Latitude Evaluation System

by P. Davini - ISAC-CNR
May 2017

p.davini@isac.cnr.it

------------------------------

WHAT IS MiLES?

MiLES is a tool for estimating properties of mid-latitude climate originally thought
for EC-Earth output and then extended to any model data.
It works on daily 500hPa geopotential height data and it produces climatological figures 
for the chosen time period. Data are interpolated on a common 2.5x2.5 grid.  
Model data are compated against ECMWF ERA-INTERIM reanalysis for a standard period (1989-2010).
It supports analysis for the 4 standard seasons.

Current version include:
1. 	2D Instantaneous Blocking: based on Davini et al (2012)
	It is a 2D version of Tibaldi and Molteni (1990) for Northern Hemisphere
	atmospheric blocking index evaluating meridional gradient reversal at 500hPa.
	Includes also Meridional Gradient Index and Blocking Intensity index.
	Full data are saved in R format and climatologies are provided also in NetCDF format

2. 	Z500 North Atlantic EOFs. Based on CDO "eofs" function.
	First 4 EOFs for NAO (over the 90W-40E 20N-85N box) and AO (20N-85N).
	Figures of linear regression are provided.
	PCs and eigenvectors, as well as the variances explained are provided in NetCDF format.

----------------

MAIN NOTES

Please be aware that this is free tool in continous development, then it may not be 
free of bugs. 

Please refer to MiLES specifing which version has been used in the acknowledgment of any publication.

Please cite "Davini, P., C. Cagnazzo, S. Gualdi, and A. Navarra, 2012:
Bidimensional Diagnostics, Variability, and Trends of Northern Hemisphere Blocking.
J. Climate, 25, 6496–6509, doi: 10.1175/JCLI-D-12-00032.1."
if you want to use the blocking index in any publication.


----------------

SOFTWARE REQUIREMENTS

a. R version >3.0
b. CDO version > 1.6.5, compiled with netCDF4
c. Convert, GhostScript and ImageMagick to produce png figures (additional)

There are 4 R packages needed to run MiLES. If everything runs fine, their installation is
performed by an automated routine that brings the user through the standard web-based installation. 
However, in case of any issues, they are included in the package and can be installed offline.

-----------------

HOW TO

MiLES is simply run executing in bash environment the "wrapper_miles.sh"

However, before running config_$MACHINE.sh should be set accordingly to your configuration.
It is extremely basic, needing only information on CDO/R paths and some folders specification.
R packages installation is the only tricky part, but they should be handled automatically by the config files.

MiLES is based on a pre-processing of data, performed by z500_prepare.sh
This script expects daily 500hPa geopotential height data in a single folder: it interpolates data on a 2.5x2.5 grid,
it select the NH only and it organizes their structure and their features in order to make them handable by MiLES.
You can use both geopotential or geopotential height data, the former will be automatically converted. 

After that, blocking analysis is performed by two different R script (blocking_fast.R an blocking_figures.R)

EOFs are produced with CDO via the eof_fast.sh script and plots are produced by the eof_figure.R script

Figures are extremely basic: they are provided in pdf but they can be converted in png or any other format with the
tool based on covert in the last lines of the wrapper_miles.sh

------------

HISTORY

v0.2 - May 2017
-support for Arctic Oscillation
-external configuration file
-psuedo-universal adaptability to any model data
-automated script for R package installing
-adaptation to geopotential/geopotential height data
-climatological blocking data are stored in NetCDF
-Now on GitHUB

v0.11 - Mar 2015

-Update to fast blocking (Blocking2-scheme) computation.

v0.1 - Oct 2014

-EOFs and 2D Blocking calculation
-Basic functions implemented
-Support for NetCDF4
-Support for 4 standard season (DJF,MAM,JJA,SON)
-ERAINTERIM comparison via netCDF files
-Parallelization Z500 extraction
-Png outputs from PDF

-----------------

TO BE DONE:

-Code consolidation
-Improve parallelizazion on seasons for R
-Find formalism for cross-year months
-Free month selection?
-Additional climatologies comparison and reorganize structure
-Polar plot support
-Additional packages: JLI and TM90 index


