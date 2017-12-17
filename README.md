# MiLES v0.4
## Mid-Latitude Evaluation System

Dec 2017

by P. Davini (ISAC-CNR, p.davini@isac.cnr.it)

Contributors: 	
J. von Hardenberg (ISAC-CNR), I. Mavilia (ISAC-CNR)

------------------------------

## WHAT IS MiLES?

**MiLES** is a tool for estimating properties of Northern Hemisphere mid-latitude climate in Global Climate Models
and Reanalysis datasets. It has been originally thought for EC-Earth output and then
it has been extended to any model data. It is based uniquely on R and CDO.
It works on daily 500hPa geopotential height data and produces NetCDF4 outputs and climatological figures 
for the chosen time period (over the for standard 4 seasons) in 3 possible output formats.
Map projection for plots can specified as well.
Before performing analysis, data are preprocessed interpolated on a common 2.5x2.5 grid using CDO.  
Model data are compared against ECMWF ERA-INTERIM reanalysis for a standard period (1979-2016) or with any 
other MiLES-generated data

Current version includes:

1. 	**1D Atmospheric Blocking**: *Tibaldi and Molteni (1990)* index for Northern Hemisphere.
        Computed at fixed latitude of 60N, with delta of -5,-2.5,0,2.5,5 deg, fiN=80N and fiS=40N.
        Full timeseries and climatologies are provided in NetCDF4 Zip format.

2. 	**2D Atmospheric blocking**: following the index by *Davini et al. (2012)*.
	It is a 2D version of *Tibaldi and Molteni (1990)* for Northern Hemisphere
	atmospheric blocking evaluating meridional gradient reversal at 500hPa.
	It includes also Meridional Gradient Index and Blocking Intensity index
	and Rossby wave orientation index, computing both Instantenous Blocking and Blocking Events frequency.
	Blocking Events definition allows the estimation of the blocking duration.
	A supplementary Instantaneous Blocking index with the GHGS2 conditon is also evaluted. 
	Full timeseries and climatologies are provided in NetCDF4 Zip format.

3. 	**Z500 Empirical Orthogonal Functions**: Based on CDO "eofs" function.
	First 4 EOFs for North Atlantic (over the 90W-40E 20N-85N box) and Northern Hemisphere (20N-85N).
	North Atlantic Oscillation, East Atlantic Pattern, and Arctic Oscillation are thus computed. 
	Figures showing linear regression of PCs on monthly Z500 are provided.
	PCs and eigenvectors, as well as the variances explained are provided in NetCDF4 Zip format.

4.	**North Atlantic Weather Regimes (beta)**: following k-means clustering of 500hPa geopotential height.
	4 weather regimes over North Atlantic (80W-40E 30N-87.5N) are evaluted using 
	anomalies from daily seasonal cycle. North Atlantic first 4 EOFs are computed to reduce 
	the phase-space dimension and then k-means clustering using Hartigan-Wong algorithm with k=4 is computed. 
	Figures report patterns and frequencies of occurrence. NetCDF4 Zip data are saved.
	*Only 4 regimes and DJF season is supported so far.*

----------------

## MAIN NOTES

Please be aware that this is free tool in continous development, then it may not be 
free of bugs. Please report any issue at p.davini@isac.cnr.it

Please refer to **MiLES** specifing which version has been used in the acknowledgment of any publication.

Please cite *"Tibaldi S, Molteni F. 1990. On the operational predictability of blocking. 
Tellus A 42(3): 343–365, doi:10.1034/j.1600- 0870.1990.t01- 2- 00003.x."*
in case you  use the 1D blocking index in any publication.

Please cite *"Davini, P., C. Cagnazzo, S. Gualdi, and A. Navarra, 2012:
Bidimensional Diagnostics, Variability, and Trends of Northern Hemisphere Blocking.
J. Climate, 25, 6496–6509, doi: 10.1175/JCLI-D-12-00032.1."*
in case you use the 2D blocking index in any publication.


----------------

## SOFTWARE REQUIREMENTS

* a. R version >3.0
* b. CDO version > 1.6.5 (1.8 at least for GRIB support), compiled with netCDF4
* c. Compiling environment (gcc)

IMPORTANT: there are 5 R packages (ncdf4, maps, PCICt, akima and mapproj) needed to run **MiLES**.
You have to run "Rscript config/installpack.R" as first step in order to install the packages.
If everything runs fine, their installation is performed by an automated
routine that brings the user through the standard web-based installation.
Packages are also included in **MiLES** and can be installed offline.
- "ncdf4" provides the interface for NetCDF files.
- "maps" provides the world maps for the plots: it needs to be at least v3.0
- "PCICt" provides the tools to handle 360-days and 365-days calendars (from model data). 
- "akima" provides the interpolation for map projections.
- "mapproj" provides a series of map projection that can be used.


If you are aware of other way to implement this 5 passages without using those packages, please contact me.

The installation of some packages requires specifically gfortran-4.8: there is an issue known on 
Mac OS X (10.11 and later at least) which requires a few turnarounds. See here for help:
http://stackoverflow.com/questions/23916219/os-x-package-installation-depends-on-gfortran-4-8

-----------------

## HOW TO

Before running **MiLES** R packages should installed (see above).

Two configuration scripts controls the program options:
1. 	*config/config_$MACHINE.sh* controls the properties of your environment. 
	It should be set accordingly to your local configuration. 
	It is a trivial configuration, needing only information on CDO/R paths and some folders definition.
2.	*config/config.R* controls the plot properties. If everything is ok, you should not touch this file.
	However, from here you can change in the properties of the plots (as figure size, palettes, axis font, etc.).
	Also output file format and map projection can be specified here if you do not use the wrapper (see later).
	Figures are extremely basic: they can be produced in pdf, png and eps format.

The simplest way to run **MiLES** is executing in bash environment "./wrapper_miles.sh". 
Options as seasons, which EOFs compute, reference dataset or file output format as well as the map projection to use
can specified at this stage: here below a list of the commands that can be set up
- "ECMWF" -> this is to call the preprocessing of Z500 ad hoc for ECMWF data structure
- "dataset_exp" -> this is simply an identifier for your experiments used by MiLES to create files and paths: if you have multiple ensemble members you should
  distinguish them from here
- "ens_list" -> ensemble list of experiments from the same dataset: set to "NO" if using a single ensemble. Ensemble "mean" will be produced by the wrapper.
- "std_clim" -> 1 to use standard ERAI 1979-2016 climatology, 0 for custom comparison. 
- "seasons" -> specify one or more of the 4 standard seasons using 3 characters 
- "tele" -> "NAO" and "AO" for standard EOFs over North Atlantic and Northern Hemisphere. Custorm regions can be specifieds as "lon1_lon2_lat1_lat2". 
- "output_file_type" -> pdf, eps or png figures format
- "map_projection" -> set "no" for standard plot (fast). Use "azequalarea" for polar plots. All projection from mapproj R package are supported.


The chain of scripts will be executed as a sequence. You can comment the script you do not need.
However, each **MiLES** script can be run autonomously from command line providing the correct sequence of arguments.
R-based script are written as functions and thus can be called inside R if needed.  

* "z500_prepare.sh". **MiLES** is based on a pre-processing of data. 
This script expects daily 500hPa geopotential height data in a single folder: it interpolates data on a 2.5x2.5 grid,
it selects the NH only and it organizes their structure and their features in order to make them handable by **MiLES**.
It produces a single NetCDF4 Zip files with all the data avaialble. A check is performed in order to avoid useless run of 
the script: if your file is corrupted you need to remove it by hand.
You can use both geopotential or geopotential height data, the former will be automatically converted.   
To simplify the analysis by R, the CDO "-a" is used to set an absolute time axis in the data.  
A new version of this file ("ecmwf_z500_prepare.sh") is provided to run with ECMWF data structure: it is clearly evident that you can think about personalize
this script in order to tailor it on your data structure.

* "eof_fast.sh" and "eof_figure.R". EOFs are computed using CDO in bash environment by the former script, while the latter
provides the figures with an R script. EOFs signs for the main EOFs are checked in order to maintain consistency with the reference dataset.

* "blocking_fast.R" and "blocking_figures.R". blocking analysis is performed by the first R script. The second provides the figures. 
Both the Davini et al. (2012) and the Tibaldi and Molteni (1990) blocking index are computed and plotted by these scripts.

* "regimes_fast.R" and "regimes_figures.R". Weather regimes analysis is performed by the first R script. The second provides the figures.
It also tries to assign the right weather regimes to its name. However please be aware that it is not always effective.

* "extra_figures_block.R". This is not called by the wrapper and it provides extra statistics, comparing several experiments with ensemble means, histogram for specific region and Taylor diagrams.

------------

## HISTORY

*v0.42 - Dec 2017*
- Inclusion of extra blocking diagnostics (Taylor diagrams, Duration-Events plots, histograms, etc.)
- Ensemble mean for blocking outputs
- Ensemble member support for blocking routine
- Bug fixing for calendar handling 
- 10-day blocking events as new output
- ECMWF data structure support
- Updated climatology (1979-2016)
- Support for Grib files

*v0.41 - Jul 2017*
- Plot bug fixing

*v0.4 - June 2017*
- Tibaldi and Molteni (1990) blocking index is now computed by blocking_fast.R
- Weather regimes based on k-means clustering over North Atlantic is now available.
- Reformulation of input Z500 files, now based on a single NetCDF file: to handle 360-days 
  and 365-day calendar package PCICt is now required.
- Polar projection support: requires mapproj and akima packages. 
- Figures updates and various bug fixing.
- Re-written wrapper to provide dynamic comparison of datasets

*v0.31 - May 2017*
- Comparison of EOFs and Blocking figures with any other MiLES-generated dataset.
- Beta-version of sign-check for main EOFs.
- Reformulation: each script is made by R function + can be run from command line.
- Change folder structure to simplify portability.
- Code consolidation and folder/variable name normalization.

*v0.3 - May 2017*
- Blocking Events definition by Davini et al. (2012) now avaiable.
- Removed dependencies from fields and spam R packages.
- Support for figures format in png, pdf or eps - by J. von Hardenberg.
- Removed dependencies on R-files saving blocking data (using now NetCDF).
- Blocking timeseries available in NetCDF.
- NetCDF4 Zip for blocking output files.
- Support for different model calendar: 30-day, Gregorian and No-Leap-Year.
- ~36x faster linear regression for EOFs (.fit.lm function).
- new ~2x faster largescale.extension.if() function.
- Improved speed in blocking for long timeseries: ~2.5x faster for 30years (predeclaration of arrays).
- Minor bugs in axis legends (removal of image.plot).
- Readme in markdown format.

*v0.2 - Apr 2017*
- Support for Arctic Oscillation.
- External unique configuration file.
- Psuedo-universal adaptability to any model data.
- Automated script for R package installing.
- Adaptation to geopotential/geopotential height data.
- Climatological blocking data are stored in NetCDF.
- Now on GitHUB.

*v0.11 - Mar 2015*

- Update to fast blocking (Blocking2-scheme) computation.

*v0.1 - Oct 2014*

- EOFs and 2D Blocking calculation.
- Basic functions implemented.
- Support for NetCDF4.
- Support for 4 standard season (DJF,MAM,JJA,SON).
- ERAINTERIM comparison via netCDF files.
- Parallelization Z500 extraction.
- Png outputs from PDF.

-----------------

