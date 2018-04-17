# MiLES v0.51
## Mid-Latitude Evaluation System

Oct 2014  - Apr 2018

by P. Davini (ISAC-CNR, p.davini@isac.cnr.it)

Acknowledgements to:
J. von Hardenberg (ISAC-CNR), I. Mavilia (ISAC-CNR), E. Arnone (ISAC-CNR)

------------------------------

## WHAT IS MiLES?

**MiLES** is a tool for estimating properties of Northern Hemisphere mid-latitude climate in Global Climate Models and Reanalysis datasets. It has been originally thought for EC-Earth output and then it has been extended to any climate model or Reanalysis datasets. 
It is based on daily 500hPa Northern Hemisphere geopotential height data and produces NetCDF4 outputs and climatological figures for the chosen time period and season.
Before performing analysis, data are preprocessed interpolated on a common 2.5x2.5 grid using CDO.  
Model data are compared against ECMWF ERA-Interim Reanalysis for a standard period (1979-2017) or with any other MiLES-generated data.

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

3. 	**Z500 Empirical Orthogonal Functions**: Based on EOFs computed by R using SVD.
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

Be aware that this is a free scientific tool in continous development, then it may not be free of bugs. Please report any issue at p.davini@isac.cnr.it

Please refer to **MiLES** specifing which version has been used in the acknowledgment of any publication.

Please cite *"Tibaldi S, Molteni F. 1990. On the operational predictability of blocking. Tellus A 42(3): 343–365, doi:10.1034/j.1600- 0870.1990.t01- 2- 00003.x."* in case you  use the 1D blocking index in any publication.

Similary, please cite *"Davini, P., C. Cagnazzo, S. Gualdi, and A. Navarra, 2012: Bidimensional Diagnostics, Variability, and Trends of Northern Hemisphere Blocking. J. Climate, 25, 6496–6509, doi: 10.1175/JCLI-D-12-00032.1."* in case you use the 2D blocking index in any publication.


----------------

## SOFTWARE REQUIREMENTS

- a. R version >3.0
- b. CDO version > 1.6.5 (1.8 at least for complete GRIB support), compiled with netCDF4
- c. Compiling environment (gcc)

IMPORTANT: there are 5 R packages (ncdf4, maps, PCICt, akima and mapproj) needed to run **MiLES**.
You have to run `Rscript config/installpack.R` as first step in order to install the packages.
If everything runs fine, their installation is performed by an automated
routine that brings the user through the standard web-based installation.
Packages are also included in **MiLES** and can be installed offline.
- _ncdf4_ provides the interface for NetCDF files.
- _maps_ provides the world maps for the plots: (version >= 3.0 )
- _PCICt_ provides the tools to handle 360-days and 365-days calendars (from model data). 
- _akima_ provides the interpolation for map projections.
- _mapproj_ provides a series of map projection that can be used.


If you are aware of other way to implement this 5 passages without using those packages, please contact me.

There are some issues on Mac Os X (10.11 and later at least) related to gfortran. It may happen that 
some packages requires specifically gfortran-4.8 (see here for help:
http://stackoverflow.com/questions/23916219/os-x-package-installation-depends-on-gfortran-4-8)
and that you may find some issue if you install gfortran via MacPorts (see here for help: 
https://stackoverflow.com/questions/29992066/rcpp-warning-directory-not-found-for-option-l-usr-local-cellar-gfortran-4-8)


-----------------

## HOW TO

Before running **MiLES** the 5 above-mentioned R packages should installed.

Two configuration scripts controls the program options:
1. 	`config/config_$MACHINE.sh` controls the properties of your environment. 
	It should be set accordingly to your local configuration.
	It is a trivial configuration, needing only information on CDO/R paths and some folders definition.
    	This also includes the directory tree for your NetCDF files and the expected input files format.
    	It's extremely important that you **create OUR OWN config file**: in this way it will not be overwritten by further pull: two `.tmpl` files for Unix and Mac Os X machines are provided.  
2.	`config/R_config.R` controls the plot properties. If everything is ok, you should not touch this file.
	However, from here you can change in the properties of the plots (as figure size, palettes, axis font, etc.).
	Also output file format and map projection can be specified here if you do not use the wrapper (see later).
	Figures are extremely basic: they can be produced in pdf, png and eps format.

The simplest way to run **MiLES** is executing in bash environment `./wrapper_miles.sh`. 
Options as seasons, which EOFs compute, reference dataset or file output format as well as the map projection to use
can specified at this stage: here below a list of the commands that can be set up
- `dataset_exp` -> this is simply an identifier for your experiments used by MiLES to create files and paths structure.
- `year1_exp` and `year2_exp` -> the years on which MiLES will run. 
- `ens_list` -> ensemble list of experiments from the same dataset: set to "NO" if using a single ensemble. In case of multiple ensemble members an extra ensemble "mean" will be produced by the wrapper only for blocking data.
- `std_clim` -> 1 to use standard ERAI 1979-2017 climatology, 0 for custom comparison. If 0, please specify the dataset you want to compare to with `dataset_ref`, `year1_ref` and `year2_ref`. 
- `seasons` -> specify one or more of the 4 standard seasons using 3 characters (DJF-MAM-JJA-SON). Use "ALL" to cover the full year. Otherwise, use 3 character for each month divided by an underscore to create your own season (e.g. "Jan_Feb_Mar"). This last functionality is under testing.
- `teles` -> A list of one or teleconnection patterns. "NAO" and "AO" for standard EOFs over North Atlantic and Northern Hemisphere. Custorm regions can be specifieds as "lon1_lon2_lat1_lat2".
- `output_file_type` -> pdf, eps or png figures format.
- `map_projection` -> set "no" for standard plot (fast). Use "azequalarea" for polar plots (default). All projection from mapproj R package are supported (but not all of them have been tested).
- `doeof`,`doblock`,`doregime` -> set to true or false in order to run some specific sections.

The chain of scripts will be executed as a sequence.
However, each **MiLES** script can be run autonomously from command line providing the correct sequence of arguments.
R-based script are written as functions and thus can be called inside R if needed.  

* `z500_prepare.sh`. **MiLES** is based on a pre-processing of data. 
This script expects geopotential height data (daily or higher frequency) in a single folder: from v0.5 it SHOULD be able to identify 500hPa data. The code interpolates data on a 2.5x2.5 grid, it selects the NH only and it organizes their structure and their features in order to make them handable by **MiLES**. It produces a single NetCDF4 Zip files with all the data available. A check is performed in order to avoid useless run of the script: if your file is corrupted you can use the `doforce` flag to overwrite it. You can use both geopotential or geopotential height data, the former will be automatically converted. To simplify the analysis by R, the CDO `-a` is used to set an absolute time axis in the data.  

* `Rbased_eof_fast.R` and `Rbased_eof_figures.R`. EOFs are computed using Singular Value Decompositon (SVD) R function by the former script, while the latter provides the figures. EOFs signs for the main EOFs are checked in order to maintain consistency with the reference dataset.

* `blocking_fast.R` and `blocking_figures.R`. blocking analysis is performed by the first R script. The second provides the figures. 
Both the Davini et al. (2012) and the Tibaldi and Molteni (1990) blocking index are computed and plotted by these scripts.

* `regimes_fast.R` and `regimes_figures.R`. Weather regimes analysis is performed by the first R script. The second provides the figures.
It also tries to assign the right weather regimes to its name. However please be aware that it is not always effective.

* `extra_figures_block.R`. This is not called by the wrapper and it provides extra statistics, comparing several experiments with ensemble means, histogram for specific region and Taylor diagrams.

## EXECUTION TIMES

MiLES is pretty fast: on iMac 2017  (MacOS High Sierra 10.13, 3.4 GHz Intel Core i5, 16GB DDR4) 30 years of analysis for a single season takes about
- EOFs: 11 seconds
- Blocking: 57 seconds
- Regimes: 25 seconds
- Figures (together): 20 seconds

Please be aware that issues may arise with large datasets (i.e. larger than 100 years) where the single file approach may be problematic. 
It is reccomended in such cases to split the analysis in different subsets. 

------------

## HISTORY

*v0.51 - Apr 2018
- Improved Netcdf conventions for output files
- Rewritten ncdf.opener function

*v0.5 - Mar 2018*
- Able to detect 500hPa level inside of any geopotential height data
- Improved wrapper with flags to control each section
- Frequency is again plotted on regimes
- Various bug fixing and consolidation
- Improved climatologies (ERAI 1979-2017)

*v0.43 - Feb 2018*
- R-based EOFs script consistent with the MiLES structure
- Rearrange structure of wrapper and config file: now $INDIR is defined in config files (increase portability!)
- Beta support for free month and season selection
- Consistent ensemble members support
- Various bug fixing for NetCDF access
- Improved functions to control path and folders for NetCDF and figures
- Faster daily anomalies computation for weather regimes script
- Variance is again plotted for EOFs
- Template files are provided for Unix and Mac Os X machines

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

