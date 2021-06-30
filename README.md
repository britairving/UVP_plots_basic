# UVP_tools | Basic tools for reading and plotting UVP data
## Introduction
This repository was developed to simplify visualizing Underwater Vision Profiler (UVP) data for the multiple users in the [McDonnell Lab](https://sites.google.com/alaska.edu/mcdonnell/) group. 

This toolbox is very much a work in progress and I welcome feedback and comments.

[UVP_submission_formatting](https://github.com/britairving/UVP_submission_formatting) may also be useful and is a repository that provides scripts to reformat UVP data downloaded from Ecotaxa for submission to data archival websites SeaBASS and BCO-DMO.

## MATLAB dependencies & requirements
This toolbox was created with Matlab R2017b and tested on ....
Matlab version R2016b or higher is required.
Users must have [m_map](https://www.eoas.ubc.ca/~rich/map.html). If any dependent functions are missing in the _utility_ folder please let me know.

## Data format requirements
UVP data must be exported in detailed ODV format from [Ecotaxa](https://ecotaxa.obs-vlfr.fr/part/). See [how-to wiki page](https://github.com/britairving/UVP_plots_basic/wiki/How-to:-exporting-UVP-data-from-Ecotaxa) for instructions. 
***
## Setting up the workflow
All example plots were generated using the data in the _testing_ folder. This data is from a 2015 repeat hydrography cruise along the p16 line in the North Pacific. More information on the cruise and data can be found in the archived datasets [BCO-DMO PAR Dataset](https://www.bco-dmo.org/dataset/787432) and here [BCO-DMO ZOO Dataset](https://www.bco-dmo.org/dataset/787966). Data were formatted for archival using the publicly available repository [UVP_submission_formatting](https://github.com/britairving/UVP_submission_formatting).

### _UVP_workflow.m_ | START HERE!
This script is the starting place where you define the project, filepaths, and variables you want to plot.
Below are the section headers (use code-folding in Matlab Editor to see section headers clearly).
* USER DEFINES _UVP_tools_ filepath. Adds the UVP_plots_basic repository folder and utility folder. 
* USER DEFINES basic plot configuration, project nickname & Ecotaxa filename 
* Read PAR file and calculate parameters from particle abundance and biovolume data
* Read ZOO file with option to limit plotting to fully validated stations
* Read CTD data (this will be implemented soon, but is not ready yet)
* Set data to ZOO OR PAR OR CTD data based on _options.plot_type_ 
* (OPTIONAL) USER hardcodes which variables you want to plot
* Prompt: Determine which variables to plot, if not hardcoded. 
* Remove two deepest depth bins in every profile 
* Prompt: Decide if want to plot sections of data, or the whole range
* PLOT | Plot standard 3 panel total particle abundance, total particle biovolume, and slope of PSD
* PLOT | Plot vertical profiles of different sizes at each station
* Loops through fields select to plot
*   NON-GRIDDED: Waterfall plot - shows evolution of profiles throughout cruise
*   NON-GRIDDED: plot_data_vs_depth(options);
*   GRIDDED: Waterfall plot - shows evolution of data throughout cruise
*   GRIDDED: 2D contour plot 
*   GRIDDED: 3D contour plot 


### _UVP_read_odv_ecotaxa_exported_par.m_
Reads in the UVP PAR file exported from Ecotaxa in detailed ODV format. Data is returned as structure _par_ with information on each field in _par_info_.

```[par_table, par_info] = UVP_read_odv_ecotaxa_exported_par(par_file);```
### _UVP_calculate_PAR_fields.m_
Calculates common variables from the standard UVP particulate output from Ecotaxa. I.e. particulate abundance [#/L] and particulate biovolume [ppm] or [mm^3/L]. New fields will be added here as the repository evolves. 

```[par,par_info] = UVP_calculate_PAR_fields(par_table,par_info);```
### _UVP_read_odv_ecotaxa_exported_zoo.m_
Reads in the UVP ZOO file exported from Ecotaxa in detailed ODV format. Data is returned as structure _zoo_ with information on each field in _zoo_info_.

```[zoo, zoo_info] = UVP_read_odv_ecotaxa_exported_zoo(zoo_file);```
### _Read_hydrographic_cruise_hyfiles.m_
This will likely be replaced in the future, with something to read in a variety of CTD and bottle data and merge with UVP data.

### _UVP_select_paramets_to_plot.m_
Prompts the user to select which variables to plot. If user selects a variable in NxM size such as NSD or biovolume, it asks again which size class or taxa to plot. 
At the end of the script, it will print formatted text to the screen so you don't have to run through this repeatedly. 

EXAMPLE: data.tot_par_abundance (i.e. "total particle abundance") and data.VSD (i.e. "particle biovolume"). Since VSD is an NxM array, you need to define which sizes you want to use.
```
*********************************************************************************
To hardcode these selections, copy and paste this into your UVP_workflow.m script
*********************************************************************************
plots.tot_par_abundance       = struct();
plots.tot_par_abundance.title = {'Abundance (102-203µm)[#/L]'};
plots.tot_par_abundance.clims = [5.35 57.86];
plots.VSD.title = {'Biovolume (102-128_µm)[ppm]'; 'Biovolume (128-161_µm)[ppm]'; 'Biovolume (161-203_µm)[ppm]'};
plots.VSD.clims = [[0.0002 0.00366];[0.000559 0.0094];[0.00058 0.0106]];
```
![](https://github.com/britairving/UVP_plots_basic/blob/master/wiki/PAR_select_variables.png)
![](https://github.com/britairving/UVP_plots_basic/blob/master/wiki/PAR_select_variables_NSD.png)
![](https://github.com/britairving/UVP_plots_basic/blob/master/wiki/ZOO_select_variables.png)
![](https://github.com/britairving/UVP_plots_basic/blob/master/wiki/ZOO_select_variables_biovol.png)
### _uvp_remove_last_two_depth_bins.m_
Removes the two deepest depth bins from each profile. Since the UVP is installed on a CTD rossette there is more time spentat the deepest point in the profile, as the winch is stopped and reverse. Because of this, the counts in the lowest depth bins can be biased to higher counts. 
``[data] = uvp_remove_last_two_depth_bins(data);``

### _select_data_sections.m_
Walks through loading sections for visualizations. Prompts the user to split the dataset given two options; 1) select the number of days or 2) click points on a map. 

### _get_bathymetry.m_
Saves high resolution coastline data using _m_gshh_f.m_ and loads & saves bathymetry using input longitude, latitude.  
Bathymetry is read in from file [project_name '_bathymetry.txt'] if possible.
If bathymetry file is not available, it is loaded using [m_map]((https://www.eoas.ubc.ca/~rich/map.html))'s _m_elev.m_.
How to download your own high resolution bathymetric data
1. EXTRACT XYZ GRID - TOPOGRAPHY OR GRAVITY : https://topex.ucsd.edu/cgi-bin/get_data.cgi
2. COPY AND SAVE WEB QUERY RESULTS TO: [project_name '_bathymetry.txt']


### _UVP_griddata.m_
Grid input Nx1 data into NxM data matricies using Matlab's _griddata.m_ function.
Allows more control of how data is gridded. See script documentation for more.
```
cfg.stophere = 1;         % 1 = stops at the end of the script and plots the gridded data so the user can see how the configuration options influence the visualization
cfg.y_step   = 100;       % [meters] UVP data output in 5m bins from Ecotaxa
cfg.y_start  = 0;         % [meters] depth of first vertical bin
cfg.method   = 'natural'; % interpolation method 'natural','linear' or 'nearest'
cfg.maskXbyZ = 1;         % 1 = masks interpolation over vast horizontal distances, 0 = no mask applied - will interpolate over everything
cfg.scalexy  = 1;         % 1 = rescale axes so the x and y dimensions span a comparable range of values (useful for fixing horizontal streaking issues)
cfg.finite   = 1;         % 1 = ignores all NaNs - may help solve poor gridding coverage
```
## Plots
See [plotting functions](https://github.com/britairving/UVP_plots_basic/wiki/Plotting-functions) wiki page for examples. 

***
# Citation/acknowledgement
If you use this toolbox, please cite as follows:
_Brita Irving, (2021), UVP_plots_basic, GitHub repository, https://github.com/britairving/UVP_plots_basic_
***
