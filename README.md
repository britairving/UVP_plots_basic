# UVP_tools | Basic tools for reading and plotting UVP data
## Introduction
This repository was developed to simplify visualizing Underwater Vision Profiler (UVP) data for the multiple users in the [McDonnell Lab](https://sites.google.com/alaska.edu/mcdonnell/) group. This toolbox is a work in progress and we hope these tools are useful for others and welcome feedback and comments.

[UVP_submission_formatting](https://github.com/britairving/UVP_submission_formatting) may also be useful and is a repository that provides scripts to reformat UVP data downloaded from Ecotaxa for submission to data archival websites SeaBASS and BCO-DMO.

## MATLAB dependencies & requirements
This toolbox was created with Matlab R2017b and tested on ....
Matlab version R2016b or higher is required.
Users must have [m_map](https://www.eoas.ubc.ca/~rich/map.html)

## Data format requirements
UVP data must be exported in detailed ODV format from [Ecotaxa](https://ecotaxa.obs-vlfr.fr/part/).
***
## Setting up the workflow
All example plots were generated using the data in the _testing_ folder. This data is from a 2015 repeat hydrography cruise along the p16 line in the North Pacific. More information on the cruise and data (not appropriate formatting for this toolbox) can be found here [BCO-DMO PAR Dataset](https://www.bco-dmo.org/dataset/787432) and here [BCO-DMO ZOO Dataset](https://www.bco-dmo.org/dataset/787966).

### _UVP_workflow.m_ | START HERE!
This script is the starting place where you define the project, filepaths, and variables you want to plot.

### _UVP_read_odv_ecotaxa_exported_par.m_

### _UVP_calculate_PAR_fields.m_

### _UVP_read_odv_ecotaxa_exported_zoo.m_

### _Read_hydrographic_cruise_hyfiles.m_

### _UVP_select_paramets_to_plot.m_

### _get_bathymetry.m_
Saves high resolution coastline data using _m_gshh_f.m_ and loads & saves bathymetry using input longitude, latitude.  
Bathymetry is read in from file [project_name '_bathymetry.txt'] if possible.
If bathymetry file is not available, it is loaded using [m_map]((https://www.eoas.ubc.ca/~rich/map.html))'s _m_elev.m_.
How to download your own high resolution bathymetric data
1. EXTRACT XYZ GRID - TOPOGRAPHY OR GRAVITY : https://topex.ucsd.edu/cgi-bin/get_data.cgi
2. COPY AND SAVE WEB QUERY RESULTS TO: [project_name '_bathymetry.txt']

## Plots
### _plot_uvp_multipanel.m_
Plot multipanel figure with all profiles of fields.
``plot_uvp_multipanel(par,par_info,{'tot_par_abundance' 'tot_par_biovolume' 'slope_b'},options);

### _plot_uvp_NSD.m_
Plot vertical profiles of different sizes at each station
``plot_uvp_NSD(par,par_info,options);

### _plot_data_vs_depth.m_
Waterfall plot of data that is intended to highlight vertical evolution of each profile.
### _UVP_griddata.m_
Grid input Nx1 data into NxM data matricies using Matlab's _griddata.m_ function.
Allows more control of how data is gridded. See script documentation for more.

### _plot_gridded_data_vs_depth_waterfall.m_
### _plot_2d_gridded_transect.m_

### _plot_3d_gridded_transect.m_

%    evolution of each profile.

***
# Citation/acknowledgement
If you use this toolbox, please cite as follows:
_Brita Irving, (2021), UVP_basic_plots, GitHub repository, https://github.com/britairving/UVP_basic_plots/_
***
