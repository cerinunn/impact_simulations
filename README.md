Welcome to impact_simulations. This code was used to:

a) View the seismograms of the artificial impacts on the Moon (1969 - 1972).

b) Set up 3D models for seismic simulations for the Moon.

c) View the simulations.

d) Create the other figures for a forthcoming paper.

----
See also https://doi.org/10.5281/zenodo.10631693

The Zenodo repository includes the parameter files and outputs (in NetCDF files).

----
DEPENDENCIES

This was the code I used to set up my local python environment - called postprocessing.
I used Miniconda - which will need to be set up first.


conda create -n postprocessing python=3.8  

conda activate postprocessing

conda install jupyter pyyaml matplotlib basemap pandas

conda install -c conda-forge obspy netcdf4 imageio

conda install conda-forge::seaborn

conda install anaconda::xarray

pip install --find-links=https://irfpy.irf.se/sdist/ irfpy.util -U
    # Moon basemap  

pip install --no-index --find-links=https://irfpy.irf.se/sdist irfpy.planets -U


----
Navigate to the directory and at the command prompt type "jupyter notebook"

NOTE I've recently had this fail with a blank page:
http://localhost:8888/tree

It also has an error message on the command line:
"Skipped non-installed server(s)"

A workaround is to do a reset on the open browser page:
(Command, Shift and the 'R' on a Mac)

See this site for an introduction to Jupyter notebooks
https://jupyter.org/try-jupyter/notebooks/?path=notebooks/Intro.ipynb


----
OVERVIEW OF THE NOTEBOOKS

#####################

postprocessing_observations.ipynb

View the observations of artificial impacts.

Figures 1, 2 and Table 1

#####################

postprocessing_simulations.ipynb

View the Simulations (some figures also include the observations)

Figures 4, 5, 6, 9, 10, 11, 12, 13

#####################

TauP_plots.ipynb

Create the TauP models and view the velocity models

Figures 3, 7, 8

#####################

quantitative_estimates.ipynb

Quantitative comparisons between M-2 and observations

Estimates of Rise Time and Characteristic Decay Time

Figures 14, 15, 16

#####################

view_all_simulations.ipynb

View all simulations from 0-60 degrees. Similar to Figure 13 in the paper, but for all simulations.

#####################

view_artificial_observations.ipynb

Alternative views of the Observations - with and without P, PS and S phases

#####################
comparison.ipynb

Compare observation at each distance with simulation M-2

Includes the observations used in Onodera et al., 2024.

#####################

slope_analysis.ipynb

View the gradient of seismic envelope

#####################

Fourier_parameter.ipynb

Assess the Fourier parameter

#####################

build_Moon_scatter.ipynb

Build heterogeneities for the models

#####################

build_Moon_surface_model.ipynb

Build Moho and surface topography

#####################

combine_netcdf.ipynb

Combine the NetCDF files

#####################

OTHER FILES

input_files/ImpactParameters.csv

Timing, location and trajectory of the Artificial Impacts.
