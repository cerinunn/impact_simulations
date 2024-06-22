Welcome to impact_simulations. This code was used to: 

a) View the seismograms of the artificial impacts on the Moon (1969 - 1972). 

b) Set up 3D models for seismic simulations for the Moon.

c) View the simulations. 

d) Create the other figures for a forthcoming paper.

----
See also https://doi.org/10.5281/zenodo.10631693 

The Zenodo repository includes the parameter files and outputs (in NetCDF files). 

----
postprocessing_observations.ipynb

Artificial Observations (with no phases, with P phase, with P, S and PS)

view_artificial_observations.ipynb 

Figures in the paper viewing the artificial impacts 

postprocessing_simulations.ipynb

Figures in the paper viewing the artificial impacts 

build_Moon_surface_model.ipynb

3D models of surface and Moho topography 

build_Moon_scatter.ipynb

3D models of random lunar scatter


----
This was the code I used to set up my local python environment - called postprocessing.
I used Miniconda - which will need to be set up first. 


conda create -n postprocessing python=3.8  
conda activate postprocessing

conda install jupyter pyyaml matplotlib basemap pandas future  
conda install -c conda-forge obspy netcdf4  
pip install opencv-python pyvista vtk

pip install --find-links=https://irfpy.irf.se/sdist/ irfpy.util -U  
    # Moon basemap  
pip install --no-index --find-links=https://irfpy.irf.se/sdist irfpy.planets -U

conda install basemap

input_files/ImpactParameters.csv 

Timing, location and trajectory of the Artificial Impacts. 
