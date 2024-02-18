#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
"""

import math
import numpy as np
import netCDF4 as nc


class Model(object):

    def __init__(self, type,  x_lim, y_lim, z_lim, elements_per_wavelength,
      dominant_freq, min_velocity, oversaturation=1, a=None, radius_lim=None, 
      inhomogeneity_lim=[0.,0.],
      linear_inhomogeneity_alpha=0.5, linear_inhomogeneity_depth=None):
        """
        Class which creates a number of 3D arrays that hold the Vp, Vs and density values in x,y,z space. The x,y,z, space can be the same size, or a subset of your simulation domain. Note that it can also be larger than the domain specified in the AxiSEM-3D simulation but anything outside of the AxiSEM-3D domain is not incorporated and AxiSEM-3D will NOT produce an error to warn you.
        :param type: 'Cartesian' or 'Spherical'
        :type type: "str"
        :param x_lim: 1D numpy array with 2 elements: [x_min, x_max] for the domain of interest.
        :type x_lim: 1D Numpy array.
        :param y_lim: 1D numpy array with 2 elements: [y_min, y_max]
        :type y_lim: 1D Numpy array.
        :param z_lim: 1D numpy array with 2 elements: [z_min, z_max]
        :type z_lim: 1D Numpy array.
        :param elements_per_wavelength: Elements per wavelength in the 1D mesh for AxiSEM3D. Controls resolution of 3D models. See notes below.
        :type elements_per_wavelength: int/float
        :param dominant_freq: Dominant frequency of 1D mesh. Controls resolution of 3D models. See notes below.
        :type dominant_freq: float
        :param min_velocity:  Minimum velocity for mesh. Controls resolution of 3D models. See notes below.
        :type min_velocity: float
        :param oversaturation: Scaling for mesh resolution. E.g. a value of 2 doubles the array resolution relative to the 1D mesh. See notes below.
        :type oversaturation: int/float or array of 3 ints for x/y/z lat/lon/depth
        :returns: Model object.
        If type=Spherical, the convention is that x = radius (out of the page), y = longitude, z = latitude, 0,0,0 at the centre of the body
        # radius (x) should be defined from center to surface.

        """

        print('initiating model_latlon')

        self.type = type.upper()
        self.x_lim = x_lim
        self.y_lim = y_lim
        self.z_lim = z_lim
        self.epw = elements_per_wavelength
        self.freq = dominant_freq
        self.min_velocity = min_velocity


        # Radius of spherical domains must be set:
        if self.type == "SPHERICAL":
            if a == None:
                raise ValueError("ERROR: Sphere radius (a) must be defined so that x_lim (depths) have context.")
            else:
                self.a = a
        
        # Calculate min wavelength from frequency and min velocity
        self.min_wavelength = min_velocity/dominant_freq
        print('self.min_wavelength',self.min_wavelength)

        # Set number of elements:
        self._set_numel(oversaturation)

        # Some details on the original 3D model (homogenous in this case)
        self.padding = np.array([0, 0, 0])
        self.unpadded_n = np.array([self.nx, self.ny, self.nz])
        self.original_shape = np.array([self.nx, self.ny, self.nz])

        # Calculate the length of each dimension:
        self.x_length = self.x_lim[1] - self.x_lim[0]
        self.y_length = self.y_lim[1] - self.y_lim[0]
        self.z_length = self.z_lim[1] - self.z_lim[0]

        # Calculate spatial step sizes
        self.dx = self.x_length / self.nx
        self.dy = self.y_length / self.ny
        self.dz = self.z_length / self.nz

        # Define a default background model arrays for Rho, Vp and Vs which are zeros 3D arrays of the correct size:
        # Note here that the model is defined originally as homogenous
        # This could be edited to take an input from the user of some pre-defined 3D array that has non-zero values
        self.bm_rho = np.zeros((self.nx, self.ny, self.nz))
        self.bm_vp  = np.zeros((self.nx, self.ny, self.nz))
        self.bm_vs  = np.zeros((self.nx, self.ny, self.nz))

        if inhomogeneity_lim != [0,0]:
            self.inhomogeneity = np.random.uniform(inhomogeneity_lim[0], inhomogeneity_lim[1], size=(self.nx, self.ny, self.nz))
            self.bm_rho = self.inhomogeneity 
            self.bm_vp = self.inhomogeneity 
            self.bm_vs = self.inhomogeneity 

        elif linear_inhomogeneity_depth is not None:
            if self.x_lim[0] > self.x_lim[1] : 
                raise ValueError("ERROR: x_lim[1] (radius) must be greater than x_lim[0]")

            print(self.x_lim[0], self.x_lim[1], self.x_lim[1]-self.x_lim[0])
            depth_list = []
            grid_radius = np.linspace(self.x_lim[0], self.x_lim[1], self.nx)        
            for x in grid_radius:
                D = self.a - x
                # if D < 0:
                #     break
                depth_list.append(D)

            # self.inhomogeneity = np.zeros((len(alpha_list), self.ny, self.nz))
            # for a1, alpha in enumerate(alpha_list):
            #     for y in range(self.ny):
            #         for z in range (self.nz):
            #             self.inhomogeneity[a1,y,z] = alpha * np.random.uniform(-1., 1., size=None)

            # random distribution from -1 to 1 
            self.inhomogeneity = np.random.uniform(-1., 1., size=(len(depth_list), self.ny, self.nz))
            # linear version of the random matrix - multiply by alpha at the top and 
            # zero at the bottom

            alpha_list = []
            for d, D in enumerate(depth_list):
                alpha = linear_inhomogeneity_alpha * (((D/1000)/-linear_inhomogeneity_depth) + 1)
                alpha_list.append(alpha)
                self.inhomogeneity[d,:,:] = self.inhomogeneity[d,:,:] * alpha
            print(alpha_list)
            print('depth list ', depth_list)
            print('grid radius ', grid_radius)

            self.bm_rho = self.inhomogeneity 
            self.bm_vp = self.inhomogeneity 
            self.bm_vs = self.inhomogeneity 

            print("")
            print('rho')
            print(self.nx, self.ny, self.nz)
            print('Spacing ', self.dx, self.dy, self.dz)
            print(self.inhomogeneity.size)
            print('rho - top layer min max')
            print(self.bm_rho[-1,:,:].min(),self.bm_rho[-1,:,:].max())
            print('rho - bottom layer min max')
            print(self.bm_rho[0,:,:].min(),self.bm_rho[0,:,:].max())

            # for x in 
            # for x, y, z in ):
            #     
            # 

        print('shape ', self.bm_rho.shape)
        print('initiated')




    def writeNetCDF(self, filename, paraview=False, add_longitude=0):
        """
        Writes 3D arrays for velocity and density to a .nc file for inclusion in AxiSEM-3D simulation.

        :param filename: Output filename. Note that the suffix '.nc' must be included.
        :type filename: Str
        :return: Outputs .nc file.
        """
        if self.type == "CARTESIAN":
            print("Writing cartesian model...")
            f = nc.Dataset(filename, 'w', format='NETCDF4')

            # Create dimension arrays
            x_array = np.linspace(-self.x_lim[1], self.x_lim[1], self.nx)
            y_array = np.linspace(-self.y_lim[1], self.y_lim[1], self.ny)
            z_array = np.linspace(-self.z_lim[0], self.z_lim[1], self.nz)

            # Create the dimensions:
            x_dim = f.createDimension('x_dim', self.nx)
            y_dim = f.createDimension('y_dim', self.ny)
            z_dim = f.createDimension('z_dim', self.nz)

            # Creating the variables:
            x   = f.createVariable('x', 'f4', ('x_dim',))
            y   = f.createVariable('y', 'f4', ('y_dim',))
            z   = f.createVariable('z', 'f4', ('z_dim',))
            v_rho = f.createVariable('rho', 'f4', ('x_dim', 'y_dim', 'z_dim',))
            v_vp  = f.createVariable('vp', 'f4', ('x_dim', 'y_dim', 'z_dim',))
            v_vs  = f.createVariable('vs', 'f4', ('x_dim', 'y_dim', 'z_dim',))

            # Assigning values to the variables:
            x[:] = x_array
            y[:] = y_array
            z[:] = z_array
            v_rho[:,:,:] = self.bm_rho
            v_vp[:,:,:]  =  self.bm_vp
            v_vs[:,:,:]  =  self.bm_vs
            print('Data written to file ', filename)
            f.close()
        elif self.type == "SPHERICAL":
            print("Writing spherical model...")

            grid_radius = np.linspace(self.x_lim[0], self.x_lim[1], self.nx)
            grid_lon = np.linspace(self.y_lim[0], self.y_lim[1], self.ny)
            # grid_depth = np.linspace(self.z_lim[0], self.z_lim[1], self.nz)
            grid_lat = np.linspace(self.z_lim[0], self.z_lim[1], self.nz)

            f = nc.Dataset(filename, 'w', format='NETCDF4')
            # Create dimension arrays
            # We now create the dimensions:
            radius = f.createDimension('radius', self.nx)
            lon = f.createDimension('lon', self.ny)
            lat = f.createDimension('lat', self.nz)
            # depth = f.createDimension('depth', self.nz)


            # Creating the variables:
            lats = f.createVariable('lat', 'f4', ('lat',))
            lats.units = 'degrees_north'
            lats.long_name = 'latitude'

            lons = f.createVariable('lon', 'f4', ('lon',))
            lons.units = 'degrees_east'
            lons.long_name = 'longitude'

            # z = f.createVariable('depth', 'f4', ('depth',))
            # z.units = 'meters'

            z = f.createVariable('radius', 'f4', ('radius',))
            z.units = 'meters'
            z.positive = 'up'

            # v_rho = f.createVariable('rho', 'f4', ('lat', 'lon', 'depth',))
            # v_vp = f.createVariable('vp', 'f4', ('lat', 'lon', 'depth',))
            # v_vs = f.createVariable('vs', 'f4', ('lat', 'lon', 'depth',))

            # swapped order around, 
            v_rho = f.createVariable('rho', 'f4', ('radius', 'lon', 'lat',))
            v_vp = f.createVariable('vp', 'f4', ('radius', 'lon', 'lat',))
            v_vs = f.createVariable('vs', 'f4', ('radius', 'lon', 'lat',))

            # Assigning values to the variables:
            lats[:] = grid_lat
            lons[:] = grid_lon + add_longitude
            if paraview == True:
                # print("PARAVIEW FLAG = TRUE: USE MODEL FOR VISUALS BUT NOT SIMULATION")
                # z[:] = self.a - grid_depth
                print('paraview == True not implemented')
                return
            elif paraview == False:
                print("PARAVIEW FLAG = FALSE: MODEL OKAY FOR SIMULATION")
                z[:] = grid_radius
            else:
                raise ValueError("paraview flag must be True or False")

            # v_rho[:, :, :] = self.bm_rho
            # v_vp[:, :, :] = self.bm_vp
            # v_vs[:, :, :] = self.bm_vs
            # copy the array (radius, lon, lat)
            v_rho[:, :, :] = self.bm_rho
            v_vp[:, :, :] = self.bm_vp
            v_vs[:, :, :] = self.bm_vs
            print('Data written to file ', filename)

            if self.type == "SPHERICAL" and paraview == False:
                self._print_inparam_model_script(filename)
            f.close()








    # Functions for updating 
    def set_bm_rho(self, bm_rho):
        """
        Update the background model for density

        :param bm_rho: 3D array to replace current array
        :type bm_rho: 3D numpy array
        :returns: None.
        """
        self.bm_rho = bm_rho

    def set_bm_vp(self, bm_vp):
        """
        Update the background model for P-wave velocity

        :param bm_vp: 3D array to replace current array
        :type bm_vp: 3D numpy array
        :returns: None.
        """
        self.bm_rho = bm_vp

    def set_bm_vs(self, bm_vs):
        """
        Update the background model for S-wave velocity

        :param bm_vs: 3D array to replace current array
        :type bm_vs: 3D numpy array
        :returns: None.
        """
        self.bm_rho = bm_vs


    def _set_numel(self, oversaturation):

        # sort oversaturation - can be an integer or array of length three
        if type(oversaturation) == int or type(oversaturation) == float:
            oversat_x = oversaturation
            oversat_y = oversaturation
            oversat_z = oversaturation

        elif np.size(np.array(oversaturation))==3:
            oversat_x = oversaturation[0]
            oversat_y = oversaturation[1]
            oversat_z = oversaturation[2]


        if self.type.upper() =="SPHERICAL":
            # x,y,z = latlon_to_cartesian(lat=self.x_lim, long=self.y_lim, depth=self.z_lim, e2=0, a=self.a)
            # self.nz = int(np.abs(oversat_z * math.ceil((x[1] - x[0]) * self.epw / self.min_wavelength)))
            # self.ny = int(np.abs(oversat_y * math.ceil((y[1] - y[0]) * self.epw / self.min_wavelength)))
            # self.nx = int(np.abs(oversat_x * math.ceil((z[1] - z[0]) * self.epw / self.min_wavelength)))

            constant = 30000.
            
            # self.nx = int(np.abs(oversat_x * math.ceil((self.x_lim[1] - self.x_lim[0]) * self.epw / self.min_wavelength)))
            # self.ny = int(np.abs(oversat_y * math.ceil((self.y_lim[1] - self.y_lim[0] * constant) * self.epw / self.min_wavelength)))
            # self.nz = int(np.abs(oversat_z * math.ceil((self.z_lim[1] - self.z_lim[0] * constant) * self.epw / self.min_wavelength)))

            self.nx = int(np.abs(oversat_x * (self.x_lim[1] - self.x_lim[0]) * self.epw / self.min_wavelength))
            self.ny = int(np.abs(oversat_y * (self.y_lim[1] - self.y_lim[0]) * constant * self.epw / self.min_wavelength))
            self.nz = int(np.abs(oversat_z * (self.z_lim[1] - self.z_lim[0]) * constant * self.epw / self.min_wavelength))


            print('x, y, z ', self.nx, self.ny, self.nz)

            # Ensure that nx, ny, nz are odd
            if self.nx % 2 == 0:
                self.nx += 1
            if self.ny % 2 == 0:
                self.ny += 1
            if self.nz % 2 == 0:
                self.nz += 1

            print('x, y, z ', self.nx, self.ny, self.nz)


        elif self.type.upper() =="CARTESIAN":
            self.nx = oversaturation*math.ceil((self.x_lim[1] - self.x_lim[0])*self.epw /self.min_wavelength)
            self.ny = oversaturation*math.ceil((self.y_lim[1] - self.y_lim[0])*self.epw /self.min_wavelength)
            self.nz = oversaturation*math.ceil((self.z_lim[1] - self.z_lim[0])*self.epw /self.min_wavelength)
        else:
            raise ValueError("model variable 'type' must be 'cartesian' or 'spherical'. ")



    def _print_inparam_model_script(self, fname):
        print("# ___________________________________________________________________________________")
        print(f"    - {fname[:-3]}:")
        print(f"        activated: true")
        print(f"        class_name: StructuredGridV3D")
        print(f"        nc_data_file: {fname}")
        print(f"        coordinates:")
        print(f"            horizontal: LATITUDE_LONGITUDE")
        print(f"            vertical: RADIUS")
        print(f"            ellipticity: FILL THIS IN - true/false")
        print(f"            depth_below_solid_surface: FILL THIS IN - true/false")
        print(f"            nc_variables: [radius, lon, lat]")
        print(f"            data_rank: [2, 1, 0]")
        print(f"            length_unit: m")
        print(f"            angle_unit: degree")
        print(f"            undulated_geometry: false")
        print(f"            whole_element_inplane: false")
        print(f"        properties:")
        print(f"            - VP:")
        print(f"                nc_var: vp")
        print(f"                factor: 1")
        print(f"                reference_kind: FILL THIS IN - ABS/REF1D/REF3D")
        print(f"            - VS:")
        print(f"                nc_var: vs")
        print(f"                factor: 1")
        print(f"                reference_kind: FILL THIS IN - ABS/REF1D/REF3D")
        print(f"            - RHO:")
        print(f"                nc_var: rho")
        print(f"                factor: 1")
        print(f"                reference_kind: FILL THIS IN - ABS/REF1D/REF3D")
        print(f"        store_grid_only_on_leaders: true")
        print("# ___________________________________________________________________________________")


if __name__ == '__main__':
    linear_inhomogeneity_depth = 20.


    MOON_RADIUS_IN_KM = 1737.1

    # We dont need to make a model that spans the whole domain, just the part we are interested in injecting a plume in: 
    radius = MOON_RADIUS_IN_KM*1000
    lat_lim = [0, 1]
    long_lim = [0, 1]


    radius_lim = [(MOON_RADIUS_IN_KM-linear_inhomogeneity_depth)*1000, MOON_RADIUS_IN_KM*1000]

    # If type=Spherical, the convention is that x = radius (out of the page), y = longitude, z = latitude, 0,0,0 at the centre of the body


    # Create our global model, including inhomogeneity
    glob_m = Model("spherical", x_lim=radius_lim, y_lim=long_lim, z_lim=lat_lim, elements_per_wavelength=1, 
                   dominant_freq=5.234333, min_velocity=4000, oversaturation=1, a=radius, 
                   linear_inhomogeneity_depth=linear_inhomogeneity_depth)
