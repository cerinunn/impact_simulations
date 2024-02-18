import numpy as np
from object_latlon import Object
from gen_scripts_latlon import latlon_to_cartesian_radius

class Ellipsoid(Object):
    def __init__(self, model, vp, vs, rho, dim, loc=None, random_mag=0):
        """
        :param model: The instance of :class:`~model.Model` object shape is injected into.
        :type  model: :class:`~model.Model`
        :param vp:    Homogenous p-wave velocity for ellipsoid.
        :type vp:   float
        :param vs: Homogenous s-wave velocity for ellipsoid.
        :type vs:   float
        :param rho: Homogenous density for ellipsoid.
        :type rho: float
        :param dim: Dimensions of the ellipsoid. If single value then no rotation and all radii are equal (sphere). If 6 elements, these must be given in the following order: [rad_x, rad_y, rad_z, theta, phi, expand_int] where the first 3 elements are the radii in each direction, theta and phi are rotation angles away from the x and z aces and expand_int is an integer value with which to scale the grid in which the shape is searched for. See notes on expand_int below.
        :type dim: single value, or 6-element list/array
        :param loc: [x,y,z] of centre of ellipsoid. (type=Spherical, the convention is that x = radius (out of the page), y = longitude, z = latitude, 0,0,0 at the centre of the body).
        :type loc: 3-element list or numpy array 
        """

        self.shape_name = "ellipsoid"
        print('making an ellipsoide')
        super().__init__(model, vp, vs, rho, dim, loc, random_mag=random_mag)


    # def _in_shape_condition(self, rot_coords):
    #     """
    #     Checks if coordinates are within ellipsoid.
    # 
    #     :param rot_coords: Rotated coordinates to be checked
    #     :type rot_coords: Numpy array or list
    #     :return: bool
    #     """
    # 
    #     # A point (x,y,z) is inside the sphere with center (cx,cy,cz) and radius r if
    #     # 
    #     #  (x - cx)^2 + (y - cy)^2 + (z - cz)^2 < r^2 
    # 
    #     rad = (rot_coords[0] ** 2) / (self.radius[0] ** 2) + (rot_coords[1] ** 2) / (self.radius[1] ** 2) + (
    #                 rot_coords[2] ** 2) / (self.radius[2] ** 2)
    #     if rad <= 1:
    #         return True
    #     else:
    #         return False

    def _in_shape_condition_simple(self,x_centre,y_centre,z_centre):
        """
        Checks if coordinates are a sphere.
    
        :param rot_coords: Rotated coordinates to be checked
        :type rot_coords: Numpy array or list
        :return: bool
        """

        # get the grids for the main model 
        grid_radius = np.linspace(self.m.x_lim[0], self.m.x_lim[1], self.m.nx)
        grid_lon = np.linspace(self.m.y_lim[0], self.m.y_lim[1], self.m.ny)
        grid_lat = np.linspace(self.m.z_lim[0], self.m.z_lim[1], self.m.nz)


    
        # A point (x,y,z) is inside the sphere with center (cx,cy,cz) and radius r if
        # 
        #  (x - cx)^2 + (y - cy)^2 + (z - cz)^2 < r^2 

        shape_grid = np.full((self.m.nx,self.m.ny,self.m.nz), False)
        # random_grid = np.full((self.m.nx,self.m.ny,self.m.nz), 0.)

        # print('test case')
        # x, y, z = latlon_to_cartesian_radius(lat=lat1, long=lon1, radius=radius1,
        #                                                       e2=0, a=self.m.a)
        # x_centre = 0 
        # y_centre = 0 
        # z_centre = 0
        # x = 
        # test_sphere = (x-x_centre)** 2 + (y-y_centre)** 2 + (z-z_centre)** 2
        # print(test_sphere)
        # from numpy import random
        for i, radius1 in enumerate(grid_radius):
            for j, lon1 in enumerate(grid_lon): 
                for k, lat1 in enumerate(grid_lat): 
                    # Get coordinates of grid
                    x, y, z = latlon_to_cartesian_radius(lat=lat1, long=lon1, radius=radius1,
                                                                          e2=0, a=self.m.a)
                    test_sphere = (x-x_centre)** 2 + (y-y_centre)** 2 + (z-z_centre)** 2
                    if test_sphere < self.dim[0]**2:
                        shape_grid[i,j,k] = True
                        # random_grid[i,j,k] = random.random() -0.5

                        
        # self.vp =  vs, rho

        print('shape_grid ',shape_grid.min(),shape_grid.max())

        self.m.bm_rho = np.where(shape_grid == True, self.rho, self.m.bm_rho)
        self.m.bm_vp = np.where(shape_grid == True, self.vp, self.m.bm_vp)
        self.m.bm_vs = np.where(shape_grid == True, self.vs, self.m.bm_vs)

        # self.m.bm_rho = np.where(random_grid != 0., random_grid, self.m.bm_rho)
        # self.m.bm_vp = np.where(random_grid != 0., random_grid, self.m.bm_vp)
        # self.m.bm_vs = np.where(random_grid != 0., random_grid, self.m.bm_vs)
        # print('Random generation within the array')


        print('self.m.bm_rho ',self.m.bm_rho.min(),self.m.bm_rho.max())

        print('Count ', np.count_nonzero((np.where(self.m.bm_rho > 2.9, True, False))))

        # copy the array, reversing the depth parameter to use radius
        # v_rho[::-1,:, : ] = self.bm_rho
        # v_vp[::-1,:, : ] = self.bm_vp
        # v_vs[::-1,:, : ] = self.bm_vs

        # test_sphere = (x_coords[0] ** 2) / (self.radius[0] ** 2) + (coords[1] ** 2) / (self.radius[1] ** 2) + (
        #             coords[2] ** 2) / (self.radius[2] ** 2)
        # if test_sphere <= self.radius[0] ** 2:
        #     return True
        # else:
        #     return False


    def set_dimensions(self, dimensions):
        """
        Set dimensions for ellipsoid.

        :param dimensions: Either single value or 6-element array/list. See constructor for details.
        :type dimensions:  float/int or 6-element array/list
        """
        if type(dimensions) == float or type(dimensions) == int or len(dimensions) == 1:
            self.radius = [dimensions, dimensions, dimensions]
            self.theta = 0
            self.phi = 0
            self.expand_int = int(1)

        elif len(np.array(dimensions))==6:
            self.dim = dimensions

            self.radius = self.dim[:3]
            self.theta = self.dim[3]
            self.phi = self.dim[4]
            self.expand_int = int(self.dim[5])
        else:
            raise ValueError("Dim/radius must have either 1 entry (sphere radius) or 5 (3 radii + theta, phi)")

        self._gen_obj()
        # self._reset_sa_centre()


