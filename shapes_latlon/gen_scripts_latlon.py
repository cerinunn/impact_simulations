import numpy as np

def latlon_to_cartesian(lat, long, depth, e2, a):
    # Convert the ellipsoidal (geographic) coordinates to x,y,z:
    latr = np.radians(lat)
    lonr = np.radians(long)
    depth = np.array(depth)
    
    #N = a / np.sqrt(1 - e2 * (np.sin(latr) ** 2))
    N = a
    X =  np.cos(latr) * np.cos(lonr) * (N - depth)
    Y =  np.cos(latr) * np.sin(lonr) *(N - depth)
    Z =  np.sin(latr) * ((1 - e2) * N - depth)
    return X, Y, Z#, N



def latlon_to_cartesian_radius(lat, long, radius, e2, a):
    # Convert the ellipsoidal (geographic) coordinates to x,y,z:
    latr = np.radians(lat)
    lonr = np.radians(long)
    radius = np.array(radius)

    #N = a / np.sqrt(1 - e2 * (np.sin(latr) ** 2))
    # N = a
    X =  np.cos(latr) * np.cos(lonr) * radius
    Y =  np.cos(latr) * np.sin(lonr) * radius
    Z =  np.sin(latr) * radius
    return X, Y, Z#, N

