import numpy as np

# Convert angles (elevation and azimuth) to a direction vector
def angle_to_dir(EL, AZ):
    d = np.pi / 180  # Conversion factor from degrees to radians
    x = np.cos(AZ * d) * np.cos(EL * d)
    y = np.sin(AZ * d) * np.cos(EL * d)
    z = np.sin(EL * d)
    return np.array([x, y, z])

# Horizon base data for Venus, Earth, and Sun
AZ_venus_HOR, EL_venus_HOR = 75.806271, 69.278559
AZ_earth_HOR, EL_earth_HOR = 94.647066, 66.507475
AZ_sun_HOR,   EL_sun_HOR   = 89.324328, 13.697882

# Compute horizon direction vectors and their average
HOR_venus = angle_to_dir(EL_venus_HOR, AZ_venus_HOR)
HOR_earth = angle_to_dir(EL_earth_HOR, AZ_earth_HOR)
HOR_sun   = angle_to_dir(EL_sun_HOR,   AZ_sun_HOR)
HOR_avrg = (HOR_venus + HOR_earth + HOR_sun) / 3
HOR_final_avrg = HOR_avrg / np.linalg.norm(HOR_avrg)
print("Average Horizons azimuth vector:", HOR_final_avrg)

# Initial coordinates and angles
xk, yk, zk = 0.11096966105001496, -0.63620006179630684, 0.76350194216964518
hx, hy, hz = HOR_final_avrg
X0, Y0, Z0 = 1.193314605098880, 4.203848106904280, -0.534337414313240

# Compute rotation angles
alpha = -np.arctan(yk / xk)
beta = np.arcsin(zk)
gamma = np.radians(-11.032961419914068)
delta = -np.arcsin(hz)
eta = np.pi / 2
iota = -np.arctan(hy / hx) + eta
print("\nDegrees:\nalpha = ", np.degrees(alpha), '\nbeta =', np.degrees(beta), '\ngamma = ', np.degrees(gamma),
      '\ndelta = ', np.degrees(delta), '\neta = ', np.degrees(eta), '\niota = ', np.degrees(iota))


def sequentialRotations(x, y, z):
    """
    Applies a series of rotations to a 3D vector using specified angles around the XYZ axes.

    -------------------------------------------------------
    Function Description:
    -------------------------------------------------------
    This function, also known as the "six-step" algorithm, performs a series of 
    rotations on a 3D vector based on input angles around the X, Y, and Z axes.
    The rotations are applied in the following order:
    1. Rotation around the Z-axis (alpha)
    2. Rotation around the Y-axis (beta)
    3. Rotation around the X-axis (gamma)
    4. Rotation around the Y-axis (delta)
    5. Rotation around the Z-axis [+ 6. Rotation around the Z-axis by 90 degrees] (iota)
    
    These rotations are used to transform the input vector (x, y, z) through a sequence
    of coordinate transformations to achieve the desired orientation.

    -------------------------------------------------------
    Parameters:
    -------------------------------------------------------
    x : float
        The x-coordinate of the initial vector.
    y : float
        The y-coordinate of the initial vector.
    z : float
        The z-coordinate of the initial vector.

    -------------------------------------------------------
    Returns:
    -------------------------------------------------------
    numpy.ndarray
        The transformed vector after applying the specified rotations, flattened to a 1D array.

    -------------------------------------------------------
    Notes:
    -------------------------------------------------------
    - The function rotates the vector around the Z, Y, and X axes in the specified sequence.

    -------------------------------------------------------
    Example:
    -------------------------------------------------------
    >>> import numpy as np
    >>> sequentialRotations(1, 0, 0)
    [ 0.27924625  0.95326278 -0.1153759 ]
    """
    # Rotation around z-axis (alpha)
    a1 = np.array([[np.cos(alpha), -np.sin(alpha), 0], [np.sin(alpha), np.cos(alpha), 0], [0, 0, 1]])
    a = a1 @ np.array([[x, y, z]]).T
    # Rotation around y-axis (beta)
    b2 = np.array([[np.cos(beta), 0, np.sin(beta)], [0, 1, 0], [-np.sin(beta), 0, np.cos(beta)]])
    b = b2 @ a
    # Rotation around x-axis (gamma)
    g1 = np.array([[1, 0, 0], [0, np.cos(gamma), -np.sin(gamma)], [0, np.sin(gamma), np.cos(gamma)]])
    g = g1 @ b
    # Rotation around y-axis (delta)
    d2 = np.array([[np.cos(delta), 0, np.sin(delta)], [0, 1, 0], [-np.sin(delta), 0, np.cos(delta)]])
    d = d2 @ g
    # Rotation around z-axis (iota - eta)
    i1 = np.array([[np.cos(iota), -np.sin(iota), 0], [np.sin(iota), np.cos(iota), 0], [0, 0, 1]])
    return (i1 @ d).flatten()

print(sequentialRotations(1,1,1))