import numpy as np
from scipy.spatial.transform import Rotation as R

# Convert degrees to radians
def deg_to_rad(degrees):
    return np.deg2rad(degrees)

# Euler angles in degrees
euler_angles = [
    (0, 0, 80.1056830575758170),
    (0, 49.7739016604574402, 0),
    (-11.0329389306380126, 0, 0),
    (0, -50.6492951614895688, 0),
    (0, 0, 2.2052577138069775)
]

# The vector to be rotated
vector = np.array([1.0, 1.0, 1.0])

# Initialize the resulting vector
resulting_vector = vector.copy()

# Apply each rotation sequentially to the vector
for angles in euler_angles:
    # Convert Euler angles to quaternion
    quaternion = R.from_euler('xyz', [deg_to_rad(angle) for angle in angles], degrees=False).as_quat()
    
    # Create a rotation object from the quaternion
    rotation = R.from_quat(quaternion)
    
    # Apply the rotation to the vector
    resulting_vector = rotation.apply(resulting_vector)

print(resulting_vector)
