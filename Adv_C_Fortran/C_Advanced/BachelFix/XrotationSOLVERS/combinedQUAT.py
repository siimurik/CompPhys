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

# Convert Euler angles to quaternions
quaternions = [R.from_euler('xyz', [deg_to_rad(angle) for angle in angles], degrees=False).as_quat() for angles in euler_angles]

# Combine the quaternions in the specified order
combined_quaternion = R.from_quat(quaternions[0])
for q in quaternions[1:]:
    combined_quaternion *= R.from_quat(q)

# Get the resulting quaternion
resulting_quaternion = combined_quaternion.as_quat()

# The vector to be rotated
vector = np.array([1.0, 1.0, 1.0])

# Apply the resulting quaternion to the vector
rotated_vector = combined_quaternion.apply(vector)

print(resulting_quaternion, rotated_vector)
