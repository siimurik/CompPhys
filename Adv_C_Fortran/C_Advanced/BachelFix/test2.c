#include <stdio.h>
#include <math.h>

#define PI 3.14159265358979323846

// Converts degrees to radians
double toRadians(double degrees) {
    return degrees * (PI / 180.0);
}

// Define a quaternion
typedef struct {
    double w, x, y, z;
} Quaternion;

// Define a 3D vector
typedef struct {
    double x, y, z;
} Vector3;

// Function to create a quaternion from an axis-angle representation
Quaternion createFromAxisAngle(double angle, double x, double y, double z) {
    Quaternion q;
    double halfAngle = angle / 2.0;
    double sinHalfAngle = sin(halfAngle);
    q.w = cos(halfAngle);
    q.x = x * sinHalfAngle;
    q.y = y * sinHalfAngle;
    q.z = z * sinHalfAngle;
    return q;
}

// Function to multiply two quaternions
Quaternion multiplyQuaternions(Quaternion q1, Quaternion q2) {
    Quaternion q;
    q.w = q1.w * q2.w - q1.x * q2.x - q1.y * q2.y - q1.z * q2.z;
    q.x = q1.w * q2.x + q1.x * q2.w + q1.y * q2.z - q1.z * q2.y;
    q.y = q1.w * q2.y - q1.x * q2.z + q1.y * q2.w + q1.z * q2.x;
    q.z = q1.w * q2.z + q1.x * q2.y - q1.y * q2.x + q1.z * q2.w;
    return q;
}

// Function to rotate a vector by a quaternion
Vector3 rotateVectorByQuaternion(Vector3 v, Quaternion q) {
    Quaternion qConjugate = { q.w, -q.x, -q.y, -q.z };
    Quaternion qVector = { 0, v.x, v.y, v.z };
    Quaternion qResult = multiplyQuaternions(multiplyQuaternions(q, qVector), qConjugate);

    Vector3 rotatedVector = { qResult.x, qResult.y, qResult.z };
    return rotatedVector;
}

int main() {
    // Given Euler angles in degrees
    double yaw1 = toRadians(80.1056830575758170);
    double pitch1 = toRadians(49.7739016604574402);
    double roll = toRadians(-11.0329389306380126);
    double pitch2 = toRadians(-50.6492951614895688);
    double yaw2 = toRadians(2.2052577138069775);

    // Original vector
    Vector3 vector = {1.0, 1.0, 1.0};

    // Create quaternions for each rotation
    Quaternion qYaw1 = createFromAxisAngle(yaw1, 0, 0, 1);
    Quaternion qPitch1 = createFromAxisAngle(pitch1, 0, 1, 0);
    Quaternion qRoll = createFromAxisAngle(roll, 1, 0, 0);
    Quaternion qPitch2 = createFromAxisAngle(pitch2, 0, 1, 0);
    Quaternion qYaw2 = createFromAxisAngle(yaw2, 0, 0, 1);

    // Combine all rotations into one quaternion
    Quaternion qCombined = multiplyQuaternions(qYaw2, multiplyQuaternions(qPitch2, multiplyQuaternions(qRoll, multiplyQuaternions(qPitch1, qYaw1))));

    // Rotate the vector using the combined quaternion
    Vector3 rotatedVector = rotateVectorByQuaternion(vector, qCombined);

    // Print the rotated vector
    printf("Rotated vector: [%f, %f, %f]\n", rotatedVector.x, rotatedVector.y, rotatedVector.z);

    return 0;
}
