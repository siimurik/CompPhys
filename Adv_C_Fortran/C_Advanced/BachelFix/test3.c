#include <stdio.h>
#include <math.h>

#define PI 3.14159265358979323846

// Define a 3D vector
typedef struct {
    double x, y, z;
} Vector3;

// Define a quaternion
typedef struct {
    double w, x, y, z;
} Quaternion;

// Function to convert degrees to radians
double toRadians(double degrees) {
    return degrees * (PI / 180.0);
}

// Function to calculate the cross product of two vectors
Vector3 crossProduct(Vector3 a, Vector3 b) {
    Vector3 cross;
    cross.x = a.y * b.z - a.z * b.y;
    cross.y = a.z * b.x - a.x * b.z;
    cross.z = a.x * b.y - a.y * b.x;
    return cross;
}

// Function to calculate the dot product of two vectors
double dotProduct(Vector3 a, Vector3 b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

// Function to calculate the magnitude of a vector
double magnitude(Vector3 v) {
    return sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
}

// Function to normalize a vector
Vector3 normalize(Vector3 v) {
    double mag = magnitude(v);
    Vector3 normalized = { v.x / mag, v.y / mag, v.z / mag };
    return normalized;
}

// Function to create a quaternion from an angle and axis
Quaternion createQuaternion(double angle, Vector3 axis) {
    Quaternion q;
    double halfAngle = angle * 0.5;
    double s = sin(halfAngle);
    q.w = cos(halfAngle);
    q.x = axis.x * s;
    q.y = axis.y * s;
    q.z = axis.z * s;
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
    // Original and target vectors
    Vector3 vectorA = {1.0, 1.0, 1.0};
    Vector3 vectorB = {-0.691712, 1.352383, 0.832223};

    // Find the axis of rotation (cross product of vectorA and vectorB)
    Vector3 axis = crossProduct(vectorA, vectorB);

    // Normalize the axis of rotation
    axis = normalize(axis);

    // Given angle between vectors in degrees
    double angleDegrees = 46.812002;

    // Convert angle to radians
    double angleRadians = toRadians(angleDegrees);

    // Create a quaternion from the angle and axis
    Quaternion q = createQuaternion(angleRadians, axis);

    // Rotate the original vector by the quaternion
    Vector3 rotatedVector = rotateVectorByQuaternion(vectorA, q);

    // Print the rotated vector
    printf("Rotated vector: [%f, %f, %f]\n", rotatedVector.x, rotatedVector.y, rotatedVector.z);

    return 0;
}
