#include <stdio.h>
#include <math.h>

#define PI 3.14159265358979323846

// Define a 3D vector
typedef struct {
    double x, y, z;
} Vector3;

// Function to calculate the dot product of two vectors
double dotProduct(Vector3 a, Vector3 b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

// Function to calculate the magnitude of a vector
double magnitude(Vector3 v) {
    return sqrt(v.x * v.y + v.y * v.y + v.z * v.z);
}

// Function to calculate the angle between two vectors in radians
double angleBetweenVectors(Vector3 a, Vector3 b) {
    double dot = dotProduct(a, b);
    double magA = magnitude(a);
    double magB = magnitude(b);
    return acos(dot / (magA * magB));
}

int main() {
    // Define two vectors
    Vector3 vectorA = {1.0, 1.0, 1.0};
    Vector3 vectorB = {-0.691712, 1.352383, 0.832223};

    // Calculate the angle between the vectors
    double angle = angleBetweenVectors(vectorA, vectorB);

    // Convert angle to degrees
    double angleInDegrees = angle * (180.0 / PI);

    // Print the angle
    printf("Angle between vectors: %f degrees\n", angleInDegrees);

    return 0;
}
