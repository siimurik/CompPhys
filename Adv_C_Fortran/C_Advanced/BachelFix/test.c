/*
Compile and execute with:
    $ gcc test.c -o test -lm
    $ ./test
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define PI 3.14159265358979323846

typedef struct {
    double w, x, y, z;
} Quaternion;

typedef struct {
    double x, y, z;
} Vector3;

// Vector operations
Vector3 vec_add(Vector3 a, Vector3 b) {
    return (Vector3){a.x + b.x, a.y + b.y, a.z + b.z};
}

Vector3 vec_sub(Vector3 a, Vector3 b) {
    return (Vector3){a.x - b.x, a.y - b.y, a.z - b.z};
}

Vector3 vec_scale(Vector3 v, double s) {
    return (Vector3){v.x * s, v.y * s, v.z * s};
}

double vec_norm(Vector3 v) {
    return sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
}

Vector3 vec_normalize(Vector3 v) {
    double norm = vec_norm(v);
    return (Vector3){v.x / norm, v.y / norm, v.z / norm};
}

Vector3 az_el_to_unitvec(double az, double el) {
    double az_rad = az * PI / 180.0;
    double el_rad = el * PI / 180.0;
    return (Vector3){
        cos(az_rad) * cos(el_rad),
        -sin(az_rad) * cos(el_rad),
        sin(el_rad)
    };
}

// Quaternion operations
Quaternion quat_from_angle_axis(double angle, Vector3 axis) {
    double half_angle = angle / 2.0;
    double s = sin(half_angle);
    return (Quaternion){cos(half_angle), axis.x * s, axis.y * s, axis.z * s};
}

Quaternion quat_multiply(Quaternion q1, Quaternion q2) {
    return (Quaternion){
        q1.w * q2.w - q1.x * q2.x - q1.y * q2.y - q1.z * q2.z,
        q1.w * q2.x + q1.x * q2.w + q1.y * q2.z - q1.z * q2.y,
        q1.w * q2.y - q1.x * q2.z + q1.y * q2.w + q1.z * q2.x,
        q1.w * q2.z + q1.x * q2.y - q1.y * q2.x + q1.z * q2.w
    };
}

Vector3 quat_rotate_vector(Quaternion q, Vector3 v) {
    Quaternion qv = {0, v.x, v.y, v.z};
    Quaternion q_conj = {q.w, -q.x, -q.y, -q.z};
    Quaternion q_result = quat_multiply(quat_multiply(q, qv), q_conj);
    return (Vector3){q_result.x, q_result.y, q_result.z};
}

// Coordinate transformation using quaternions
Vector3 coord_turn(Vector3 v) {
    double radParam =  PI / 180.0;
    double alpha    =  80.1056830575758170 * radParam;
    double beta     =  49.7739016604574402 * radParam;
    double gamma    = -11.0329389306380126 * radParam;  // value found using rollCALC.c
    double delta    = -50.6492951614895688 * radParam;
    double iota     =   2.2052577138069775 * radParam;

    Quaternion qz1 = quat_from_angle_axis(alpha, (Vector3){0, 0, 1});
    Quaternion qy  = quat_from_angle_axis(beta,  (Vector3){0, 1, 0});
    Quaternion qx  = quat_from_angle_axis(gamma, (Vector3){1, 0, 0});
    Quaternion qy2 = quat_from_angle_axis(delta, (Vector3){0, 1, 0});
    Quaternion qz2 = quat_from_angle_axis(iota,  (Vector3){0, 0, 1});

    Quaternion q_combined = quat_multiply(qz2, quat_multiply(qy2, quat_multiply(qx, quat_multiply(qy, qz1))));

    return quat_rotate_vector(q_combined, v);
}

void print_vector(const char* label, Vector3 vec) {
    printf("%s: [%.6f, %.6f, %.6f]\n", label, vec.x, vec.y, vec.z);
}

int main(){

    Vector3 v1 = {1.0, 1.0, 1.0};
    Vector3 v2 = coord_turn(v1);

    print_vector("v2", v2);

    return 0;
}
