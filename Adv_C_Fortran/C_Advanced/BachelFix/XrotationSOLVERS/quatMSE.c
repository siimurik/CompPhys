#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define PI 3.14159265358979323846

// Quaternion utilities
typedef struct {
    double w, x, y, z;
} Quaternion;

// Vector utilities
typedef struct {
    double x, y, z;
} Vector3;

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

// Coordinate transformation
Vector3 coord_turn(double angle, Vector3 v, double alpha, double beta, double delta, double iota) {
    double gamma = angle * PI / 180.0;

    Quaternion qz1 = quat_from_angle_axis(alpha, (Vector3){0, 0, 1});
    Quaternion qy = quat_from_angle_axis(beta, (Vector3){0, 1, 0});
    Quaternion qx = quat_from_angle_axis(gamma, (Vector3){1, 0, 0});
    Quaternion qy2 = quat_from_angle_axis(delta, (Vector3){0, 1, 0});
    Quaternion qz2 = quat_from_angle_axis(iota, (Vector3){0, 0, 1});
    
    Quaternion q_combined = quat_multiply(qz2, quat_multiply(qy2, quat_multiply(qx, quat_multiply(qy, qz1))));
    
    return quat_rotate_vector(q_combined, v);
}

// Compute errors
void compute_errors(double *angles, int num_angles, Vector3 vectors[], Vector3 targets[], double alpha, double beta, double delta, double iota, double *errors) {
    for (int i = 0; i < num_angles; i++) {
        double angle = angles[i];
        double total_error = 0.0;

        for (int j = 0; j < 4; j++) {
            Vector3 coord = coord_turn(angle, vectors[j], alpha, beta, delta, iota);
            double error = vec_norm(vec_sub(targets[j], coord));
            total_error += error * error;
        }
        
        errors[i] = total_error;
    }
}

// Convert azimuth and elevation to unit vectors
Vector3 az_el_to_unitvec(double az, double el) {
    double az_rad = az * PI / 180.0;
    double el_rad = el * PI / 180.0;
    return (Vector3){
        cos(az_rad) * cos(el_rad),
        -sin(az_rad) * cos(el_rad),
        sin(el_rad)
    };
}

int main() {
    // Camera vectors
    Vector3 v9190v = {3.82564, 4.18553, 1.46657};
    Vector3 v9196v = {3.70499, 4.15347, 1.55078};
    Vector3 cam_avrg = vec_scale(vec_add(v9190v, v9196v), 0.5);
    
    // Venus and Earth location vectors in IM system
    Vector3 veen_loc = {10.857, -30.625, 94.405};
    Vector3 maa_loc = {2.077, -25.212, 62.602};
    
    Vector3 veen_vec = vec_sub(veen_loc, cam_avrg);
    Vector3 veen_unitvec = vec_normalize(veen_vec);
    
    Vector3 maa_vec = vec_sub(maa_loc, cam_avrg);
    Vector3 maa_unitvec = vec_normalize(maa_vec);
    
    // Antenna and flag vectors
    Vector3 loc_36 = {11.234, 1.629, 3.688};
    Vector3 loc_230 = {7.739, 14.612, 0.488};
    Vector3 antenn_vec = vec_sub(loc_36, loc_230);
    Vector3 antenn_unitvec = vec_normalize(antenn_vec);
    
    Vector3 loc_17 = {7.502, 7.397, 2.302};
    Vector3 loc_229 = {5.468, 14.882, 0.452};
    Vector3 lipp_vec = vec_sub(loc_17, loc_229);
    Vector3 lipp_unitvec = vec_normalize(lipp_vec);
    
    // Sun unit vector
    Vector3 paike_vec = vec_scale(vec_add(antenn_unitvec, lipp_unitvec), 0.5);
    Vector3 paike_unitvec = vec_normalize(paike_vec);
    
    // Average vector
    Vector3 kesk_vec = vec_scale(vec_add(vec_add(maa_unitvec, veen_unitvec), paike_unitvec), 1.0/3.0);
    Vector3 kesk_unitvec = vec_normalize(kesk_vec);
    
    // HORIZONS
    Vector3 HOR_veen_unit = az_el_to_unitvec(75.800407, 69.286656);
    Vector3 HOR_maa_unit = az_el_to_unitvec(94.645414, 66.507913);
    Vector3 HOR_paike_unit = az_el_to_unitvec(89.324328, 13.697882);
    
    Vector3 HOR_kesk_vec = vec_scale(vec_add(vec_add(HOR_veen_unit, HOR_maa_unit), HOR_paike_unit), 1.0/3.0);
    Vector3 HOR_kesk_unitvec = vec_normalize(HOR_kesk_vec);
    
    // Rotation parameters
    Vector3 keskmine = {0.024404, -0.633554, 0.773313};
    Vector3 maa = {-0.032284, -0.397313, 0.917115};
    Vector3 veenus = {0.086761, -0.342886, 0.935362};
    Vector3 paike = {0.011457, -0.971490, 0.236802};
    
    Vector3 vectors[] = {
        {0.11096966105001496, -0.63620006179630684, 0.76350194216964518},
        {-0.024897, -0.433276, 0.900917},
        {0.071308, -0.349863, 0.934083},
        {0.253980, -0.939030, 0.231770}
    };
    
    Vector3 targets[] = {keskmine, maa, veenus, paike};
    
    double alpha = -atan(-0.63620006179630684 / 0.11096966105001496);
    double beta = asin(0.76350194216964518);
    double delta = -asin(0.77331312210251590);
    double iota = atan(-0.63355444896004787 / 0.024404);
    
    double start_angle = -11.032938932;
    double end_angle   = -11.032938930;
    double step_size = 1e-16;
    int num_angles = (int)((end_angle - start_angle) / step_size);
    
    double *angles = (double *)malloc(num_angles * sizeof(double));
    double *errors = (double *)malloc(num_angles * sizeof(double));
    
    for (int i = 0; i < num_angles; i++) {
        angles[i] = start_angle + i * step_size;
    }
    
    clock_t start_time = clock();
    
    compute_errors(angles, num_angles, vectors, targets, alpha, beta, delta, iota, errors);
    
    clock_t end_time = clock();
    double elapsed_time = (double)(end_time - start_time) / CLOCKS_PER_SEC;
    
    int min_index = 0;
    for (int i = 1; i < num_angles; i++) {
        if (errors[i] < errors[min_index]) {
            min_index = i;
        }
    }
    
    double min_value = errors[min_index];
    double best_angle = angles[min_index];
    
    printf("Angle Range: [%f, %f] degrees\n", start_angle, end_angle);
    printf("Required Accuracy: %e degrees\n", step_size);
    printf("Minimum Value: %e\n", min_value);
    printf("Best Angle: %22.16f\n", best_angle);
    printf("\nElapsed time: %.3e seconds aka %.3g minutes.\n", elapsed_time, elapsed_time/60.0);
    
    free(angles);
    free(errors);
    
    return 0;
}
