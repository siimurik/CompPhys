#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define PI 3.14159265358979323846
#define MAX_ITERATIONS 10000
#define INITIAL_TEMP 1000
#define COOLING_RATE 0.99
#define TOLERANCE 1e-9

typedef struct {
    double x, y, z;
} Vector3;

typedef struct {
    double m[3][3];
} Matrix3x3;

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

Vector3 mat_vec_mult(Matrix3x3 mat, Vector3 vec) {
    Vector3 result;
    result.x = mat.m[0][0] * vec.x + mat.m[0][1] * vec.y + mat.m[0][2] * vec.z;
    result.y = mat.m[1][0] * vec.x + mat.m[1][1] * vec.y + mat.m[1][2] * vec.z;
    result.z = mat.m[2][0] * vec.x + mat.m[2][1] * vec.y + mat.m[2][2] * vec.z;
    return result;
}

Matrix3x3 rotation_matrix_z(double angle) {
    double c = cos(angle), s = sin(angle);
    Matrix3x3 mat = {{
        {c, -s, 0},
        {s, c, 0},
        {0, 0, 1}
    }};
    return mat;
}

Matrix3x3 rotation_matrix_y(double angle) {
    double c = cos(angle), s = sin(angle);
    Matrix3x3 mat = {{
        {c, 0, s},
        {0, 1, 0},
        {-s, 0, c}
    }};
    return mat;
}

Matrix3x3 rotation_matrix_x(double angle) {
    double c = cos(angle), s = sin(angle);
    Matrix3x3 mat = {{
        {1, 0, 0},
        {0, c, -s},
        {0, s, c}
    }};
    return mat;
}

Vector3 coord_turn(double angle, Vector3 v, double alpha, double beta, double delta, double iota) {
    double gamma = angle * PI / 180.0;
    Matrix3x3 Rz = rotation_matrix_z(alpha);
    Matrix3x3 Ry = rotation_matrix_y(beta);
    Matrix3x3 Rx = rotation_matrix_x(gamma);
    Matrix3x3 Rdelta = rotation_matrix_y(delta);
    Matrix3x3 Riota = rotation_matrix_z(iota);
    
    Vector3 result = mat_vec_mult(Riota, mat_vec_mult(Rdelta, mat_vec_mult(Rx, mat_vec_mult(Ry, mat_vec_mult(Rz, v)))));
    return result;
}

double compute_total_error(double angle, Vector3 vectors[], Vector3 targets[], double alpha, double beta, double delta, double iota) {
    double total_error = 0.0;
    for (int j = 0; j < 4; j++) {
        Vector3 coord = coord_turn(angle, vectors[j], alpha, beta, delta, iota);
        double error = vec_norm(vec_sub(targets[j], coord));
        total_error += error * error;
    }
    return total_error;
}

double random_double(double min, double max) {
    return min + (double)rand() / RAND_MAX * (max - min);
}

double simulated_annealing(Vector3 vectors[], Vector3 targets[], double alpha, double beta, double delta, double iota, double start_angle, double end_angle) {
    double current_angle = random_double(start_angle, end_angle);
    double best_angle = current_angle;
    double current_temp = INITIAL_TEMP;
    double best_error = compute_total_error(current_angle, vectors, targets, alpha, beta, delta, iota);
    
    for (int i = 0; i < MAX_ITERATIONS; i++) {
        double new_angle = random_double(start_angle, end_angle);
        double new_error = compute_total_error(new_angle, vectors, targets, alpha, beta, delta, iota);
        
        if (new_error < best_error || exp((best_error - new_error) / current_temp) > random_double(0, 1)) {
            current_angle = new_angle;
            best_error = new_error;
            best_angle = new_angle;
        }
        
        current_temp *= COOLING_RATE;
        
        if (current_temp < TOLERANCE) break;
    }
    
    return best_angle;
}

int main() {
    Vector3 v9190v = {3.82564, 4.18553, 1.46657};
    Vector3 v9196v = {3.70499, 4.15347, 1.55078};
    Vector3 cam_avrg = vec_scale(vec_add(v9190v, v9196v), 0.5);
    
    Vector3 veen_loc = {10.857, -30.625, 94.405};
    Vector3 maa_loc = {2.077, -25.212, 62.602};
    
    Vector3 veen_vec = vec_sub(veen_loc, cam_avrg);
    Vector3 veen_unitvec = vec_normalize(veen_vec);
    
    Vector3 maa_vec = vec_sub(maa_loc, cam_avrg);
    Vector3 maa_unitvec = vec_normalize(maa_vec);
    
    Vector3 loc_36 = {11.234, 1.629, 3.688};
    Vector3 loc_230 = {7.739, 14.612, 0.488};
    Vector3 antenn_vec = vec_sub(loc_36, loc_230);
    Vector3 antenn_unitvec = vec_normalize(antenn_vec);
    
    Vector3 loc_17 = {7.502, 7.397, 2.302};
    Vector3 loc_229 = {5.468, 14.882, 0.452};
    Vector3 lipp_vec = vec_sub(loc_17, loc_229);
    Vector3 lipp_unitvec = vec_normalize(lipp_vec);
    
    Vector3 paike_vec = vec_scale(vec_add(antenn_unitvec, lipp_unitvec), 0.5);
    Vector3 paike_unitvec = vec_normalize(paike_vec);
    
    Vector3 kesk_vec = vec_scale(vec_add(vec_add(maa_unitvec, veen_unitvec), paike_unitvec), 1.0/3.0);
    Vector3 kesk_unitvec = vec_normalize(kesk_vec);
    
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
    
    double alpha = -atan(0.63620006179630684 / 0.11096966105001496);
    double beta = asin(0.76350194216964518);
    double delta = -asin(0.77331312210251590);
    double iota = atan(-0.63355444896004787 / 0.024404);
    
    double start_angle = -11.034;
    double end_angle = -11.03;

    srand(time(NULL));
    
    clock_t start_time = clock();
    
    double best_angle = simulated_annealing(vectors, targets, alpha, beta, delta, iota, start_angle, end_angle);
    
    clock_t end_time = clock();
    double elapsed_time = (double)(end_time - start_time) / CLOCKS_PER_SEC;
    
    double min_value = compute_total_error(best_angle, vectors, targets, alpha, beta, delta, iota);
    
    printf("Angle Range: [%f, %f] degrees\n", start_angle, end_angle);
    printf("Minimum Value: %f\n", min_value);
    printf("Best Angle: %f\n", best_angle);
    printf("\nElapsed time: %f seconds aka %f minutes.\n", elapsed_time, elapsed_time / 60.0);
    
    return 0;
}
