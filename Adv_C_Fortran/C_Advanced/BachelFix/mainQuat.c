/*
Compile and execute with:
    $ gcc main.c -o main -lm
    $ ./main
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/stat.h>

#define PI 3.14159265358979323846  // Define constant PI for mathematical calculations
#define CLOCK_MONOTONIC 1

void ensure_directories_exist();
void angle_to_dir(double EL, double AZ, double dir[3]);
void print_vector(const char* label, double vec[3]);
void normalize_vector(double vec[3]);
void extrinsic_to_intrinsic(double x, double y, double z, double result[3]);
void sequential_rotations(double x, double y, double z, double result[3]);
double azimuth(double z);
double elevation(double x);
void vec_to_angles(double x, double y, double z, double angles[2]);
void centralize(double x, double y, double z, double result[3]);
void load_data_from_csv(const char* filename, double*** data, int* rows, int* cols);
void free_data(double** data, int rows);

typedef struct {
    double w, x, y, z;
} Quaternion;

Quaternion quaternion_from_axis_angle(double x, double y, double z, double angle);
Quaternion quaternion_multiply(Quaternion q1, Quaternion q2);
void quaternion_rotate_vector(Quaternion q, double v[3], double result[3]);

int main() {
    // Timing start
    struct timespec start, stop;
    clock_gettime(CLOCK_MONOTONIC, &start); // Get starting time

    ensure_directories_exist(); // Ensure the output directories exist

    // Define horizontal angles for celestial bodies
    double AZ_venus_HOR = 75.806271, EL_venus_HOR = 69.278559;
    double AZ_earth_HOR = 94.647066, EL_earth_HOR = 66.507475;
    double AZ_sun_HOR   = 89.324328, EL_sun_HOR   = 13.697882;

    double HOR_venus[3], HOR_earth[3], HOR_sun[3], HOR_avrg[3], HOR_final_avrg[3];

    // Convert angles to direction vectors
    angle_to_dir(EL_venus_HOR, AZ_venus_HOR, HOR_venus);
    angle_to_dir(EL_earth_HOR, AZ_earth_HOR, HOR_earth);
    angle_to_dir(EL_sun_HOR, AZ_sun_HOR, HOR_sun);

    // Compute the average direction vector
    for (int i = 0; i < 3; ++i) {
        HOR_avrg[i] = (HOR_venus[i] + HOR_earth[i] + HOR_sun[i]) / 3.0;
    }
    normalize_vector(HOR_avrg); // Normalize the average vector
    for (int i = 0; i < 3; ++i) {
        HOR_final_avrg[i] = HOR_avrg[i];
    }
    print_vector("Average Horizons azimuth vector", HOR_final_avrg);

    // Load data from a CSV file
    double** data;
    int rows, cols;
    load_data_from_csv("InputFiles/coord_orien_all.csv", &data, &rows, &cols);

    // Print initial data
    printf("\nInitial data:\nID           X            Y            Z            RX           RY           RZ\n");
    for (int i = 0; i < 5; ++i) {
        printf("%d", (int)data[i][0]);
        for (int j = 1; j < cols; ++j) {
            printf("%12.5f ", data[i][j]);
        }
        printf("\n");
    }

    // Compute elevation and azimuth for rotation
    double modRX[rows], modRZ[rows];
    for (int i = 0; i < rows; ++i) {
        modRX[i] = elevation(data[i][4]);
        modRZ[i] = azimuth(data[i][6]);
    }

    double nurgad[rows][2];
    for (int i = 0; i < rows; ++i) {
        nurgad[i][0] = modRX[i];
        nurgad[i][1] = modRZ[i];
    }

    // Convert rotation angles to direction vectors
    double dirANG[rows][3];
    for (int i = 0; i < rows; ++i) {
        angle_to_dir(nurgad[i][0], nurgad[i][1], dirANG[i]);
    }

    // Apply sequential rotations to transform position
    double trans_pos_old[rows][3];
    for (int i = 0; i < rows; ++i) {
        sequential_rotations(dirANG[i][0], dirANG[i][1], dirANG[i][2], trans_pos_old[i]);
    }

    // Convert extrinsic angles to intrinsic system
    double init_dataANG[rows][3];
    for (int i = 0; i < rows; ++i) {
        extrinsic_to_intrinsic(data[i][4], data[i][5], data[i][6], init_dataANG[i]);
    }

    // Apply sequential rotations to intrinsic angles
    double transfANG[rows][3];
    for (int i = 0; i < rows; ++i) {
        sequential_rotations(init_dataANG[i][0], init_dataANG[i][1], init_dataANG[i][2], transfANG[i]);
    }

    // Convert transformed vectors to angles
    double newANG[rows][2];
    for (int i = 0; i < rows; ++i) {
        vec_to_angles(transfANG[i][0], transfANG[i][1], transfANG[i][2], newANG[i]);
    }

    // Apply sequential rotations to original coordinates
    double coord_temp[rows][3];
    for (int i = 0; i < rows; ++i) {
        sequential_rotations(data[i][1], data[i][2], data[i][3], coord_temp[i]);
    }

    // Centralize coordinates
    double newCOORD[rows][3];
    for (int i = 0; i < rows; ++i) {
        centralize(coord_temp[i][0], coord_temp[i][1], coord_temp[i][2], newCOORD[i]);
    }

    // Print transformed values
    printf("\nTransformed values:\nID           X            Y            Z            EL           AZ\n");
    for (int i = 0; i < 5; ++i) {
        printf("%d", (int)data[i][0]);
        for (int j = 0; j < 3; ++j) {
            printf("%12.5f ", newCOORD[i][j]);
        }
        for (int j = 0; j < 2; ++j) {
            printf("%12.5f ", newANG[i][j]);
        }
        printf("\n");
    }

    // Save transformed data to CSV files
    FILE* output_file1 = fopen("OutputFiles/camera_posANDdir.csv", "w");
    fprintf(output_file1, "ID, X_(m), Y_(m), Z_(m), EL_new_(DEG), AZ_new_(DEG)\n");
    for (int i = 0; i < rows; ++i) {
        fprintf(output_file1, "%d", (int)data[i][0]);
        for (int j = 0; j < 3; ++j) {
            fprintf(output_file1, ",%22.16f", newCOORD[i][j]);
        }
        for (int j = 0; j < 2; ++j) {
            fprintf(output_file1, ",%22.16f", newANG[i][j]);
        }
        fprintf(output_file1, "\n");
    }
    fclose(output_file1);

    FILE* output_file2 = fopen("OutputFiles/camera_dirVECS.csv", "w");
    fprintf(output_file2, "Vx, Vy, Vz\n");
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < 3; ++j) {
            fprintf(output_file2, "%22.16f", init_dataANG[i][j]);
            if (j < 2) fprintf(output_file2, ",");
        }
        fprintf(output_file2, "\n");
    }
    fclose(output_file2);

    // Timing end
    clock_gettime(CLOCK_MONOTONIC, &stop); // Get end time

    // Calculate the elapsed time in seconds
    double time_taken = (stop.tv_sec - start.tv_sec) * 1e9;
    time_taken = (time_taken + (stop.tv_nsec - start.tv_nsec)) * 1e-9;
    printf("\nExecution time: %.3e seconds.\n", time_taken);

    free_data(data, rows); // Free allocated memory

    return 0;
}

//=================================================================================
//=============================End of main code====================================
//=================================================================================


/**
 * -------------------------------------------------------
 * Function Description:
 * -------------------------------------------------------
 * This function checks whether a directory named "OutputFiles" exists. If not, it creates it.
 */
void ensure_directories_exist() {
    struct stat st = {0};

    if (stat("OutputFiles", &st) == -1) {
        mkdir("OutputFiles", 0700); // Create directory with full read/write permissions
    }
}

/**
 * -------------------------------------------------------
 * Function Description:
 * -------------------------------------------------------
 * Converts azimuth and elevation angles to a direction vector.
 * 
 * Parameters:
 * - EL: Elevation angle in degrees.
 * - AZ: Azimuth angle in degrees.
 * - dir: Array to store the resulting direction vector.
 */
void angle_to_dir(double EL, double AZ, double dir[3]) {
    double el_rad = EL * PI / 180.0; // Convert elevation to radians
    double az_rad = AZ * PI / 180.0; // Convert azimuth to radians

    dir[0] = cos(el_rad) * cos(az_rad);
    dir[1] = cos(el_rad) * sin(az_rad);
    dir[2] = sin(el_rad);
}

/**
 * -------------------------------------------------------
 * Function Description:
 * -------------------------------------------------------
 * Prints a vector with a label.
 * 
 * Parameters:
 * - label: The label for the vector.
 * - vec: The vector to be printed.
 */
void print_vector(const char* label, double vec[3]) {
    printf("%s: [%.6f, %.6f, %.6f]\n", label, vec[0], vec[1], vec[2]);
}

/**
 * -------------------------------------------------------
 * Function Description:
 * -------------------------------------------------------
 * Normalizes a 3D vector.
 * 
 * Parameters:
 * - vec: The vector to be normalized.
 */
void normalize_vector(double vec[3]) {
    double norm = sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);
    vec[0] /= norm;
    vec[1] /= norm;
    vec[2] /= norm;
}

/**
 * -------------------------------------------------------
 * Function Description:
 * -------------------------------------------------------
 * Converts extrinsic to intrinsic rotations.
 * 
 * Parameters:
 * - x, y, z: Extrinsic rotation angles.
 * - result: Array to store the intrinsic rotation results.
 */
void extrinsic_to_intrinsic(double x, double y, double z, double result[3]) {
    // Using quaternions for more robust rotations
    Quaternion qx = quaternion_from_axis_angle(1, 0, 0, x);
    Quaternion qy = quaternion_from_axis_angle(0, 1, 0, y);
    Quaternion qz = quaternion_from_axis_angle(0, 0, 1, z);

    Quaternion q = quaternion_multiply(qz, quaternion_multiply(qy, qx));

    result[0] = q.x;
    result[1] = q.y;
    result[2] = q.z;
}

/**
 * -------------------------------------------------------
 * Function Description:
 * -------------------------------------------------------
 * Applies sequential rotations to a vector using quaternions.
 * 
 * Parameters:
 * - x, y, z: Rotation angles.
 * - result: Array to store the rotated vector.
 */
void sequential_rotations(double x, double y, double z, double result[3]) {
    double vec[3] = {x, y, z};

    // Using quaternions for more robust rotations
    Quaternion qx = quaternion_from_axis_angle(1, 0, 0, x);
    Quaternion qy = quaternion_from_axis_angle(0, 1, 0, y);
    Quaternion qz = quaternion_from_axis_angle(0, 0, 1, z);

    Quaternion q = quaternion_multiply(qz, quaternion_multiply(qy, qx));

    quaternion_rotate_vector(q, vec, result);
}

/**
 * -------------------------------------------------------
 * Function Description:
 * -------------------------------------------------------
 * Computes the azimuth angle from a given vector.
 * 
 * Parameters:
 * - z: The z-component of the vector.
 * 
 * Returns:
 * - The azimuth angle in degrees.
 */
double azimuth(double z) {
    return atan2(z, sqrt(1 - z * z)) * 180.0 / PI;
}

/**
 * -------------------------------------------------------
 * Function Description:
 * -------------------------------------------------------
 * Computes the elevation angle from a given vector.
 * 
 * Parameters:
 * - x: The x-component of the vector.
 * 
 * Returns:
 * - The elevation angle in degrees.
 */
double elevation(double x) {
    return asin(x) * 180.0 / PI;
}

/**
 * -------------------------------------------------------
 * Function Description:
 * -------------------------------------------------------
 * Converts a direction vector to elevation and azimuth angles.
 * 
 * Parameters:
 * - x, y, z: The components of the direction vector.
 * - angles: Array to store the resulting elevation and azimuth angles.
 */
void vec_to_angles(double x, double y, double z, double angles[2]) {
    angles[0] = atan2(y, x) * 180.0 / PI;  // Azimuth
    angles[1] = atan2(z, sqrt(x * x + y * y)) * 180.0 / PI;  // Elevation
}

/**
 * -------------------------------------------------------
 * Function Description:
 * -------------------------------------------------------
 * Centralizes coordinates by subtracting the mean.
 * 
 * Parameters:
 * - x, y, z: The original coordinates.
 * - result: Array to store the centralized coordinates.
 */
void centralize(double x, double y, double z, double result[3]) {
    result[0] = x - 0.0;  // Replace 0.0 with the mean if needed
    result[1] = y - 0.0;  // Replace 0.0 with the mean if needed
    result[2] = z - 0.0;  // Replace 0.0 with the mean if needed
}

/**
 * -------------------------------------------------------
 * Function Description:
 * -------------------------------------------------------
 * Loads data from a CSV file.
 * 
 * Parameters:
 * - filename: The name of the CSV file.
 * - data: Pointer to store the loaded data.
 * - rows: Pointer to store the number of rows.
 * - cols: Pointer to store the number of columns.
 */
void load_data_from_csv(const char* filename, double*** data, int* rows, int* cols) {
    FILE* file = fopen(filename, "r");
    if (!file) {
        perror("Could not open file");
        exit(EXIT_FAILURE);
    }

    char line[1024];
    int row = 0;
    int col = 0;
    int header_processed = 0;

    // First pass: determine the number of rows and columns
    while (fgets(line, sizeof(line), file)) {
        if (!header_processed) {
            // Skip the header row
            header_processed = 1;
            // Count the number of columns based on the first row
            char* token = strtok(line, ",");
            while (token) {
                col++;
                token = strtok(NULL, ",");
            }
            continue;
        }
        row++;
    }

    *rows = row;
    *cols = col;

    // Allocate memory for the data array
    *data = (double**)malloc(*rows * sizeof(double*));
    for (int i = 0; i < *rows; i++) {
        (*data)[i] = (double*)malloc(*cols * sizeof(double));
    }

    rewind(file);

    // Skip the header row
    fgets(line, sizeof(line), file);

    // Second pass: populate the data array
    row = 0;
    while (fgets(line, sizeof(line), file)) {
        col = 0;
        char* token = strtok(line, ",");
        while (token) {
            (*data)[row][col] = atof(token);
            token = strtok(NULL, ",");
            col++;
        }
        row++;
    }

    fclose(file);
}

/**
 * -------------------------------------------------------
 * Function Description:
 * -------------------------------------------------------
 * Frees the allocated memory for data.
 * 
 * Parameters:
 * - data: The data to be freed.
 * - rows: The number of rows in the data.
 */
void free_data(double** data, int rows) {
    for (int i = 0; i < rows; ++i) {
        free(data[i]);
    }
    free(data);
}

/**
 * -------------------------------------------------------
 * Function Description:
 * -------------------------------------------------------
 * Creates a quaternion from an axis-angle representation.
 * 
 * Parameters:
 * - x, y, z: The axis components.
 * - angle: The rotation angle in degrees.
 * 
 * Returns:
 * - The resulting quaternion.
 */
Quaternion quaternion_from_axis_angle(double x, double y, double z, double angle) {
    Quaternion q;
    double half_angle = angle * PI / 360.0;
    double s = sin(half_angle);
    q.w = cos(half_angle);
    q.x = x * s;
    q.y = y * s;
    q.z = z * s;
    return q;
}

/**
 * -------------------------------------------------------
 * Function Description:
 * -------------------------------------------------------
 * Multiplies two quaternions.
 * 
 * Parameters:
 * - q1, q2: The quaternions to be multiplied.
 * 
 * Returns:
 * - The resulting quaternion.
 */
Quaternion quaternion_multiply(Quaternion q1, Quaternion q2) {
    Quaternion q;
    q.w = q1.w * q2.w - q1.x * q2.x - q1.y * q2.y - q1.z * q2.z;
    q.x = q1.w * q2.x + q1.x * q2.w + q1.y * q2.z - q1.z * q2.y;
    q.y = q1.w * q2.y - q1.x * q2.z + q1.y * q2.w + q1.z * q2.x;
    q.z = q1.w * q2.z + q1.x * q2.y - q1.y * q2.x + q1.z * q2.w;
    return q;
}

/**
 * -------------------------------------------------------
 * Function Description:
 * -------------------------------------------------------
 * Rotates a vector by a quaternion.
 * 
 * Parameters:
 * - q: The quaternion.
 * - v: The vector to be rotated.
 * - result: Array to store the rotated vector.
 */
void quaternion_rotate_vector(Quaternion q, double v[3], double result[3]) {
    Quaternion qv = {0, v[0], v[1], v[2]};
    Quaternion q_conjugate = {q.w, -q.x, -q.y, -q.z};
    Quaternion qv_rotated = quaternion_multiply(quaternion_multiply(q, qv), q_conjugate);

    result[0] = qv_rotated.x;
    result[1] = qv_rotated.y;
    result[2] = qv_rotated.z;
}
