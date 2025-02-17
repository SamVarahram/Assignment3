// Copilot was used to help better comment the code and add standard error prints
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include <sys/time.h>

// Function prototypes updated to use double arrays for positions and velocities
void get_force_on_body(const int nstars, const double G, const double e0, double* pos_x, double* pos_y, double* mass, double* Fx, double* Fy);
void update_velocity_and_position(const double stepsize, const int nstars, double* vel_x, double* vel_y, double* pos_x, double* pos_y, double* mass, double* Fx, double* Fy);
static double get_wall_seconds();

int main(int argc, char* argv[]) {
    // Start the timer
    double start_time = get_wall_seconds();
    // Check if the correct number of arguments are given
    if(argc != 6) {
        fprintf(stderr, "Usage: %s <Number of stars> <input file> <Number of timesteps> <size of timesteps> <graphics>\n", argv[0]);  return 1;
    }

    // Read in command line arguments
    // Make everything const for optimization
    const int nstars = atoi(argv[1]);
    const char* input_file = argv[2];
    const int nsteps = atoi(argv[3]);
    const double stepsize = atof(argv[4]);
    const int graphics = atoi(argv[5]); // 0 or 1 for false or true
    if (graphics) {
        fprintf(stderr, "Graphics not supported\n");
        return 1;
    }

    const double G = 100./nstars;
    const double e0 = 1e-3; // Softening factor 10^-3

    // Allocate separate arrays for positions and velocities
    double* pos_x      = (double*) malloc(nstars * sizeof(double));
    double* pos_y      = (double*) malloc(nstars * sizeof(double));
    double* mass       = (double*) malloc(nstars * sizeof(double));
    double* vel_x      = (double*) malloc(nstars * sizeof(double));
    double* vel_y      = (double*) malloc(nstars * sizeof(double));
    double* brightness = (double*) malloc(nstars * sizeof(double));
    if (!pos_x || !pos_y || !mass || !vel_x || !vel_y || !brightness) {
        fprintf(stderr, "Error allocating memory\n");
        return 1;
    }

    // Read binary input file. File structure:
    // position x, position y, mass, velocity x, velocity y, brightness (all doubles)
    FILE* file = fopen(input_file, "rb");
    if(file == NULL) {
        fprintf(stderr, "Error opening file\n");
        return 1;
    }
    for (int i = 0; i < nstars; i++) {
        if (fread(&pos_x[i], sizeof(double), 1, file) != 1 ||
            fread(&pos_y[i], sizeof(double), 1, file) != 1 ||
            fread(&mass[i], sizeof(double), 1, file) != 1 ||
            fread(&vel_x[i], sizeof(double), 1, file) != 1 ||
            fread(&vel_y[i], sizeof(double), 1, file) != 1 ||
            fread(&brightness[i], sizeof(double), 1, file) != 1) {
            fprintf(stderr, "Error reading file\n");
            return 1;
        }
    }
    fclose(file);

    double* Fx = (double*) malloc(nstars * sizeof(double));
    double* Fy = (double*) malloc(nstars * sizeof(double));
    for(int i = 0; i < nstars; i++){
        Fx[i] = 0; Fy[i] = 0;
    }

    // Main simulation loop
    for (int time = 0; time < nsteps; time++) {
        get_force_on_body(nstars, G, e0, pos_x, pos_y, mass, Fx, Fy);
        update_velocity_and_position(stepsize, nstars, vel_x, vel_y, pos_x, pos_y, mass, Fx, Fy);
    }

    // Write output binary file with same structure as input
    FILE* output = fopen("result.gal", "wb");
    if(output == NULL) {
        fprintf(stderr, "Error opening file\n");
        return 1;
    }
    for (int i = 0; i < nstars; i++) {
        fwrite(&pos_x[i], sizeof(double), 1, output);
        fwrite(&pos_y[i], sizeof(double), 1, output);
        fwrite(&mass[i], sizeof(double), 1, output);
        fwrite(&vel_x[i], sizeof(double), 1, output);
        fwrite(&vel_y[i], sizeof(double), 1, output);
        fwrite(&brightness[i], sizeof(double), 1, output);
    }
    fclose(output);

    // Free the memory
    free(Fx);
    free(Fy);
    free(pos_x);
    free(pos_y);
    free(mass);
    free(vel_x);
    free(vel_y);
    free(brightness);
    
    // Print the time taken
    double end_time = get_wall_seconds();
    printf("Time taken: %f\n", end_time - start_time);
    return 0;

}

// Function to calculate the force on a body
inline void get_force_on_body(const int nstars, const double G, const double e0, double* pos_x, double* pos_y, double* mass, double* Fx, double* Fy) {
    for (int i = 0; i < nstars; i++) {
        double sumx = 0;
        double sumy = 0;
        for (int j = 0; j < nstars; j++) {
            if (i != j) {

                double dx = pos_x[i] - pos_x[j];
                double dy = pos_y[i] - pos_y[j];
                double r = sqrt(dx * dx + dy * dy) + e0;
                double inv_r3 = 1.0 / (r * r * r);
                double temp = mass[j] * inv_r3;
                sumx += temp * dx;
                sumy += temp * dy;
            }
        }
        Fx[i] = -G * mass[i] * sumx;
        Fy[i] = -G * mass[i] * sumy;
    }
}

// Function to calculate the acceleration of a body
inline void update_velocity_and_position(const double stepsize, const int nstars, double* vel_x, double* vel_y, double* pos_x, double* pos_y, double* mass, double* Fx, double* Fy) {
    for (int i = 0; i < nstars; i++) {
        vel_x[i] += stepsize * Fx[i] / mass[i];
        vel_y[i] += stepsize * Fy[i] / mass[i];
        pos_x[i] += stepsize * vel_x[i];
        pos_y[i] += stepsize * vel_y[i];
    }
}

// Taken from Lab6 Task 02!
static double get_wall_seconds() {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    double seconds = tv.tv_sec + (double)tv.tv_usec / 1000000;
    return seconds;
  }