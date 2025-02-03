// Copilot was used to help better comment the code
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>

typedef struct {
    double x;
    double y;
} Vector2D;

// Function prototypes
Vector2D get_force_on_body(const int nstars, const int G, const float e0, int i, double* pos_x, double* pos_y, double* mass);
double get_rij(int i, int j, double* pos_x, double* pos_y);
Vector2D get_position_vector(int i, int j, double* pos_x, double* pos_y);


int main(int argc, char* argv[]) {
    if(argc != 6) {
        fprintf(stderr, "Usage: %s <Number of stars> <input file> <Number of timesteps> <size of timesteps> <graphics>\n", argv[0]);  return 1;
    }
    // Read in command line arguments
    // Make everything const for optimization
    const int nstars = atoi(argv[1]);
    const char* input_file = argv[2];
    const int nsteps = atoi(argv[3]);
    const int stepsize = atoi(argv[4]);
    const bool graphics = atoi(argv[5]); // 0 or 1 for false or true

    const int G = 100/nstars;
    const double e0 = 0.001; // Softening factor 10^-3

    // Read in the input file
    // Input file has structure:
    /*particle 0 position x
    particle 0 position y
    particle 0 mass
    particle 0 velocity x
    particle 0 velocity y
    particle 0 brightness
    particle 1 position x*/

    // Create arrays to store the data
    double* pos_x = (double*) malloc(nstars * sizeof(double));
    double* pos_y = (double*) malloc(nstars * sizeof(double));
    double* mass = (double*) malloc(nstars * sizeof(double));
    double* vel_x = (double*) malloc(nstars * sizeof(double));
    double* vel_y = (double*) malloc(nstars * sizeof(double));
    bool* brightness = (double*) malloc(nstars * sizeof(double));

    // Read in the data
    FILE* file = fopen(input_file, "r");
    if(file == NULL) {
        fprintf(stderr, "Error opening file\n");
        return 1;
    }
    for (int i = 0; i < nstars; i++) {
        fscanf(file, "%lf", &pos_x[i]);
        fscanf(file, "%lf", &pos_y[i]);
        fscanf(file, "%lf", &mass[i]);
        fscanf(file, "%lf", &vel_x[i]);
        fscanf(file, "%lf", &vel_y[i]);
        fscanf(file, "%d", &brightness[i]);
    }
}
// Function to calculate the force on a body
Vector2D get_force_on_body(const int nstars, const int G, const float e0, int i, double* pos_x, double* pos_y, double* mass) {
    Vector2D F;
    F.x = -G * mass[i];
    F.y = -G * mass[i]; 
    for (int j = 0; j < nstars; j++) {
        if (i != j) {
            double rij = get_rij(i, j, pos_x, pos_y);
            Vector2D vector = get_position_vector(i, j, pos_x, pos_y);
            double temp = mass[j] / pow(rij + e0, 3);
            F.x *= temp * vector.x;
            F.y *= temp * vector.y;
        }
    }
    return F;
}

// Function to calculate the distance between two bodies
double get_rij(int i, int j, double* pos_x, double* pos_y) {
    return sqrt(pow(pos_x[i] - pos_x[j], 2) + pow(pos_y[i] - pos_y[j], 2));
}

// Function to calculate the position vector between two bodies
Vector2D get_position_vector(int i, int j, double* pos_x, double* pos_y) {
    Vector2D rij;
    rij.x = pos_x[i] - pos_x[j];
    rij.y = pos_y[i] - pos_y[j];
    return rij;
}
