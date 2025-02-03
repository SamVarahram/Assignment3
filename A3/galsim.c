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
Vector2D get_force_on_body(const int nstars, const int G, const float e0, int i, Vector2D* position, double* mass);
double get_rij(int i, int j, Vector2D* position);
Vector2D get_position_vector(int i, int j, Vector2D* position);
void update_velocity_and_position(int i, const int stepsize, Vector2D* velocity, Vector2D* position, Vector2D* F, double* mass);



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

    // Create arrays to store the data and Vector2D to store the forces
    Vector2D* position = (Vector2D*) malloc(nstars * sizeof(Vector2D));
    double* mass = (double*) malloc(nstars * sizeof(double));
    Vector2D* velocity = (Vector2D*) malloc(nstars * sizeof(Vector2D));
    bool* brightness = (double*) malloc(nstars * sizeof(double));

    if (position == NULL || mass == NULL || velocity == NULL || brightness == NULL) {
        fprintf(stderr, "Error allocating memory\n");
        return 1;
    }

    // Read in the data
    FILE* file = fopen(input_file, "r");
    if(file == NULL) {
        fprintf(stderr, "Error opening file\n");
        return 1;
    }

    // Read in the input file
    // Input file has structure:
    /*particle 0 position x
    particle 0 position y
    particle 0 mass
    particle 0 velocity x
    particle 0 velocity y
    particle 0 brightness
    particle 1 position x*/
    // The input file is binary
    for (int i = 0; i < nstars; i++) {
        fread(&position[i], sizeof(Vector2D), 1, file);
        fread(&mass[i], sizeof(double), 1, file);
        fread(&velocity[i], sizeof(Vector2D), 1, file);
        fread(&brightness[i], sizeof(bool), 1, file);
    }
    fclose(file);

    // Output the data in a binary file
    FILE* output = fopen("output.bin", "wb");
    if(output == NULL) {
        fprintf(stderr, "Error opening file\n");
        return 1;
    }
    for (int i = 0; i < nstars; i++) {
        fwrite(&position[i], sizeof(Vector2D), 1, output);
        fwrite(&mass[i], sizeof(double), 1, output);
        fwrite(&velocity[i], sizeof(Vector2D), 1, output);
        fwrite(&brightness[i], sizeof(bool), 1, output);
    }

    // Free the memory
    free(position);
    free(mass);
    free(velocity);
    free(brightness);
    // Jag bytte till test

}
// Function to calculate the force on a body
Vector2D get_force_on_body(const int nstars, const int G, const float e0, int i, Vector2D* position, double* mass) {
    Vector2D F;
    F.x = -G * mass[i];
    F.y = -G * mass[i]; 
    for (int j = 0; j < nstars; j++) {
        if (i != j) {
            double rij = get_rij(i, j, position);
            Vector2D vector = get_position_vector(i, j, position);
            double temp = mass[j] / pow(rij + e0, 3);
            F.x *= temp * vector.x;
            F.y *= temp * vector.y;
        }
    }
    return F;
}

// Function to calculate the distance between two bodies
double get_rij(int i, int j, Vector2D* position) {
    return sqrt(pow(position[i].x - position[j].x, 2) + pow(position[i].y - position[j].y, 2));
}

// Function to calculate the position vector between two bodies
Vector2D get_position_vector(int i, int j, Vector2D* position) {
    Vector2D rij;
    rij.x = position[i].x - position[j].x;
    rij.y = position[i].y - position[j].y;
    return rij;
}

// Function to calculate the acceleration of a body
void update_velocity_and_position(int i, const int stepsize, Vector2D* velocity, Vector2D* position, Vector2D* F, double* mass) {
    velocity[i+1].x += stepsize * F[i].x / mass[i];
    velocity[i+1].y += stepsize * F[i].y / mass[i];
    position[i+1].x += stepsize * velocity[i+1].x;
    position[i+1].y += stepsize * velocity[i+1].y;
}

int main_update(int nstars, int G, float e0, int i, Vector2D* position, double* mass, Vector2D* velocity, int stepsize) {
    Vector2D F = get_force_on_body(nstars, G, e0, i, position, mass);
    update_velocity_and_position(i, stepsize, velocity, position, &F, mass);
    return 0;
}