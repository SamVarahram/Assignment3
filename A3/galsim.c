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
    const int graphics = atoi(argv[5]); // 0 or 1 for false or true
    if (graphics) {
        fprintf(stderr, "Graphics not supported\n");
        return 1;
    }

    const int G = 100/nstars;
    const double e0 = 1e-3; // Softening factor 10^-3

    // Create arrays to store the data and Vector2D to store the forces
    Vector2D* position = (Vector2D*) malloc(nstars * sizeof(Vector2D));
    double* mass = (double*) malloc(nstars * sizeof(double));
    Vector2D* velocity = (Vector2D*) malloc(nstars * sizeof(Vector2D));
    double* brightness = (double*) malloc(nstars * sizeof(double));

    if (position == NULL || mass == NULL || velocity == NULL || brightness == NULL) {
        fprintf(stderr, "Error allocating memory\n");
        return 1;
    }
    // Initialize arrays to zero
    memset(position, 0, nstars * sizeof(Vector2D));
    memset(mass, 0, nstars * sizeof(double));
    memset(velocity, 0, nstars * sizeof(Vector2D));
    memset(brightness, 0, nstars * sizeof(double));

    // Read in the data
    FILE* file = fopen(input_file, "rb");
    if(file == NULL) {
        fprintf(stderr, "Error opening file\n");
        free(position);
        free(mass);
        free(velocity);
        free(brightness);
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
        if (fread(&position[i], sizeof(Vector2D), 1, file) != 1 ||
            fread(&mass[i], sizeof(double), 1, file) != 1 ||
            fread(&velocity[i], sizeof(Vector2D), 1, file) != 1 ||
            fread(&brightness[i], sizeof(double), 1, file) != 1) {
            fprintf(stderr, "Error reading file\n");
            return 1;
        }
        // Print values for debugging
        // printf("Star %d: Position (%f, %f), Mass %f, Velocity (%f, %f), Brightness %f\n",
        //    i, position[i].x, position[i].y, mass[i], velocity[i].x, velocity[i].y, brightness[i]);
    }
    fclose(file);




    Vector2D F = {0.0, 0.0};
    // Loop over the timesteps
    for (int time = 0; time < nsteps; time++) {
        for (int i = 0; i < nstars; i++) {
            F = get_force_on_body(nstars, G, e0, i, position, mass);
            update_velocity_and_position(i, stepsize, velocity, position, &F, mass);
        }
    }
    


    // Output the data in a binary file
    FILE* output = fopen("result.gal", "wb");
    if(output == NULL) {
        fprintf(stderr, "Error opening file\n");
        return 1;
    }
    for (int i = 0; i < nstars; i++) {
        fwrite(&position[i], sizeof(Vector2D), 1, output);
        fwrite(&mass[i], sizeof(double), 1, output);
        fwrite(&velocity[i], sizeof(Vector2D), 1, output);
        fwrite(&brightness[i], sizeof(double), 1, output);
    }
    fclose(output);

    // Free the memory
    free(position);
    free(mass);
    free(velocity);
    free(brightness);
    

}
// Function to calculate the force on a body
Vector2D get_force_on_body(const int nstars, const int G, const float e0, int i, Vector2D* position, double* mass) {
    Vector2D F = {0.0, 0.0}; 
    for (int j = 0; j < nstars; j++) {
        if (i != j) {

            double dx = position[i].x - position[j].x;
            double dy = position[i].y - position[j].y;
            double rij = sqrt(dx * dx + dy * dy);
            
            
            double temp = mass[j] / ((rij + e0) * (rij + e0) * (rij + e0)); // pow(rij + e0, 3)
            F.x += temp * dx;
            F.y += temp * dy;
        }
    }
    F.x *= -G * mass[i];
    F.y *= -G * mass[i];
    return F;
}

// Function to calculate the acceleration of a body
void update_velocity_and_position(int i, const int stepsize, Vector2D* velocity, Vector2D* position, Vector2D* F, double* mass) {
    velocity[i].x += stepsize * F->x / mass[i];
    velocity[i].y += stepsize * F->y / mass[i];
    position[i].x += stepsize * velocity[i].x;
    position[i].y += stepsize * velocity[i].y;
}