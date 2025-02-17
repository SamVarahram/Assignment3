// Copilot was used to help better comment the code and add standard error prints
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include <sys/time.h>

// Struct to store the position of a body
typedef struct {
    double x;
    double y;
} Vector2D;

// Function prototypes
void get_force_on_body(const int nstars, const double G, const double e0, Vector2D* position, double* mass, double* Fx, double* Fy);
void update_velocity_and_position(const double stepsize, const int nstars, Vector2D* velocity, Vector2D* position, double* mass, double* Fx, double* Fy);
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
    }
    fclose(file);




    double* Fx = (double*) malloc(nstars * sizeof(double));
    double* Fy = (double*) malloc(nstars * sizeof(double));
    memset(Fx, 0, nstars * sizeof(double));
    memset(Fy, 0, nstars * sizeof(double));

    // Loop over the timesteps
    for (int time = 0; time < nsteps; time++) {
        get_force_on_body(nstars, G, e0, position, mass, Fx, Fy);
        update_velocity_and_position(stepsize, nstars, velocity, position, mass, Fx, Fy);
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
    free(Fx);
    free(Fy);
    free(position);
    free(mass);
    free(velocity);
    free(brightness);
    
    // Print the time taken
    double end_time = get_wall_seconds();
    printf("Time taken: %f\n", end_time - start_time);
    return 0;

}

// Function to calculate the force on a body
void get_force_on_body(const int nstars, const double G, const double e0, Vector2D* position, double* mass, double* Fx, double* Fy) {
    for (int i = 0; i < nstars; i++) {
        double sumx = 0;
        double sumy = 0;
        for (int j = 0; j < nstars; j++) {
            if (i != j) {

                double dx = position[i].x - position[j].x;
                double dy = position[i].y - position[j].y;
                double rije0 = sqrt(dx * dx + dy * dy) + e0;
                double pow_rije0 = 1.0 / (rije0 * rije0 * rije0);
                
                double temp = mass[j] * pow_rije0; // pow(rij + e0, 3)
                sumx += temp * dx;
                sumy += temp * dy;
            }
        }
        Fx[i] = sumx * -G * mass[i];
        Fy[i] = sumy * -G * mass[i];
    }
}

// Function to calculate the acceleration of a body
void update_velocity_and_position(const double stepsize, const int nstars, Vector2D* velocity, Vector2D* position, double* mass, double* Fx, double* Fy) {
    for (int i = 0; i < nstars; i++) {
        velocity[i].x += stepsize * Fx[i] / mass[i];
        velocity[i].y += stepsize * Fy[i] / mass[i];
        position[i].x += stepsize * velocity[i].x;
        position[i].y += stepsize * velocity[i].y;
    }
}

// Taken from Lab6 Task 02!
static double get_wall_seconds() {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    double seconds = tv.tv_sec + (double)tv.tv_usec / 1000000;
    return seconds;
  }