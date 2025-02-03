#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

int main(int argc, char** argv[]) {
    if(argc != 6) {
        fprintf(stderr, "Usage: %s <Number of stars> <input file> <Number of timesteps> <size of timesteps> <graphics>\n", argv[0]);  return 1;
    }
    int nstars = atoi(argv[1]);
    char** input_file = argv[2];
    int nsteps = atoi(argv[3]);
    int stepsize = atoi(argv[4]);
    bool graphics = atoi(argv[5]); // 0 or 1 for false or true
}

