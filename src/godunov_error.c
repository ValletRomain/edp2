#include "godunov.h"
#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>

#include "raler.c"
#include "function.c"

int main(int argc, char * argv[]){

    if (argc != 3){
        raler(0, "usage : godunov_error init_file output_path");
    }

    godunov_error gderr = {0};

    godunov_error_init_file(&gderr, argv[1]);
    godunov_error_compute(&gderr);
    godunov_error_plot(&gderr, argv[2]);
    godunov_error_free(&gderr);

    exit(0);
}
