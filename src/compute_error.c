#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>

#include "function.c"

int main(int argc, char * argv[]){

    if (argc != 3){
        raler(0, "usage : godunov_error init_file output_path");
    }

    parameters_error parerr = {0};
    
    parameters_error_init_file(&parerr, argv[1]);
    parameters_error_compute(&parerr);
    parameters_error_plot(&parerr, argv[2]);
    parameters_error_free(&parerr);

    exit(0);
}
