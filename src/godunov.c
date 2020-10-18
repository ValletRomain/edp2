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
        raler(0, "usage : godunov init_file output_path");
    }
    
    godunov gd = {0};
    
    godunov_init_file(&gd, argv[1]);
    godunov_solve(&gd, 1);
    godunov_plot(&gd, argv[2]);
    godunov_free(&gd);

    exit(0);
}
