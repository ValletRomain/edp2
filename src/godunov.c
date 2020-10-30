/* Application godunov : resolve the schema of godunov
usage : godunov init_file output_path
init_file : file of parameters of problem
output_path : path where the output folder is create and fill the result of problem

The new output folder contains :
- the path is output_path/name_file with name_file is the name of file of init_file
- plot.dat contains the absissa, the numeric solution and if keept_solution=1 the exact solution
- plotcom.gnu a gnuplot script to create graphe.png by plot.dat
- graphe.png the graphe of numeric solution and if keept_solution=1 the graphe of exate solution
*/

#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>

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
