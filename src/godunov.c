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

    int c, errflag=0, aflag=0;
    extern char *optarg;
    extern int optind;

    while ( (c=getopt(argc, argv, "a"))!=-1 ){
        switch (c)
        {
        case 'a':
            aflag++;
            break;
        
        default:
            errflag++;
            break;
        }
    }

    if ((argc-optind) != 2)
        errflag++;
    
    if (errflag)
        raler(0, "usage : godunov [a] path_input path_output");

    char * path_input = malloc(CHEMIN_MAX);
    char * path_output = malloc(CHEMIN_MAX);
    strcpy(path_input, argv[optind]);
    optind++;
    strcpy(path_output, argv[optind]);


    godunov gd = {0};
    
    godunov_init_file(&gd, path_input, path_output, aflag);
    godunov_solve(&gd, 1);
    godunov_plot(&gd);
    godunov_free(&gd);

    exit(0);
}
