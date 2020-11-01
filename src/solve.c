/* Application godunov : resolve the schema of godunov
usage : godunov -v path_input path_output
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


    parameters par = {0};
    
    parameters_init_file(&par, path_input, path_output, aflag);
    parameters_solve(&par, 1);
    parameters_plot(&par);
    parameters_free(&par);

    exit(0);
}
