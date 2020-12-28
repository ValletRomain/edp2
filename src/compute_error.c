#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>

#include "function.c"

int main(int argc, char * argv[]){

    int c, errflag=0, gflag=0, rflag=0, mflag=0;
    extern char *optarg;
    extern int optind;

    while ( (c=getopt(argc, argv, "grm"))!=-1 ){
        switch (c) {        
        case 'g':
            gflag++;
            break;

        case 'r':
            rflag++;
            break;

        case 'm':
            mflag++;
            break;

        default:
            errflag++;
            break;
        }
    }

    if ((gflag==0) && (rflag==0) && (mflag==0))
        errflag++;

    if ((argc-optind) != 2)
        errflag++;
    
    if (errflag)
        raler(0, "usage : compute_error {g|r|gr} path_input path_output");

    char * path_input = malloc(CHEMIN_MAX);
    char * path_output = malloc(CHEMIN_MAX);
    strcpy(path_input, argv[optind]);
    optind++;
    strcpy(path_output, argv[optind]);

    parameters_error parerr = {0};

    parameters_error_init_file(&parerr, path_input, path_output, gflag, rflag, mflag);

    printf("len_N={%d}, g=%d, r=%d, m=%d\n", parerr.len_liste_N, parerr.option_godunov,
                    parerr.option_rusanov, parerr.option_muscl);

    parameters_error_compute(&parerr);
    parameters_error_plot(&parerr);
    parameters_error_free(&parerr);

    exit(0);
}
