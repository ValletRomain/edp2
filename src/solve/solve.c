/*
Solve the problem with a method (godunov, rusanov or MUSCL). Plot the result.

Usage of application :
    solve {g|r|m|e} path_input path_output

path_input : string of path of input file. The input file is filled by the parameters of problem. This file is structure :
    - option_godunov : equation to resolve
    - xmin : inferior born of spatial interval
    - xmax : superior born of spatial interval
    - cfl : cfl of problem
    - N : number of spatial point
    - tmax : temporal born
    You can see examples in folder input.

path_output : string of path of output directory where the result is put. This directory is composed by :
    - parameters : file that sums up the parameters of problem
    - plot.dat : result of problem
    - plotcom.gnu : serve to plot the data of plot.dat
    - plot.png : graphic of results in relation to spatial at tmax.

option :
    - g solve with method of Godunov
    - r solve with method of Rusanov
    - m solve with method of MUSCL
    - e compute the exact solution
*/

#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>

#include "function.c"

int main(int argc, char * argv[]){

    int c, errflag=0, aflag=0, gflag=0, rflag=0, mflag=0, eflag=0;
    extern char *optarg;
    extern int optind;

    aflag = 0;

    while ( (c=getopt(argc, argv, "grme"))!=-1 ){
        switch (c)
        {
/*      case 'a':
            aflag++;
            break;*/
        
        case 'g':
            gflag++;
            break;

        case 'r':
            rflag++;
            break;

        case 'm':
            mflag++;
            break;
        
        case 'e':
            eflag++;
            break;

        default:
            errflag++;
            break;
        }
    }

    if ((gflag==0) && (rflag==0) && (mflag==0) && (eflag==0))
        errflag++;

    if ((argc-optind) != 2)
        errflag++;
    
    if (errflag)
        raler(0, "usage : solve {g|r|gr} [a] path_input path_output");

    char * path_input = malloc(CHEMIN_MAX);
    char * path_output = malloc(CHEMIN_MAX);
    strcpy(path_input, argv[optind]);
    optind++;
    strcpy(path_output, argv[optind]);


    parameters par = {0};
    
    parameters_init_file(&par,
                        path_input, path_output,
                        gflag, rflag, mflag, eflag);
    
    if (gflag)
        godunov_solve(&par, 1);
    if (rflag)
        rusanov_solve(&par, 1);
    //if (mflag)
    //    muscl_solve(&par, 1);

    w_to_hu(&par);

    parameters_plot(&par);
    
    parameters_free(&par);

    exit(0);
}
