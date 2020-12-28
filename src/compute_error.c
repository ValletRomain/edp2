/*
Compute the error between result of a method (godunov, rusanov or MUSCL) and exact solution
for different N (number of points of spatial discretization) and compute the duration of method
for different N. Plot the results.

Usage of application :
    compute_error {g|r|m} path_input path_output

path_input : string of path of input file. The input file is filled by the parameters of problem. This file is structure :
    - option_error : norm used to calculate the error
    - option_godunov : equation to resolve
    - xmin : inferior born of spatial interval
    - xmax : superior born of spatial interval
    - cfl : cfl of problem
    - len_liste_N : number of compute of error
    - liste_N : liste of N (number of spatial point) used to compute the error
    - tmax : temporal born
    You can see examples in folder input.

path_output : string of path of output directory where the result is put. This directory is composed by :
    - parameters : file that sums up the parameters of problem
    - plot.dat : result of problem
    - plotcom.gnu : serve to plot the data of plot.dat
    - error.png : graphic of error with different N
    - time.png : graphic of duration with different N

option :
    - g calculate the problem with method of Godunov
    - r calculate the problem with method of Rusanov
    - m calculate the problem with method of MUSCL
*/

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
        raler(0, "usage : compute_error {g|r|m} path_input path_output");

    char * path_input = malloc(CHEMIN_MAX);
    char * path_output = malloc(CHEMIN_MAX);
    strcpy(path_input, argv[optind]);
    optind++;
    strcpy(path_output, argv[optind]);

    parameters_error parerr = {0};

    parameters_error_init_file(&parerr, path_input, path_output, gflag, rflag, mflag);
    parameters_error_compute(&parerr);
    parameters_error_plot(&parerr);
    parameters_error_free(&parerr);

    exit(0);
}
