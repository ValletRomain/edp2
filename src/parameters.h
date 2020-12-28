#ifndef _PARAMETERS_H
#define _PARAMETERS_H

typedef struct parameters{

    // Name of different file or folder
    char * name_file; //name of file of input

    char * path_input; // path of input file
    char * complete_path_output; // path of output folder (output_path/name_file)

    // Option
    int option_solexacte; //keept the solution (=1) or not (=0)
    int option_animation;
    int option_godunov;
    int option_rusanov;
    int option_muscl;
    char * option_equation; // equation is resolve (burgers, transport...)

    // Parametre du probleme
    int N; // m dimension of problem, N number of space point
    double xmin, xmax;
    double cfl;
    double dt, dx;
    double tmax;

    // Parameters of equation
    double (*plambda_ma)(double);

    // Flux for godunov
    double (*pfluxnum_gd)(double, double);

    // Flux for rusanov
    double (*pfluxnum_ru)(double, double);

    // Border condition
    double (*pboundary_spatial)(double);
    double (*pboundary_temporal_left)(double, double);
    double (*pboundary_temporal_right)(double, double);

    // Solution exacte
    double (*psolexacte)(double, double);


    // Resultats du probleme
    unsigned long time_gd;
    unsigned long time_ru;
    unsigned long time_muscl;
    
    double *xi; // centre des milieux des cellules
    
    // Godunov
    double *un; // solution a l'instant n
    double *unp1; // solution a l'instant n+1

    // Rusanov
    double *vn; // solution a l'instant n
    double *vnp1; // solution a l'instant n+1

    // MUSCL
    double *wn;
    double *wnp1;

    // solution exact
    double *sol;

    int int_tnow_gd;
    int int_tnow_ru;

} parameters;

typedef struct parameters_error{

    // Name of different file or folder
    char * name_file;

    char * path_input; // path of input file
    char * complete_path_output; // path of output folder (output_path/name_file)

    char * option_error;
    char * option_equation;
    int option_godunov;
    int option_rusanov;
    int option_muscl;

    // Parametre du probleme
    //int m; // nombre de variables conservatives, nombre de cellules
    double dt, dx; // pas de temps, pas d'espace
    double xmin, xmax; // bornes de l'intervalles
    double cfl; // vmax : dt/dx
    double tmax;

    int len_liste_N;
    int * liste_N;

    double (*perror)(int, double*, double*); // error (norm_L1, L2, inf...)

    double * liste_error_gd;
    unsigned long * liste_time_gd;

    double * liste_error_ru;
    unsigned long * liste_time_ru;

    double *liste_error_muscl;
    unsigned long *liste_time_muscl;

} parameters_error;

void parameters_init(parameters *par,
                    int option_godunov, int option_rusanov, int option_muscl, int option_solexacte,
                    double xmin, double xmax, double cfl, double tmax,
                    int N,
                    char * option_equation);

void parameters_init_file(parameters *par,
                    char * path_input, char * path_output,
                    int option_godunov, int option_rusanov, int option_muscl, int option_solexacte);

void parameters_plot(parameters *par);

void parameters_free(parameters *par);

void godunov_solve(parameters *par, int option_visual);

void rusanov_solve(parameters *par, int option_visual);

void muscl_solve(parameters *par, int option_visual);

void parameters_error_init(parameters_error *pperr,
                        int option_godunov, int option_rusanov, int option_muscl,
                        double xmin, double xmax, double cfl, double tmax,
                        int len_liste_N, int * liste_N,
                        char * option_error, char * option_equation);

void parameters_error_init_file(parameters_error *pperr,
                        char * path_input, char * path_output,
                        int option_godunov, int option_rusanov, int option_muscl);

void parameters_error_plot(parameters_error *pperr);

void parameters_error_free(parameters_error * pperr);

void parameters_error_compute(parameters_error *pperr);

#endif