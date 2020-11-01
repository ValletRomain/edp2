#ifndef _PARAMETERS_H
#define _PARAMETERS_H

typedef struct parameters{

    // Name of different file or folder
    char * name_file; //name of file of input

    char * path_input; // path of input file
    char * complete_path_output; // path of output folder (output_path/name_file)

    // Option
    int keept_solexacte; //keept the solution (=1) or not (=0)
    int option_animation;

    // Parametre du probleme
    int m, N; // m dimension of problem, N number of space point
    double xmin, xmax;
    double cfl;
    double dt, dx;
    double tmax;
    char * option_godunov; // equation is resolve (burgers, transport...)

    void (*pfluxnum)(double*, double*, double*);
    double (*plambda_ma)(double*);

    void (*pboundary_spatial)(double, double*);
    void (*pboundary_temporal_left)(double, double, double*);
    void (*pboundary_temporal_right)(double, double, double*);

    void (*psolexacte)(double, double, double*);


    // Resultats du probleme
    unsigned long time;
    
    double *xi; // centre des milieux des cellules
    double *un; // solution a l'instant n
    double *unp1; // solution a l'instant n+1
    double *sol; // solution exact

    int int_tnow;

} parameters;

typedef struct parameters_error{

    // Name of different file or folder
    char * name_file;

    // Parametre du probleme
    int m; // nombre de variables conservatives, nombre de cellules
    double dt, dx; // pas de temps, pas d'espace
    double xmin, xmax; // bornes de l'intervalles
    double cfl; // vmax : dt/dx
    double tmax;
    char * option_error;
    char * option_godunov;

    int len_liste_N;
    int * liste_N;

    double (*perror)(int, double*, double*); // error (norm_L1, L2, inf...)

    double * liste_error;
    unsigned long * liste_time;

} parameters_error;

void parameters_init(parameters *par,
                    int keept_solexacte, int option_animation,
                    double xmin, double xmax, double cfl, double tmax,
                    int m, int N,
                    char * option_godunov);

void parameters_init_file(parameters *par, char * path_input, char * path_output, int option_animation);

void parameters_plot(parameters *par);

void parameters_free(parameters *par);

void parameters_solve(parameters *par, int option_visual);

void parameters_error_init(parameters_error *pperr,
                        char * name_file,
                        double xmin, double xmax, double cfl, double tmax,
                        int m, int len_liste_N, int * liste_N,
                        char * option_error, char * option_godunov);

void parameters_error_init_file(parameters_error *pperr, char * name_input);

void parameters_error_plot(parameters_error *pperr, char * output_path);

void parameters_error_free(parameters_error * pperr);

void parameters_error_compute(parameters_error *pperr);

#endif