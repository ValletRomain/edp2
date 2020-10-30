#ifndef _GODUNOV_H
#define _GODUNOV_H

typedef struct godunov{

    // Name of different file or folder
    char * name_file; //name of file of input

    // Option
    int keept_solexacte; //keept the solution (=1) or not (=0)

    // Parametre du probleme
    int m, N; // m dimension of problem, N number of space point
    double xmin, xmax;
    double cfl;
    double dt, dx;
    double tmax;
    char * option; // equation is resolve (burgers, transport...)

    void (*pfluxnum)(double*, double*, double, double*);
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

} godunov;

typedef struct godunov_error{

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

    double (*perror)(int, int, double*, double*); // error (norm_L1, L2, inf...)

    double * liste_error;
    unsigned long * liste_time;

} godunov_error;

void godunov_init(godunov *pgd,
                    char * name_file, int keept_solution,
                    double xmin, double xmax, double cfl, double tmax,
                    int m, int N,
                    char * option,
                    int option_visual);

void godunov_init_file(godunov *pgd, char * name_input);

void godunov_plot(godunov *pgd, char * output_path);

void godunov_free(godunov *pgd);

void godunov_solve(godunov *pgd, int option_visual);

void godunov_error_init(godunov_error *pgderr,
                        char * name_file,
                        double xmin, double xmax, double cfl, double tmax,
                        int m, int len_liste_N, int * liste_N,
                        char * option_error, char * option_godunov);

void godunov_error_init_file(godunov_error *pgderr, char * name_input);

void godunov_error_plot(godunov_error *pgderr, char * output_path);

void godunov_error_free(godunov_error * pgd);

void godunov_error_compute(godunov_error *pgderr);

#endif