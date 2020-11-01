#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>
#include <dirent.h>
#include <sys/dir.h>
#include <sys/stat.h>

#include "function_annex.c"

#define CHEMIN_MAX 512


//-----------------------------------------------------------------------------
// Compute of schema of Godunov
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Input

void parameters_init(parameters *ppar,
                    int keept_solexacte, int option_animation,
                    double xmin, double xmax, double cfl, double tmax,
                    int m, int N,
                    char * option_godunov){
    // Initialize the object pointed by ppar with the arguments
    
    ppar->keept_solexacte = keept_solexacte;
    ppar->option_animation = option_animation;

    ppar->xmin = xmin;
    ppar->xmax = xmax;
    ppar->m = m;
    ppar->N = N;
    ppar->cfl = cfl;
    ppar->tmax = tmax;
    give_parameters(ppar, option_godunov);
    ppar->option_godunov = malloc(CHEMIN_MAX);
    strcpy(ppar->option_godunov, option_godunov);

    ppar->dx = (xmax - xmin) / N;

    ppar->xi = malloc((N+2) * sizeof(double) * m);
    ppar->un = malloc((N+2) * sizeof(double) * m);
    ppar->unp1 = malloc((N+2) * sizeof(double) * m);
    ppar->sol = malloc((N+2) * sizeof(double)*m);

    for (int i = 0; i < N + 2; i++){

        ppar->xi[i] = xmin + ppar->dx/2 + (i-1)*ppar->dx;
        ppar->pboundary_spatial(ppar->xi[i], ppar->un + i*m);

        if (keept_solexacte){
            ppar->psolexacte(ppar->xi[i], tmax, ppar->sol + i*m);
        }
    }
}

void parameters_init_file(parameters *ppar, char * path_input, char * path_output, int option_animation){
    // Initialize the object pointed by ppar with the file of path name_input

    FILE * file = NULL;
    char * line = malloc(CHEMIN_MAX);
    char * str = malloc(CHEMIN_MAX);
    const char * separators = " \n";
    const char * separators1 = "/";

    //--------------------------------------------------------
    // Sauvegarde du nom du fichier
    
    char * name_file = malloc(CHEMIN_MAX);
    strcpy(name_file, path_input);
    strcpy(str, path_input);

    strtok(str, separators1); 
    while ( (str=strtok(NULL, separators1)) != NULL){
        strcpy(name_file, str);
    }

    //--------------------------------------------------------
    // Lecture du fichier

    file = fopen(path_input, "r");

    // La consigne + le saut de ligne
    fgets(line, CHEMIN_MAX, file);
    printf("Initialisation <- %s\n", path_input);
    fgets(line, CHEMIN_MAX, file);

    // keept_solexacte :
    fgets(line, CHEMIN_MAX, file);

    str = strtok(line, separators);
    str = strtok(NULL, separators);
    int keept_solexacte = atoi(str);

    // option :
    fgets(line, CHEMIN_MAX, file);

    str = strtok(line, separators);
    str = strtok(NULL, separators);
    char * option = malloc(CHEMIN_MAX);
    strcpy(option, str);

    // xmin :
    fgets(line, CHEMIN_MAX, file);

    str = strtok(line, separators);
    str = strtok(NULL, separators);
    double xmin = atof(str);

    // xmax :
    fgets(line, CHEMIN_MAX, file);

    str = strtok(line, separators);
    str = strtok(NULL, separators);
    double xmax = atof(str);

    // cfl :
    fgets(line, CHEMIN_MAX, file);

    str = strtok(line, separators);
    str = strtok(NULL, separators);
    double cfl = atof(str);

    // m :
    fgets(line, CHEMIN_MAX, file);

    str = strtok(line, separators);
    str = strtok(NULL, separators);
    int m = atoi(str);

    // N :
    fgets(line, CHEMIN_MAX, file);

    str = strtok(line, separators);
    str = strtok(NULL, separators);
    int N = atoi(str);

    // tmax :
    fgets(line, CHEMIN_MAX, file);

    str = strtok(line, separators);
    str = strtok(NULL, separators);
    double tmax = atof(str);

    fclose(file);

    //--------------------------------------------------------
    // Initialization of ppar
    
    // Save of path_input
    ppar->path_input = malloc(CHEMIN_MAX);
    strcpy(ppar->path_input, path_input);
    
    // Save of name_input
    ppar->name_file = malloc(CHEMIN_MAX);
    strcpy(ppar->name_file, name_file);
    
    // Creation of path of output (ppar->complete_output_path)
    ppar->complete_path_output = malloc(CHEMIN_MAX);
    strcpy(ppar->complete_path_output, path_output);
    strcat(ppar->complete_path_output, ppar->name_file);

    // Initialization of other parameters
    parameters_init(ppar, keept_solexacte, option_animation, xmin, xmax, cfl, tmax, m, N, option);

    //--------------------------------------------------------
    // Creation of folder

    // Creation of folder output
    mkdir(ppar->complete_path_output, ACCESSPERMS);
    
    // Creation of parameters
    par_create_parameters(ppar);

    if (option_animation){
        ppar->int_tnow = 0;
        
        // Creation of folder plots in output
        char * path_plots = malloc(CHEMIN_MAX);
        strcpy(path_plots, ppar->complete_path_output);
        strcat(path_plots, "/plots");

        mkdir(path_plots, ACCESSPERMS);

        free(path_plots);

        // Creation of plots0.dat
        par_create_plots(ppar);
    }

    printf("Creation of folfer -> %s\n", ppar->complete_path_output);

    free(name_file);
    free(line);

    printf("Fin Initilization\n");
}


//-----------------------------------------------------------------------------
// Output

void parameters_plot(parameters *ppar){
    // Create and fill the folder of output (of path out_path)
    // with the result of ppar
    
    // Creation of file plot.dat
    par_create_plot(ppar);
    
    // Creation and execution of file plotcom.gnu
    par_create_execute_gnu(ppar);

    // Creation of animation
    if (ppar->option_animation){
        par_create_animation(ppar);
    }

    printf("Fin Plot\n");
}


//-----------------------------------------------------------------------------
// Free

void parameters_free(parameters *ppar){
    // Liberate tables the of ppar

    free(ppar->xi);
    free(ppar->un);
    free(ppar->unp1);
    free(ppar->sol);

    printf("Fin Liberation parameters\n");
}


//-----------------------------------------------------------------------------
// Solveur

void parameters_solve(parameters *ppar, int option_visual){
    // Solve the problem of ppar
    // option_visual give visuality on terminal
    
    int m = ppar->m;

    if (option_visual){
        printf("Debut Resolution\n");
    }

    time_t begin = time(NULL);
    
    double tnow = 0;
    while(tnow < ppar->tmax){
        
        // calcul de la vitesse max
        double vmax = 0;
        for (int i = 0; i < ppar->N + 2; i++){
            double vloc = fabs(ppar->plambda_ma(ppar->un + m * i));
            vmax = vmax > vloc ? vmax : vloc;
        }
        
        ppar->dt = ppar->cfl * ppar->dx / vmax;
        for(int i = 1; i < ppar->N+1; i++){
            double flux[m];
            ppar->pfluxnum(ppar->un + i*m, ppar->un + (i+1)*m, flux);
            for(int iv = 0; iv < m; iv++){
                ppar->unp1[i*m + iv] = ppar->un[i*m + iv] - ppar->dt/ppar->dx * flux[iv];
            }
            ppar->pfluxnum(ppar->un + (i - 1) * m, ppar->un + i * m, flux);
            for(int iv = 0; iv < m;iv++){
                ppar->unp1[i * m + iv] += ppar->dt / ppar->dx * flux[iv];
            }
        }
        // mise à jour
        tnow += ppar->dt;

        if (option_visual){
            printf("tnow = %f vmax = %f tmax = %f\n", tnow, vmax, ppar->tmax);
        }
        
        // conditions aux limites
        ppar->pboundary_temporal_left(ppar->xmin, tnow, ppar->unp1);
        ppar->pboundary_temporal_right(ppar->xmax, tnow, ppar->unp1 + (ppar->N+1)*m);

        memcpy(ppar->un, ppar->unp1, (ppar->N + 2) * m *sizeof(double));

        if (ppar->option_animation){
            ppar->int_tnow++;
            par_create_plots(ppar);
        }
    }
    time_t end = time(NULL);

    ppar->time = (unsigned long) difftime(end, begin);

    if (option_visual){
        printf("Fin Resolution\n");
    }
}


//-----------------------------------------------------------------------------
// Calcul des erreurs
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Input

void parameters_error_init(parameters_error *pparerr,
                        char * name_file,
                        double xmin, double xmax, double cfl, double tmax,
                        int m, int len_liste_N, int * liste_N,
                        char * option_error, char * option_godunov){
    // Initialize the parameters_error with the arguments
    
    pparerr->name_file = name_file;

    pparerr->xmin = xmin;
    pparerr->xmax = xmax;
    pparerr->cfl = cfl;
    pparerr->m = m;
    pparerr->len_liste_N = len_liste_N;
    pparerr->liste_N = liste_N;
    pparerr->tmax = tmax;

    pparerr->option_error = malloc(CHEMIN_MAX);
    pparerr->option_godunov = malloc(CHEMIN_MAX);
    strcpy(pparerr->option_error, option_error);
    strcpy(pparerr->option_godunov, option_godunov);

    give_error_parameters(pparerr, option_error);

    pparerr->liste_error = malloc(len_liste_N * sizeof(double));
    pparerr->liste_time = malloc(len_liste_N * sizeof(unsigned long));

    printf("Fin Initialisation\n");
}

void parameters_error_init_file(parameters_error *pparerr, char * name_input){
    // Initialize of parameters_error pparerr with file of path name_input

    FILE * file = NULL;
    char * line = malloc(CHEMIN_MAX);
    char * str = malloc(CHEMIN_MAX);
    const char * separators = " \n";
    const char * separators1 = "/";
    
    //--------------------------------------------------------
    // Sauvegarde du nom du fichier
    
    char * name_file = malloc(CHEMIN_MAX);
    strcpy(name_file, name_input);
    strcpy(str, name_input);

    strtok(str, separators1); 
    while ( (str=strtok(NULL, separators1)) != NULL){
        strcpy(name_file, str);
    }

    //--------------------------------------------------------
    // Lecture du fichier

    file = fopen(name_input, "r");

    // La consigne + le saut de ligne
    fgets(line, CHEMIN_MAX, file);
    printf("Initialisation <- %s\n", name_input);
    fgets(line, CHEMIN_MAX, file);

    // option_error :
    fgets(line, CHEMIN_MAX, file);

    str = strtok(line, separators);
    str = strtok(NULL, separators);
    char * option_error = malloc(CHEMIN_MAX);
    strcpy(option_error, str);

    // option_godunov :
    fgets(line, CHEMIN_MAX, file);

    str = strtok(line, separators);
    str = strtok(NULL, separators);
    char * option_godunov = malloc(CHEMIN_MAX);
    strcpy(option_godunov, str);

    // xmin :
    fgets(line, CHEMIN_MAX, file);

    str = strtok(line, separators);
    str = strtok(NULL, separators);
    double xmin = atof(str);

    // xmax :
    fgets(line, CHEMIN_MAX, file);

    str = strtok(line, separators);
    str = strtok(NULL, separators);
    double xmax = atof(str);

    // cfl :
    fgets(line, CHEMIN_MAX, file);

    str = strtok(line, separators);
    str = strtok(NULL, separators);
    double cfl = atof(str);

    // m :
    fgets(line, CHEMIN_MAX, file);

    str = strtok(line, separators);
    str = strtok(NULL, separators);
    int m = atoi(str);

    // len_liste_N :
    fgets(line, CHEMIN_MAX, file);

    str = strtok(line, separators);
    str = strtok(NULL, separators);
    int len_liste_N = atoi(str);

    // liste_N :
    int * liste_N = malloc(len_liste_N * sizeof(int));

    fgets(line, CHEMIN_MAX, file);
    str = strtok(line, separators);
    for (int i=0; i<len_liste_N; i++){
        str = strtok(NULL, separators);
        liste_N[i] = atoi(str);
    }

    // tmax :
    fgets(line, CHEMIN_MAX, file);

    str = strtok(line, separators);
    str = strtok(NULL, separators);
    double tmax = atof(str);

    fclose(file);

    //--------------------------------------------------------
    // Initialisation of pparerr 
    
    parameters_error_init(pparerr, name_file,
                    xmin, xmax, cfl, tmax,
                    m, len_liste_N, liste_N,
                    option_error, option_godunov);

    free(option_error);
    free(option_godunov);
}


//-----------------------------------------------------------------------------
// Ouput

void parameters_error_plot(parameters_error *pparerr, char * output_path){
    // Creation of folder (of path output_path) of result of pparerr

    // Creation du dossier
    char * output_path_final = malloc(CHEMIN_MAX);
    strcpy(output_path_final, output_path);
    strcat(output_path_final, pparerr->name_file);
    
    mkdir(output_path_final, ACCESSPERMS);

    // Creation du fichier parameters
    parerr_create_parameters(pparerr, output_path_final);

    // Creation du fichier plot.dat
    parerr_create_plot(pparerr, output_path_final);

    // Creation du fichier plotcom.gnu si il n'existe pas
    parerr_create_execute_gnu(pparerr, output_path_final);

    printf("Fin Plot -> %s\n", output_path_final);

    free(output_path_final);

}


//-----------------------------------------------------------------------------
// Free

void parameters_error_free(parameters_error * pgd){
    // Liberate the tables of pgd

    free(pgd->liste_N);
    free(pgd->liste_error);
    free(pgd->liste_time);

    printf("Fin Liberation\n");
}


//-----------------------------------------------------------------------------
// Calcul

void parameters_error_compute(parameters_error *pparerr){
    // Compute the error between the numeric solution and exact solution of
    // same problem with different N (give by pparerr->liste_N)

    parameters par;
    
    for (int i=0; i<pparerr->len_liste_N; i++){

        parameters_init(&par, 1, 0,
                        pparerr->xmin, pparerr->xmax, pparerr->cfl, pparerr->tmax,
                        pparerr->m, pparerr->liste_N[i],
                        pparerr->option_godunov);

        parameters_solve(&par, 0);

        pparerr->liste_error[i] = pparerr->perror((par.N+2)*par.m, par.un, par.sol);
        pparerr->liste_time[i] = par.time;

        printf("Compute error for N=%d error=%f time=%ld s\n", par.N, pparerr->liste_error[i], pparerr->liste_time[i]);   
    }

    printf("Fin calcul des erreurs\n");

}