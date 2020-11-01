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

void godunov_init(godunov *pgd,
                    int keept_solexacte, int option_animation,
                    double xmin, double xmax, double cfl, double tmax,
                    int m, int N,
                    char * option_godunov){
    // Initialize the object pointed by pgd with the arguments
    
    pgd->keept_solexacte = keept_solexacte;
    pgd->option_animation = option_animation;

    pgd->xmin = xmin;
    pgd->xmax = xmax;
    pgd->m = m;
    pgd->N = N;
    pgd->cfl = cfl;
    pgd->tmax = tmax;
    godunov_parameters(pgd, option_godunov);
    pgd->option_godunov = malloc(CHEMIN_MAX);
    strcpy(pgd->option_godunov, option_godunov);

    pgd->dx = (xmax - xmin) / N;

    pgd->xi = malloc((N+2) * sizeof(double) * m);
    pgd->un = malloc((N+2) * sizeof(double) * m);
    pgd->unp1 = malloc((N+2) * sizeof(double) * m);
    pgd->sol = malloc((N+2) * sizeof(double)*m);

    for (int i = 0; i < N + 2; i++){

        pgd->xi[i] = xmin + pgd->dx/2 + (i-1)*pgd->dx;
        pgd->pboundary_spatial(pgd->xi[i], pgd->un + i*m);

        if (keept_solexacte){
            pgd->psolexacte(pgd->xi[i], tmax, pgd->sol + i*m);
        }
    }
}

void godunov_init_file(godunov *pgd, char * path_input, char * path_output, int option_animation){
    // Initialize the object pointed by pgd with the file of path name_input

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
    // Initialization of pgd
    
    // Save of path_input
    pgd->path_input = malloc(CHEMIN_MAX);
    strcpy(pgd->path_input, path_input);
    
    // Save of name_input
    pgd->name_file = malloc(CHEMIN_MAX);
    strcpy(pgd->name_file, name_file);
    
    // Creation of path of output (pgd->complete_output_path)
    pgd->complete_path_output = malloc(CHEMIN_MAX);
    strcpy(pgd->complete_path_output, path_output);
    strcat(pgd->complete_path_output, pgd->name_file);

    // Initialization of other parameters
    godunov_init(pgd, keept_solexacte, option_animation, xmin, xmax, cfl, tmax, m, N, option);

    //--------------------------------------------------------
    // Creation of folder

    // Creation of folder output
    mkdir(pgd->complete_path_output, ACCESSPERMS);
    
    // Creation of parameters
    gd_create_parameters(pgd);

    if (option_animation){
        pgd->len_U = 0;
        
        // Creation of folder plots in output
        char * path_plots = malloc(CHEMIN_MAX);
        strcpy(path_plots, pgd->complete_path_output);
        strcat(path_plots, "/plots");

        mkdir(path_plots, ACCESSPERMS);

        free(path_plots);

        // Creation of plots0.dat
        gd_create_plots(pgd);
    }

    printf("Creation of folfer -> %s\n", pgd->complete_path_output);

    free(name_file);
    free(line);

    printf("Fin Initilization\n");
}


//-----------------------------------------------------------------------------
// Output

void godunov_plot(godunov *pgd){
    // Create and fill the folder of output (of path out_path)
    // with the result of pgd
    
    // Creation of file plot.dat
    gd_create_plot(pgd);
    
    // Creation and execution of file plotcom.gnu
    gd_create_execute_gnu(pgd);

    // Creation of animation
    if (pgd->option_animation){
        gd_create_animation(pgd);
    }

    printf("Fin Plot godunov\n");
}


//-----------------------------------------------------------------------------
// Free

void godunov_free(godunov *pgd){
    // Liberate tables the of pgd

    free(pgd->xi);
    free(pgd->un);
    free(pgd->unp1);
    free(pgd->sol);

    printf("Fin Liberation godunov\n");
}


//-----------------------------------------------------------------------------
// Solveur

void godunov_solve(godunov *pgd, int option_visual){
    // Solve the godunov problem of pgd
    // option_visual give visuality on terminal
    
    int m = pgd->m;

    if (option_visual){
        printf("Debut Resolution\n");
    }

    time_t begin = time(NULL);
    
    double tnow = 0;
    while(tnow < pgd->tmax){
        
        // calcul de la vitesse max
        double vmax = 0;
        for (int i = 0; i < pgd->N + 2; i++){
            double vloc = fabs(pgd->plambda_ma(pgd->un + m * i));
            vmax = vmax > vloc ? vmax : vloc;
        }
        
        pgd->dt = pgd->cfl * pgd->dx / vmax;
        for(int i = 1; i < pgd->N+1; i++){
            double flux[m];
            pgd->pfluxnum(pgd->un + i*m, pgd->un + (i+1)*m, flux);
            for(int iv = 0; iv < m; iv++){
                pgd->unp1[i*m + iv] = pgd->un[i*m + iv] - pgd->dt/pgd->dx * flux[iv];
            }
            pgd->pfluxnum(pgd->un + (i - 1) * m, pgd->un + i * m, flux);
            for(int iv = 0; iv < m;iv++){
                pgd->unp1[i * m + iv] += pgd->dt / pgd->dx * flux[iv];
            }
        }
        // mise à jour
        tnow += pgd->dt;

        if (option_visual){
            printf("tnow = %f vmax = %f tmax = %f\n", tnow, vmax, pgd->tmax);
        }
        
        // conditions aux limites
        pgd->pboundary_temporal_left(pgd->xmin, tnow, pgd->unp1);
        pgd->pboundary_temporal_right(pgd->xmax, tnow, pgd->unp1 + (pgd->N+1)*m);

        memcpy(pgd->un, pgd->unp1, (pgd->N + 2) * m *sizeof(double));

        if (pgd->option_animation){
            pgd->len_U++;
            gd_create_plots(pgd);
        }
    }
    time_t end = time(NULL);

    pgd->time = (unsigned long) difftime(end, begin);

    if (option_visual){
        printf("Fin Resolution godunov\n");
    }
}


//-----------------------------------------------------------------------------
// Calcul des erreurs
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Input

void godunov_error_init(godunov_error *pgderr,
                        char * name_file,
                        double xmin, double xmax, double cfl, double tmax,
                        int m, int len_liste_N, int * liste_N,
                        char * option_error, char * option_godunov){
    // Initialize the godunov_error with the arguments
    
    pgderr->name_file = name_file;

    pgderr->xmin = xmin;
    pgderr->xmax = xmax;
    pgderr->cfl = cfl;
    pgderr->m = m;
    pgderr->len_liste_N = len_liste_N;
    pgderr->liste_N = liste_N;
    pgderr->tmax = tmax;

    pgderr->option_error = malloc(CHEMIN_MAX);
    pgderr->option_godunov = malloc(CHEMIN_MAX);
    strcpy(pgderr->option_error, option_error);
    strcpy(pgderr->option_godunov, option_godunov);

    godunov_error_parameters(pgderr, option_error);

    pgderr->liste_error = malloc(len_liste_N * sizeof(double));
    pgderr->liste_time = malloc(len_liste_N * sizeof(unsigned long));

    printf("Fin Initialisation\n");
}

void godunov_error_init_file(godunov_error *pgderr, char * name_input){
    // Initialize of godunov_error pgderr with file of path name_input

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
    // Initialisation of pgderr 
    
    godunov_error_init(pgderr, name_file,
                    xmin, xmax, cfl, tmax,
                    m, len_liste_N, liste_N,
                    option_error, option_godunov);

    free(option_error);
    free(option_godunov);
}


//-----------------------------------------------------------------------------
// Ouput

void godunov_error_plot(godunov_error *pgderr, char * output_path){
    // Creation of folder (of path output_path) of result of pgderr

    // Creation du dossier
    char * output_path_final = malloc(CHEMIN_MAX);
    strcpy(output_path_final, output_path);
    strcat(output_path_final, pgderr->name_file);
    
    mkdir(output_path_final, ACCESSPERMS);

    // Creation du fichier parameters
    gderr_create_parameters(pgderr, output_path_final);

    // Creation du fichier plot.dat
    gderr_create_plot(pgderr, output_path_final);

    // Creation du fichier plotcom.gnu si il n'existe pas
    gderr_create_execute_gnu(pgderr, output_path_final);

    printf("Fin Plot godunov_error -> %s\n", output_path_final);

    free(output_path_final);

}


//-----------------------------------------------------------------------------
// Free

void godunov_error_free(godunov_error * pgd){
    // Liberate the tables of pgd

    free(pgd->liste_N);
    free(pgd->liste_error);
    free(pgd->liste_time);

    printf("Fin Liberation godunov_error\n");

}


//-----------------------------------------------------------------------------
// Calcul

void godunov_error_compute(godunov_error *pgderr){
    // Compute the error between the numeric solution and exact solution of
    // same problem with different N (give by pgderr->liste_N)

    godunov gd;
    
    for (int i=0; i<pgderr->len_liste_N; i++){

        godunov_init(&gd, 1, 0,
                        pgderr->xmin, pgderr->xmax, pgderr->cfl, pgderr->tmax,
                        pgderr->m, pgderr->liste_N[i],
                        pgderr->option_godunov);

        godunov_solve(&gd, 0);

        pgderr->liste_error[i] = pgderr->perror((gd.N+2)*gd.m, gd.un, gd.sol);
        pgderr->liste_time[i] = gd.time;

        printf("Compute error for N=%d error=%f time=%ld s\n", gd.N, pgderr->liste_error[i], pgderr->liste_time[i]);   
    }

    printf("Fin calcul des erreurs godunov_error\n");

}