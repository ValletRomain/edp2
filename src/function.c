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
                    int option_solexacte, int option_godunov, int option_rusanov, int option_muscl,
                    double xmin, double xmax, double cfl, double tmax,
                    int N,
                    char * option_equation){
    // Initialize the object pointed by ppar with the arguments

    ppar->option_solexacte = option_solexacte;
    ppar->option_godunov = option_godunov;
    ppar->option_rusanov = option_rusanov;
    ppar->option_muscl = option_muscl;

    ppar->xmin = xmin;
    ppar->xmax = xmax;
    ppar->N = N;
    ppar->cfl = cfl;
    ppar->tmax = tmax;
    give_parameters(ppar, option_equation);
    ppar->option_equation = malloc(CHEMIN_MAX);
    strcpy(ppar->option_equation, option_equation);

    ppar->dx = (xmax - xmin) / N;

    ppar->xi = malloc((N+2) * sizeof(double));

    if (option_godunov){
        ppar->un = malloc((N+2) * sizeof(double));
        ppar->unp1 = malloc((N+2) * sizeof(double));
    }
    if (option_rusanov){
        ppar->vn = malloc((N+2) * sizeof(double));
        ppar->vnp1 = malloc((N+2) * sizeof(double));
    }
    if (option_muscl){
        ppar->wn = malloc((N+2) * sizeof(double));
        ppar->wnp1 = malloc((N+2) * sizeof(double));        
    }
    if (option_solexacte)
        ppar->sol = malloc((N+2) * sizeof(double));

    for (int i = 0; i < N + 2; i++){

        ppar->xi[i] = xmin + ppar->dx/2 + (i-1)*ppar->dx;

        if (option_godunov)
            ppar->un[i] = ppar->pboundary_spatial(ppar->xi[i]);

        if (option_rusanov)
            ppar->vn[i] = ppar->pboundary_spatial(ppar->xi[i]);

        if (option_muscl)
            ppar->wn[i] = ppar->pboundary_spatial(ppar->xi[i]);

        if (option_solexacte)
            ppar->sol[i] = ppar->psolexacte(ppar->xi[i], tmax);
        
    }
}

void parameters_init_file(parameters *ppar, char * path_input, char * path_output, int option_godunov, int option_rusanov, int option_muscl){
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

    // option_solexacte :
    fgets(line, CHEMIN_MAX, file);

    str = strtok(line, separators);
    str = strtok(NULL, separators);
    int option_solexacte = atoi(str);

    // option_equation :
    fgets(line, CHEMIN_MAX, file);

    str = strtok(line, separators);
    str = strtok(NULL, separators);
    char * option_equation = malloc(CHEMIN_MAX);
    strcpy(option_equation, str);

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
    parameters_init(ppar,
                    option_solexacte, option_godunov, option_rusanov, option_muscl,
                    xmin, xmax, cfl, tmax,
                    N, option_equation);

    //--------------------------------------------------------
    // Creation of folder

    // Creation of folder output (before remove the folder is the folder exist)
    if (access(ppar->complete_path_output, F_OK) == 0)
        remove_directory(ppar->complete_path_output);
    mkdir(ppar->complete_path_output, ACCESSPERMS);
    
    // Creation of parameters
    par_create_parameters(ppar);

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
    
    /*
    // Creation of animation
    if (ppar->option_animation){
        if (ppar->option_godunov)
            par_create_animation(ppar, 0);
        if (ppar->option_godunov)
            par_create_animation(ppar, 1);
    }
    */

    printf("Fin Plot\n");
}


//-----------------------------------------------------------------------------
// Free

void parameters_free(parameters *ppar){
    // Liberate tables the of ppar

    free(ppar->xi);
    if (ppar->option_godunov){
        free(ppar->un);
        free(ppar->unp1);
    }
    if (ppar->option_rusanov){
        free(ppar->vn);
        free(ppar->vnp1);
    }
    if (ppar->option_muscl){
        free(ppar->wn);
        free(ppar->wnp1);
    }
    if (ppar->option_solexacte)
        free(ppar->sol);

    printf("Fin Liberation parameters\n");
}


//-----------------------------------------------------------------------------
// Solveur

void godunov_solve(parameters *ppar, int option_visual){
    // Solve the problem of ppar
    // option_visual give visuality on terminal

    if (option_visual){
        printf("Debut Resolution godunov\n");
    }

    time_t begin = time(NULL);
    
    double tnow = 0;
    while(tnow < ppar->tmax){
        
        // calcul de la vitesse max
        double vmax = 0;
        for (int i = 0; i < ppar->N + 2; i++){
            double vloc = fabs(ppar->plambda_ma(ppar->un[i]));
            vmax = vmax > vloc ? vmax : vloc;
        }
        
        ppar->dt = ppar->cfl * ppar->dx / vmax;
        for(int i = 1; i < ppar->N+1; i++){
            double flux;
            flux = ppar->pfluxnum_gd(ppar->un[i], ppar->un[i+1]);
            ppar->unp1[i] = ppar->un[i] - ppar->dt/ppar->dx * flux;

            flux = ppar->pfluxnum_gd(ppar->un[i-1], ppar->un[i]);
            ppar->unp1[i] += ppar->dt / ppar->dx * flux;
        }
        // mise à jour
        tnow += ppar->dt;

        if (option_visual){
            printf("tnow = %f vmax = %f tmax = %f\n", tnow, vmax, ppar->tmax);
        }
        
        // conditions aux limites
        ppar->unp1[0] = ppar->pboundary_temporal_left(ppar->xmin, tnow);
        ppar->unp1[ppar->N+1] = ppar->pboundary_temporal_right(ppar->xmax, tnow);

        memcpy(ppar->un, ppar->unp1, (ppar->N + 2) *sizeof(double));

    }
    time_t end = time(NULL);

    ppar->time_gd = (unsigned long) difftime(end, begin);

    if (option_visual){
        printf("Fin Resolution godunov\n");
    }
}

void rusanov_solve(parameters *ppar, int option_visual){
    // Solve the problem of ppar
    // option_visual give visuality on terminal

    if (option_visual){
        printf("Debut Resolution rusanov\n");
    }

    time_t begin = time(NULL);
    
    double tnow = 0;
    while(tnow < ppar->tmax){
        
        // calcul de la vitesse max
        double vmax = 0;
        for (int i = 0; i < ppar->N + 2; i++){
            double vloc = fabs(ppar->plambda_ma(ppar->vn[i]));
            vmax = vmax > vloc ? vmax : vloc;
        }
        
        ppar->dt = ppar->cfl * ppar->dx / vmax;
        for(int i = 1; i < ppar->N+1; i++){
            double flux;
            flux = ppar->pfluxnum_ru(ppar->vn[i], ppar->vn[i+1]);
            ppar->vnp1[i] = ppar->vn[i] - ppar->dt/ppar->dx * flux;
            
            flux = ppar->pfluxnum_ru(ppar->vn[i-1], ppar->vn[i]);
            ppar->vnp1[i] += ppar->dt / ppar->dx * flux;
        }
        // mise à jour
        tnow += ppar->dt;

        if (option_visual){
            printf("tnow = %f vmax = %f tmax = %f\n", tnow, vmax, ppar->tmax);
        }
        
        // conditions aux limites
        ppar->vnp1[0] = ppar->pboundary_temporal_left(ppar->xmin, tnow);
        ppar->vnp1[ppar->N+1] = ppar->pboundary_temporal_right(ppar->xmax, tnow);

        memcpy(ppar->vn, ppar->vnp1, (ppar->N + 2) * sizeof(double));

    }
    time_t end = time(NULL);

    ppar->time_ru = (unsigned long) difftime(end, begin);

    if (option_visual)
        printf("Fin Resolution rusanov\n");
}

void muscl_solve(parameters *ppar, int option_visual){
    
    // Solve the problem of ppar
    // option_visual give visuality on terminal

    if (option_visual){
        printf("Debut Resolution MUSCL\n");
    }

    time_t begin = time(NULL);
    
    double tnow = 0;
    while(tnow < ppar->tmax){
        
        // calcul de la vitesse max
        double vmax = 0;
        for (int i = 0; i < ppar->N + 2; i++){
            double vloc = fabs(ppar->plambda_ma(ppar->wn[i]));
            vmax = vmax > vloc ? vmax : vloc;
        }
        
        ppar->dt = ppar->cfl * ppar->dx / vmax;

        for(int i = 1; i < ppar->N+1; i++){
            double flux;
            flux = ppar->pfluxnum_gd(w_half_l(ppar->wn[i-1], ppar->wn[i], ppar->wn[i+1]),
                                        w_half_p(ppar->wn[i], ppar->wn[i+1], ppar->wn[i+2]));
            ppar->wnp1[i] = ppar->wn[i] - ppar->dt/ppar->dx * flux;
            
            flux = ppar->pfluxnum_gd(w_half_l(ppar->wn[i-2], ppar->wn[i-1], ppar->wn[i]),
                                        w_half_p(ppar->wn[i-1], ppar->wn[i], ppar->wn[i+1]));
            ppar->wnp1[i] += ppar->dt / ppar->dx * flux;
        }
        // mise à jour
        tnow += ppar->dt;

        if (option_visual){
            printf("tnow = %f vmax = %f tmax = %f\n", tnow, vmax, ppar->tmax);
        }
        
        // conditions aux limites
        ppar->wnp1[0] = ppar->pboundary_temporal_left(ppar->xmin, tnow);
        ppar->wnp1[ppar->N+1] = ppar->pboundary_temporal_right(ppar->xmax, tnow);

        memcpy(ppar->wn, ppar->wnp1, (ppar->N + 2) * sizeof(double));

    }
    time_t end = time(NULL);

    ppar->time_muscl = (unsigned long) difftime(end, begin);

    if (option_visual){
        printf("Fin Resolution MUSCL\n");
    }

}


//-----------------------------------------------------------------------------
// Calcul des erreurs
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Input

void parameters_error_init(parameters_error *pparerr,
                        int option_godunov, int option_rusanov, int option_muscl,
                        double xmin, double xmax, double cfl, double tmax,
                        int len_liste_N, int * liste_N,
                        char * option_error, char * option_equation){
    // Initialize the parameters_error with the arguments

    pparerr->option_godunov = option_godunov;
    pparerr->option_rusanov = option_rusanov;
    pparerr->option_muscl = option_muscl;

    pparerr->xmin = xmin;
    pparerr->xmax = xmax;
    pparerr->cfl = cfl;
    pparerr->len_liste_N = len_liste_N;
    pparerr->liste_N = liste_N;
    pparerr->tmax = tmax;

    pparerr->option_error = malloc(CHEMIN_MAX);
    pparerr->option_equation = malloc(CHEMIN_MAX);
    strcpy(pparerr->option_error, option_error);
    strcpy(pparerr->option_equation, option_equation);

    give_error_parameters(pparerr, option_error);

    if (option_godunov){
        pparerr->liste_error_gd = malloc(len_liste_N * sizeof(double));
        pparerr->liste_time_gd = malloc(len_liste_N * sizeof(unsigned long));
    }

    if (option_rusanov){
        pparerr->liste_error_ru = malloc(len_liste_N * sizeof(double));
        pparerr->liste_time_ru = malloc(len_liste_N * sizeof(unsigned long));
    }

    if (option_muscl){
        pparerr->liste_error_muscl = malloc(len_liste_N * sizeof(double));
        pparerr->liste_time_muscl = malloc(len_liste_N * sizeof(unsigned long));
    }

    printf("Fin Initialisation\n");
}

void parameters_error_init_file(parameters_error *pparerr,
                                char * path_input, char * path_output,
                                int option_godunov, int option_rusanov, int option_muscl){
    // Initialize of parameters_error pparerr with file of path path_input

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

    // option_error :
    fgets(line, CHEMIN_MAX, file);

    str = strtok(line, separators);
    str = strtok(NULL, separators);
    char * option_error = malloc(CHEMIN_MAX);
    strcpy(option_error, str);

    // option_equation :
    fgets(line, CHEMIN_MAX, file);

    str = strtok(line, separators);
    str = strtok(NULL, separators);
    char * option_equation = malloc(CHEMIN_MAX);
    strcpy(option_equation, str);

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
    // Initialization of ppar
    
    // Save of path_input
    pparerr->path_input = malloc(CHEMIN_MAX);
    strcpy(pparerr->path_input, path_input);
    
    // Save of name_input
    pparerr->name_file = malloc(CHEMIN_MAX);
    strcpy(pparerr->name_file, name_file);
    
    // Creation of path of output (ppar->complete_output_path)
    pparerr->complete_path_output = malloc(CHEMIN_MAX);
    strcpy(pparerr->complete_path_output, path_output);
    strcat(pparerr->complete_path_output, pparerr->name_file);

    // Initialization of other parameters
    parameters_error_init(pparerr,
                            option_godunov, option_rusanov, option_muscl,
                            xmin, xmax, cfl, tmax,
                            len_liste_N, liste_N,
                            option_error, option_equation);

    //--------------------------------------------------------
    // Creation of folder

    // Creation of folder output (before remove the folder is the folder exist)
    if (access(pparerr->complete_path_output, F_OK) == 0)
        remove_directory(pparerr->complete_path_output);
    mkdir(pparerr->complete_path_output, ACCESSPERMS);
    
    // Creation of parameters
    parerr_create_parameters(pparerr);

    printf("Creation of directory -> %s\n", pparerr->complete_path_output);

    printf("Fin Initialisation\n");
}


//-----------------------------------------------------------------------------
// Ouput

void parameters_error_plot(parameters_error *pparerr){
    // Creation of folder (of path output_path) of result of pparerr

    // Creation du fichier plot.dat
    parerr_create_plot(pparerr);

    // Creation du fichier plotcom.gnu si il n'existe pas
    parerr_create_execute_gnu(pparerr);

    printf("Fin Plot -> %s\n", pparerr->complete_path_output);
}


//-----------------------------------------------------------------------------
// Free

void parameters_error_free(parameters_error * pgd){
    // Liberate the tables of pgd

    free(pgd->liste_N);
    //free(pgd->liste_error_gd);
    //free(pgd->liste_time_ru);

    printf("Fin Liberation\n");
}


//-----------------------------------------------------------------------------
// Calcul

void parameters_error_compute(parameters_error *pparerr){
    // Compute the error between the numeric solution and exact solution of
    // same problem with different N (give by pparerr->liste_N)

    parameters par;
    
    for (int i=0; i<pparerr->len_liste_N; i++){

        parameters_init(&par,
                        1, pparerr->option_godunov, pparerr->option_rusanov, pparerr->option_muscl,
                        pparerr->xmin, pparerr->xmax, pparerr->cfl, pparerr->tmax,
                        pparerr->liste_N[i],
                        pparerr->option_equation);

        if (pparerr->option_godunov){
            godunov_solve(&par, 0);
            pparerr->liste_error_gd[i] = pparerr->perror((par.N+2), par.un, par.sol);
            pparerr->liste_time_gd[i] = par.time_gd;
            printf("Compute error godunov for N=%d error=%f time=%ld s\n", par.N, pparerr->liste_error_gd[i], pparerr->liste_time_gd[i]);   
        }
        
        if (pparerr->option_rusanov){
            rusanov_solve(&par, 0);
            pparerr->liste_error_ru[i] = pparerr->perror((par.N+2), par.vn, par.sol);
            pparerr->liste_time_ru[i] = par.time_ru;
            printf("Compute error rusanov for N=%d error=%f time=%ld s\n", par.N, pparerr->liste_error_ru[i], pparerr->liste_time_ru[i]);   
        } 

        if (pparerr->option_muscl){
            muscl_solve(&par, 0);
            pparerr->liste_error_muscl[i] = pparerr->perror((par.N+2), par.vn, par.sol);
            pparerr->liste_time_muscl[i] = par.time_ru;
            printf("Compute error MUSCL for N=%d error=%f time=%ld s\n", par.N, pparerr->liste_error_muscl[i], pparerr->liste_time_muscl[i]);   
        } 
    }

    printf("Fin calcul des erreurs\n");

}