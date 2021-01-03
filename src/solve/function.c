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
                    int option_godunov, int option_rusanov, int option_muscl,int option_solexacte,
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
    
    int m = ppar->m;

    if (option_godunov){
        ppar->un = malloc((N+2) * m * sizeof(double));
        ppar->unp1 = malloc((N+2) * m  * sizeof(double));
    }
    if (option_rusanov){
        ppar->vn = malloc((N+2) * m * sizeof(double));
        ppar->vnp1 = malloc((N+2) * m * sizeof(double));
    }
    if (option_muscl){
        ppar->wn = malloc((N+2) * m * sizeof(double));
        ppar->wnp1 = malloc((N+2) * m * sizeof(double));     
    }
    if (option_solexacte){
        ppar->sol = malloc((N+2) * m * sizeof(double));
    }

    for (int i = 0; i < N + 2; i++){
        
        ppar->xi[i] = xmin + ppar->dx/2 + (i-1)*ppar->dx;
        
        if (option_godunov)
            ppar->pboundary_spatial(ppar->xi + i, ppar->un + i*m);
        
        if (option_rusanov)
            ppar->pboundary_spatial(ppar->xi + i, ppar->vn + i*m);

        if (option_muscl)
            ppar->pboundary_spatial(ppar->xi + i, ppar->wn + i*m);

        if (option_solexacte)
            ppar->pboundary_spatial(ppar->xi + i, ppar->sol + i*m);
        
    }
}

void parameters_init_file(parameters *ppar, char * path_input, char * path_output, int option_godunov, int option_rusanov, int option_muscl, int option_solexacte){
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
                    option_godunov, option_rusanov, option_muscl, option_solexacte,
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

    int m = ppar->m;

    time_t begin = time(NULL);
    
    double tnow = 0;
    while(tnow < ppar->tmax){
        
        // calcul de la vitesse max
        double vmax = 0;
        for (int i = 0; i < ppar->N + 2; i++){
            double vloc = fabs(ppar->plambda_ma(ppar->un + i*m));
            vmax = vmax > vloc ? vmax : vloc;
        }
        
        ppar->dt = ppar->cfl * ppar->dx / vmax;
        for(int i = 1; i < ppar->N+1; i++){
            double flux[ppar->m];
            ppar->pfluxnum_gd(ppar->un + i*m, ppar->un + (i+1)*m, flux);
            for (int iv=0; iv < ppar->m; iv++)
                ppar->unp1[i*m + iv] = ppar->un[i*m + iv] - ppar->dt/ppar->dx * flux[iv];

            ppar->pfluxnum_gd(ppar->un + (i-1)*m, ppar->un + i*m, flux);
            for (int iv=0; iv < ppar->m; iv++)
                ppar->unp1[i*m + iv] += ppar->dt / ppar->dx * flux[iv];
        }
        // mise à jour
        tnow += ppar->dt;

        if (option_visual){
            printf("tnow = %f vmax = %f tmax = %f\n", tnow, vmax, ppar->tmax);
        }
        
        // conditions aux limites
        ppar->pboundary_temporal_left(ppar->xmin, tnow, ppar->unp1);
        ppar->pboundary_temporal_right(ppar->xmax, tnow, ppar->unp1 + ppar->N+1 * m);

        memcpy(ppar->un, ppar->unp1, (ppar->N + 2) * m *sizeof(double));

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
    
    int m = ppar->m;

    time_t begin = time(NULL);
    
    double tnow = 0;
    while(tnow < ppar->tmax){
        
        // calcul de la vitesse max
        double vmax = 0;
        for (int i = 0; i < ppar->N + 2; i++){
            double vloc = fabs(ppar->plambda_ma(ppar->vn + i*m));
            vmax = vmax > vloc ? vmax : vloc;
        }
        
        ppar->dt = ppar->cfl * ppar->dx / vmax; 
        for(int i = 1; i < ppar->N+1; i++){
            double flux[ppar->m];
            ppar->pfluxnum_ru(ppar->vn + i*m, ppar->vn + (i+1)*m, flux);
            for (int iv=0; iv < ppar->m; iv++)
                ppar->vnp1[i*m + iv] = ppar->vn[i*m + iv] - ppar->dt/ppar->dx * flux[iv];
            
            ppar->pfluxnum_ru(ppar->vn + (i-1)*m, ppar->vn + i*m, flux);
            for (int iv=0; iv < ppar->m; iv++)
                ppar->vnp1[i*m + iv] += ppar->dt / ppar->dx * flux[iv];
        } 
        // mise à jour
        tnow += ppar->dt;

        if (option_visual){
            printf("tnow = %f vmax = %f tmax = %f\n", tnow, vmax, ppar->tmax);
        }
        
        // conditions aux limites
        ppar->pboundary_temporal_left(ppar->xmin, tnow, ppar->vnp1);
        ppar->pboundary_temporal_right(ppar->xmax, tnow, ppar->vnp1 + ppar->N+1 * m);

        memcpy(ppar->vn, ppar->vnp1, (ppar->N + 2) * m * sizeof(double));

    }
    time_t end = time(NULL);

    ppar->time_ru = (unsigned long) difftime(end, begin);

    if (option_visual)
        printf("Fin Resolution rusanov\n");
}

void w_to_hu(parameters *ppar){
    // Convert w=(h hu) to (h u)

    int m = ppar->m;

    for (int i=0; i < ppar->N+2; i++){

        if (ppar->option_godunov){
            ppar->un[i*m + 1] = ppar->un[i*m+0] / ppar->un[i*m+1];
        }
        if (ppar->option_rusanov){
            ppar->vn[i*m + 1] = ppar->vn[i*m+0] / ppar->vn[i*m+1];
        }
        if (ppar->option_muscl){
            ppar->wn[i*m + 1] = ppar->wn[i*m+0] / ppar->wn[i*m+1];
        }
    }
}