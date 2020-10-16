#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>

#include "parameters.c"

#define CHEMIN_MAX 512


//-----------------------------------------------------------------------------
// Compute of schema of Godunov
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Input

void godunov_init(godunov *pgd,
                    double xmin, double xmax, double cfl,
                    int m, int N, double tmax,
                    char * option){

    pgd->xmin = xmin;
    pgd->xmax = xmax;
    pgd->m = m;
    pgd->N = N;
    pgd->cfl = cfl;
    pgd->tmax = tmax;
    godunov_parameters(pgd, option);

    pgd->dx = (xmax - xmin) / N;
    pgd->dt = pgd->dx * cfl / pgd->pspeed();

    pgd->xi = malloc((N+2) * sizeof(double) * m);
    pgd->un = malloc((N+2) * sizeof(double) * m);
    pgd->unp1 = malloc((N+2) * sizeof(double) * m);
    pgd->sol = malloc((N+2) * sizeof(double)*m);

    for (int i = 0; i < N + 2; i++){

        pgd->xi[i] = xmin + pgd->dx/2 + (i-1)*pgd->dx;
        pgd->psolexacte(pgd->xi[i], 0, pgd->un + i*m);
        pgd->psolexacte(pgd->xi[i], tmax, pgd->sol + i*m);
    }

    if (0){
        printf("Fin Initialisation godunov\n");
    }
}

void godunov_init_file(godunov *pgd, char * name_input){
    // pgd pointeur de la structure godunov regroupant les différentes variables du problème
    // name_input le nom du fichier regroupant les variables du problème

    // Initialisation des variables

    FILE * file = NULL;
    char * line = malloc(CHEMIN_MAX);
    char * str = malloc(CHEMIN_MAX);
    const char * separators = " ";
    
    //--------------------------------------------------------
    // Lecture du fichier

    file = fopen(name_input, "r");

    // La consigne + le saut de ligne
    fgets(line, CHEMIN_MAX, file);
    printf("Initialisation de %s", line);
    fgets(line, CHEMIN_MAX, file);

    // option :
    fgets(line, CHEMIN_MAX, file);

    str = strtok(line, separators);
    str = strtok(NULL, separators);
    char * option= str;

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
    // Calcul des autres valeurs

    godunov_init(pgd, xmin, xmax, cfl, m, N, tmax, option);

}


//-----------------------------------------------------------------------------
// Output

void godunov_plot(godunov *pgd, char * output_path){
    // ATTENTION, ne prend pas en compte les dimension m>1


    char * name_output;

    // Si le dossier n'existe pas on en refait un                           <<<<---------- A FAIRE
    // Et si un fichier .gnu n'exite pas on en cree un

    // Creation du fichier parameters
    name_output = malloc(CHEMIN_MAX);
    strcpy(name_output, output_path);
    strcat(name_output, "/");
    strcat(name_output, "parameters");

    FILE *fic = fopen(name_output, "w");

    fprintf(fic, "Parametres %s:\n\n", output_path);
    fprintf(fic, "N %d\n", pgd->N);
    fprintf(fic, "m %d\n", pgd->m);
    fprintf(fic, "dt %f\n", pgd->dt);
    fprintf(fic, "dx %f\n", pgd->dx);
    fprintf(fic, "xmin %f\n", pgd->xmin);
    fprintf(fic, "xmax %f\n", pgd->xmax);
    fprintf(fic, "cfl %f\n", pgd->cfl);
    fprintf(fic, "time %ld\n", pgd->time);

    fclose(fic);
    free(name_output);

    // Creation du fichier .dat
    name_output = malloc(CHEMIN_MAX);
    strcpy(name_output, output_path);
    strcat(name_output, "/");
    strcat(name_output, "plot.dat");

    fic = fopen(name_output, "w");

    for (int i = 0; i < pgd->N+2; i++){

        fprintf(fic, "%f %f %f\n", pgd->xi[i], pgd->sol[i], pgd->un[i]);

    }

    fclose(fic);
    free(name_output);
    
    // Creation du fichier plotcom.gnu si il n'existe pas
    name_output = malloc(CHEMIN_MAX);
    strcpy(name_output, output_path);
    strcat(name_output, "/");
    strcat(name_output, "plotcom.gnu");

    fic = fopen(name_output, "a");
    
    fclose(fic);
    free(name_output);

    // Execution de la commande gnuplot
    char * name_command = malloc(CHEMIN_MAX);
    strcpy(name_command, "gnuplot ");
    strcat(name_command, output_path);
    strcat(name_command, "/");
    strcat(name_command, "plotcom.gnu");
    
    int status = system(name_command);
    assert(status == EXIT_SUCCESS);

    printf("Fin Plot godunov\n");
}

//-----------------------------------------------------------------------------
// Free

void godunov_free(godunov *pgd){

    free(pgd->xi);
    free(pgd->un);
    free(pgd->unp1);
    free(pgd->sol);

    printf("Fin Liberation godunov\n");
}

//-----------------------------------------------------------------------------
// Solveur

void godunov_solve(godunov *pgd, int option_visual){
    
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
            double vloc = pgd->plambda_ma(pgd->un + m * i);
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
        int i=0;
        pgd->psolexacte(pgd->xi[i], tnow, pgd->unp1 + i * m);
        i = pgd->N + 1;
        pgd->psolexacte(pgd->xi[i], tnow, pgd->unp1 + i * m);

        memcpy(pgd->un, pgd->unp1, (pgd->N + 2) * m *sizeof(double));
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
                        double xmin, double xmax, double cfl,
                        int m, int len_liste_N, int * liste_N,
                        double tmax, char * option_error, char * option_godunov){

    pgderr->xmin = xmin;
    pgderr->xmax = xmax;
    pgderr->cfl = cfl;
    pgderr->m = m;
    pgderr->len_liste_N = len_liste_N;
    pgderr->liste_N = liste_N;
    pgderr->tmax = tmax;
    pgderr->option_error = option_error;
    pgderr->option_godunov = option_godunov;
    godunov_error_parameters(pgderr, option_error);

    pgderr->liste_error = malloc(len_liste_N * sizeof(double));
    pgderr->liste_time = malloc(len_liste_N * sizeof(unsigned long));
}

void godunov_error_init_file(godunov_error *pgderr, char * name_input){
    // pgd pointeur de la structure godunov regroupant les différentes variables du problème
    // name_input le nom du fichier regroupant les variables du problème

    // Initialisation des variables

    FILE * file = NULL;
    char * line = malloc(CHEMIN_MAX);
    char * str = malloc(CHEMIN_MAX);
    const char * separators = " \n";
    
    //--------------------------------------------------------
    // Lecture du fichier
    file = fopen(name_input, "r");

    // La consigne + le saut de ligne
    fgets(line, CHEMIN_MAX, file);
    printf("Initialisation de %s", line);
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
    // Initialisation
    
    godunov_error_init(pgderr,
                    xmin, xmax, cfl,
                    m, len_liste_N, liste_N,
                    tmax, option_error, option_godunov);

    free(option_error);
    free(option_godunov);
}


//-----------------------------------------------------------------------------
// Ouput

void godunov_error_plot(godunov_error *pgderr, char * output_path){

    char * name_output;

    // Si le dossier n'existe pas on en refait un                           <<<<---------- A FAIRE
    // Et si un fichier .gnu n'exite pas on en cree un

    // Creation du fichier parameters
    name_output = malloc(CHEMIN_MAX);
    strcpy(name_output, output_path);
    strcat(name_output, "/");
    strcat(name_output, "parameters");

    FILE *fic = fopen(name_output, "w");

    fprintf(fic, "Parametres %s:\n\n", output_path);
    fprintf(fic, "len_liste_N %d\n", pgderr->len_liste_N);
    fprintf(fic, "liste_N ");
    for (int i=0; i<pgderr->len_liste_N; i++){
        fprintf(fic, "%d ", pgderr->liste_N[i]);
    }
    fprintf(fic, "\n");
    fprintf(fic, "m %d\n", pgderr->m);
    fprintf(fic, "xmin %f\n", pgderr->xmin);
    fprintf(fic, "xmax %f\n", pgderr->xmax);
    fprintf(fic, "cfl %f\n", pgderr->cfl);

    fclose(fic);
    free(name_output);

    // Creation du fichier error.dat
    name_output = malloc(CHEMIN_MAX);
    strcpy(name_output, output_path);
    strcat(name_output, "/");
    strcat(name_output, "error.dat");

    fic = fopen(name_output, "w");

    for (int i = 0; i < pgderr->len_liste_N; i++){

        fprintf(fic, "%d %f\n", pgderr->liste_N[i], pgderr->liste_error[i]);

    }

    fclose(fic);
    free(name_output);
    
    // Creation du fichier time.dat
    name_output = malloc(CHEMIN_MAX);
    strcpy(name_output, output_path);
    strcat(name_output, "/");
    strcat(name_output, "time.dat");

    fic = fopen(name_output, "w");

    for (int i = 0; i < pgderr->len_liste_N; i++){

        fprintf(fic, "%d %ld\n", pgderr->liste_N[i], pgderr->liste_time[i]);

    }

    fclose(fic);
    free(name_output);

    // Creation du fichier plotcom.gnu si il n'existe pas
    name_output = malloc(CHEMIN_MAX);
    strcpy(name_output, output_path);
    strcat(name_output, "/");
    strcat(name_output, "plotcom.gnu");

    fic = fopen(name_output, "a");
    
    fclose(fic);
    free(name_output);

    // Execution de la commande gnuplot
    char * name_command = malloc(CHEMIN_MAX);
    strcpy(name_command, "gnuplot ");
    strcat(name_command, output_path);
    strcat(name_command, "/");
    strcat(name_command, "plotcom.gnu");
    
    int status = system(name_command);
    assert(status == EXIT_SUCCESS);
    
    printf("Fin Plot godunov_error\n");

}


//-----------------------------------------------------------------------------
// Free

void godunov_error_free(godunov_error * pgd){

    free(pgd->liste_N);
    free(pgd->liste_error);
    free(pgd->liste_time);

    printf("Fin Liberation godunov_error\n");

}


//-----------------------------------------------------------------------------
// Calcul

void godunov_error_compute(godunov_error *pgderr){
    // Calcul l'erreur en norme L1

    godunov gd;
    
    for (int i=0; i<pgderr->len_liste_N; i++){
        
        godunov_init(&gd, pgderr->xmin, pgderr->xmax, pgderr->cfl,
                        pgderr->m, pgderr->liste_N[i],
                        pgderr->tmax, pgderr->option_godunov);

        godunov_solve(&gd, 0);

        pgderr->liste_error[i] = pgderr->perror(gd.N+2, gd.m, gd.un, gd.sol);
        pgderr->liste_time[i] = gd.time;

        printf("Compute error for N=%d error=%f time=%ld\n", gd.N, pgderr->liste_error[i], pgderr->liste_time[i]); // <<<<----- L'ERREUR EST LA !!!!!!!!!!!
        
    }

    printf("Fin calcul des erreurs godunov_error\n");

}