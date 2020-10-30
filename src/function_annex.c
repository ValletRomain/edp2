#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>
#include <dirent.h>

#include "raler.c"
#include "godunov.h"
#include "parameters.c"

#define CHEMIN_MAX 512
#define BORDER 0.1

//-----------------------------------------------------------------------------
// Manage of parameters
//-----------------------------------------------------------------------------

void godunov_parameters(godunov * pgd, char * option){

    if (strcmp(option,"transport_1d_1") == 0){
        pgd->pfluxnum = fluxnum_trans1;
        pgd->plambda_ma = lambda_ma_trans1;
        pgd->pboundary_spatial = boundary_spatial_trans1;
        pgd->pboundary_temporal_left = boundary_temporal_left_trans1;
        pgd->pboundary_temporal_right = boundary_temporal_right_trans1;

        if (pgd->keept_solexacte){
            pgd->psolexacte = solexacte_trans1;
        }
    }
    else if (strcmp(option,"burgers_1d_1") == 0){
        pgd->pfluxnum = fluxnum_burgers1;
        pgd->plambda_ma = lambda_ma_burgers1;
        pgd->pboundary_spatial = boundary_spatial_burgers1;
        pgd->pboundary_temporal_left = boundary_temporal_left_burgers1;
        pgd->pboundary_temporal_right = boundary_temporal_right_burgers1;

        if (pgd->keept_solexacte){
            pgd->psolexacte = solexacte_burgers1;
        }
    }
    else if (strcmp(option,"burgers_1d_2") == 0){
        pgd->pfluxnum = fluxnum_burgers2;
        pgd->plambda_ma = lambda_ma_burgers2;
        pgd->pboundary_spatial = boundary_spatial_burgers2;
        pgd->pboundary_temporal_left = boundary_temporal_left_burgers2;
        pgd->pboundary_temporal_right = boundary_temporal_right_burgers2;

        if (pgd->keept_solexacte){
            pgd->psolexacte = solexacte_burgers2;
        }
    }
    else {
        raler(0,"The option \"%s\" dot not exist", option);
    }
    
}

void godunov_error_parameters(godunov_error * pgderr, char * option_error){

    if (option_error = "norm_L1"){
        pgderr->perror = norm_L1;
    }
    else if (option_error = "norm_L2"){
        pgderr->perror = norm_L2;
    }
    else {
        raler(0,"The option \"%s\" dot not exist", option_error);
    }
}


//-----------------------------------------------------------------------------
// Fonction output
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Fonction pour créer les fichiers outputs godunov

void gd_create_parameters(godunov * pgd, char * output_path){

    char * name_file = malloc(CHEMIN_MAX);
    strcpy(name_file, output_path);
    strcat(name_file, "/");
    strcat(name_file, "parameters");
    
    FILE *fic = fopen(name_file, "w");
    
    fprintf(fic, "Parametres %s :\n\n", pgd->name_file);
    fprintf(fic, "keept_solexacte %d\n", pgd->keept_solexacte);
    fprintf(fic, "N %d\n", pgd->N);
    fprintf(fic, "m %d\n", pgd->m);
    fprintf(fic, "dx %f\n", pgd->dx);
    fprintf(fic, "tmax %f\n", pgd->tmax);
    fprintf(fic, "xmin %f\n", pgd->xmin);
    fprintf(fic, "xmax %f\n", pgd->xmax);
    fprintf(fic, "cfl %f\n", pgd->cfl);
    fprintf(fic, "time %ld\n", pgd->time);
    
    fclose(fic);
    free(name_file);
}

void gd_create_plot(godunov * pgd, char * output_path){

    char * name_file = malloc(CHEMIN_MAX);
    strcpy(name_file, output_path);
    strcat(name_file, "/");
    strcat(name_file, "plot.dat");

    FILE *fic = fopen(name_file, "w");

    if (pgd->keept_solexacte){
        for (int i = 0; i < pgd->N+2; i++){
            fprintf(fic, "%f %f %f\n", pgd->xi[i], pgd->un[i], pgd->sol[i]);
        }
    }
    else {
        for (int i = 0; i < pgd->N+2; i++){
            fprintf(fic, "%f %f\n", pgd->xi[i], pgd->un[i]);
        }
    }

    fclose(fic);
    free(name_file);
}

void gd_create_execute_gnu(godunov * pgd, char * output_path){

    // Create of border

    char * name_file = malloc(CHEMIN_MAX);
    strcpy(name_file, output_path);
    strcat(name_file, "/");
    strcat(name_file, "plotcom.gnu");

    FILE *fic = fopen(name_file, "w");
    
    fprintf(fic, "set terminal pngcairo\n");
    fprintf(fic, "set output \'%s/graphe.png\'\n\n", output_path);
    fprintf(fic, "set title \"Resolution de %s tmax=%f\"\n", pgd->option, pgd->tmax);
    fprintf(fic, "set xlabel \"x\"\n");
    fprintf(fic, "set ylabel \"u\"\n\n");
    fprintf(fic, "stats \'%s/plot.dat\' using 1:2 nooutput\n", output_path);
    fprintf(fic, "set xrange [STATS_min_x:STATS_max_x]\n");
    fprintf(fic, "set yrange [STATS_min_y - %f * (STATS_max_y-STATS_min_y): STATS_max_y + %f * (STATS_max_y-STATS_min_y)]\n\n", BORDER, BORDER);
    fprintf(fic, "plot \'%s/plot.dat\' using 1:2 title \"solution numerique\" w lp pt 0", output_path);
    if (pgd->keept_solexacte){
        fprintf(fic, ", \'%s/plot.dat\' using 1:3 title \"soluton exacte\" w lp pt 0", output_path);
    }
    
    

    fclose(fic);
    free(name_file);

    // Execution de la commande gnuplot
    char * name_command = malloc(CHEMIN_MAX);
    strcpy(name_command, "gnuplot ");
    strcat(name_command, output_path);
    strcat(name_command, "/");
    strcat(name_command, "plotcom.gnu");
    
    int status = system(name_command);
    assert(status == EXIT_SUCCESS);

    free(name_command);
}

//-----------------------------------------------------------------------------
// Fonction pour créer les fichiers outputs godunov_error

void gderr_create_parameters(godunov_error * pgderr, char * output_path){

    char * name_output = malloc(CHEMIN_MAX);
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
}

void gderr_create_plot(godunov_error * pgderr, char * output_path){

    char * name_output = malloc(CHEMIN_MAX);
    strcpy(name_output, output_path);
    strcat(name_output, "/");
    strcat(name_output, "plot.dat");

    FILE *fic = fopen(name_output, "w");

    for (int i = 0; i < pgderr->len_liste_N; i++){

        fprintf(fic, "%d %f %ld\n", pgderr->liste_N[i], pgderr->liste_error[i], pgderr->liste_time[i]);

    }

    fclose(fic);
    free(name_output);
}

void gderr_create_execute_gnu(godunov_error * pgderr, char * output_path){

    char * name_file = malloc(CHEMIN_MAX);
    strcpy(name_file, output_path);
    strcat(name_file, "/");
    strcat(name_file, "plotcom.gnu");

    FILE *fic = fopen(name_file, "w");
    
    fprintf(fic, "set terminal pngcairo\n\n");
    fprintf(fic, "# Graphic of error\n");
    fprintf(fic, "set output \'%s/error.png\'\n\n", output_path);
    fprintf(fic, "set title \"Erreur du schéma de Godunov pour %s en %s\"\n", pgderr->option_godunov, pgderr->option_error);
    fprintf(fic, "set xlabel \"N\"\n");
    fprintf(fic, "set ylabel \"error\"\n\n");
    fprintf(fic, "set logscale x 10\n");
    fprintf(fic, "stats \'%s/plot.dat\' using 1:2 nooutput\n", output_path);
    fprintf(fic, "set xrange [STATS_min_x:STATS_max_x]\n");
    fprintf(fic, "set yrange [0: STATS_max_y + %f * (STATS_max_y-STATS_min_y)]\n\n", BORDER);
    fprintf(fic, "plot \'%s/plot.dat\' using 1:2 title \"error\" w lp\n\n", output_path);
    fprintf(fic, "reset\n\n");
    fprintf(fic, "# Graphic of time\n");
    fprintf(fic, "set output \'%s/time.png\'\n\n", output_path);
    fprintf(fic, "set title \"Duree\"\n");
    fprintf(fic, "set xlabel \"N\"\n");
    fprintf(fic, "set ylabel \"time (s)\"\n\n");
    //fprintf(fic, "set logscale x 10\n");
    fprintf(fic, "set autoscale y\n");
    fprintf(fic, "stats \'%s/plot.dat\' using 1:3 nooutput\n", output_path);
    fprintf(fic, "set xrange [STATS_min_x:STATS_max_x]\n");
    //fprintf(fic, "set yrange [0: 1 + STATS_max_y + %f * (STATS_max_y-STATS_min_y)]\n\n", BORDER);
    fprintf(fic, "plot \'%s/plot.dat\' using 1:3 title \"time\" w lp", output_path);
    
    fclose(fic);
    free(name_file);

    // Execution de la commande gnuplot
    char * name_command = malloc(CHEMIN_MAX);
    strcpy(name_command, "gnuplot ");
    strcat(name_command, output_path);
    strcat(name_command, "/");
    strcat(name_command, "plotcom.gnu");
    
    int status = system(name_command);
    assert(status == EXIT_SUCCESS);

    free(name_command);
}