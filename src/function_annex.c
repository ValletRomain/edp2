#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>
#include <dirent.h>

#include "raler.c"
#include "parameters.h"
#include "parameters_equation.c"

#define CHEMIN_MAX 512
#define BORDER 0.1

//-----------------------------------------------------------------------------
// Manage of parameters
//-----------------------------------------------------------------------------

void give_parameters(parameters * ppar, char * option){
    // Distrib the function pointers for parameters ppar depending on option
    // If option do not exist, a alert message is send

    if (strcmp(option,"transport_1d_1") == 0){
        ppar->pfluxnum = fluxnum_trans1;
        ppar->plambda_ma = lambda_ma_trans1;
        ppar->pboundary_spatial = boundary_spatial_trans1;
        ppar->pboundary_temporal_left = boundary_temporal_left_trans1;
        ppar->pboundary_temporal_right = boundary_temporal_right_trans1;

        if (ppar->keept_solexacte){
            ppar->psolexacte = solexacte_trans1;
        }
    }
    else if (strcmp(option,"burgers_1d_1") == 0){
        ppar->pfluxnum = fluxnum_burgers;
        ppar->plambda_ma = lambda_ma_burgers;
        ppar->pboundary_spatial = boundary_spatial_burgers1;
        ppar->pboundary_temporal_left = boundary_temporal_left_burgers1;
        ppar->pboundary_temporal_right = boundary_temporal_right_burgers1;

        if (ppar->keept_solexacte){
            ppar->psolexacte = solexacte_burgers1;
        }
    }
    else if (strcmp(option,"burgers_1d_2") == 0){
        ppar->pfluxnum = fluxnum_burgers;
        ppar->plambda_ma = lambda_ma_burgers;
        ppar->pboundary_spatial = boundary_spatial_burgers2;
        ppar->pboundary_temporal_left = boundary_temporal_left_burgers2;
        ppar->pboundary_temporal_right = boundary_temporal_right_burgers2;

        if (ppar->keept_solexacte){
            ppar->psolexacte = solexacte_burgers2;
        }
    }
    else {
        raler(0,"The option \"%s\" dot not exist", option);
    }
    
}

void give_error_parameters(parameters_error * pparerr, char * option_error){
    // Distrib the function pointers for parameters_error pparerr depending on
    // option
    // If option do not exist, a alert message is send

    if (option_error = "norm_L1"){
        pparerr->perror = norm_L1;
    }
    else if (option_error = "norm_L2"){
        pparerr->perror = norm_L2;
    }
    else if (option_error = "norm_inf")
        pparerr->perror = norm_inf;
    else {
        raler(0,"The option \"%s\" dot not exist", option_error);
    }
}


//-----------------------------------------------------------------------------
// Fonction output
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Fonction pour créer les fichiers outputs godunov

void par_create_parameters(parameters * ppar){
    // Create file parameters for parameters object ppar

    char * name_file = malloc(CHEMIN_MAX);
    strcpy(name_file, ppar->complete_path_output);
    strcat(name_file, "/");
    strcat(name_file, "parameters");
    
    FILE *fic = fopen(name_file, "w");
    
    fprintf(fic, "Parametres %s :\n\n", ppar->name_file);
    fprintf(fic, "keept_solexacte %d\n", ppar->keept_solexacte);
    fprintf(fic, "N %d\n", ppar->N);
    fprintf(fic, "m %d\n", ppar->m);
    fprintf(fic, "dx %f\n", ppar->dx);
    fprintf(fic, "tmax %f\n", ppar->tmax);
    fprintf(fic, "xmin %f\n", ppar->xmin);
    fprintf(fic, "xmax %f\n", ppar->xmax);
    fprintf(fic, "cfl %f\n", ppar->cfl);
    fprintf(fic, "time %ld\n", ppar->time);
    
    fclose(fic);
    free(name_file);
}

void par_create_plot(parameters * ppar){
    // Create file plot.dat for parameters object ppar
    // Give exact solution if ppar->keept_solution=1

    char * name_file = malloc(CHEMIN_MAX);
    strcpy(name_file, ppar->complete_path_output);
    strcat(name_file, "/");
    strcat(name_file, "plot.dat");

    FILE *fic = fopen(name_file, "w");

    if (ppar->keept_solexacte){
        for (int i = 0; i < ppar->N+2; i++){
            fprintf(fic, "%f %f %f\n", ppar->xi[i], ppar->un[i], ppar->sol[i]);
        }
    }
    else {
        for (int i = 0; i < ppar->N+2; i++){
            fprintf(fic, "%f %f\n", ppar->xi[i], ppar->un[i]);
        }
    }

    fclose(fic);
    free(name_file);
}

void par_create_plots(parameters * ppar){
    // Create the ploti.dat with i=ppar->int_tnow

    char * name_plot = malloc(CHEMIN_MAX);
    strcpy(name_plot, ppar->complete_path_output);
    char * inter = malloc(CHEMIN_MAX);
    sprintf(inter, "/plots/plot%d.dat", ppar->int_tnow);
    strcat(name_plot, inter);

    free(inter);

    FILE * fic = fopen(name_plot, "w");

    for (int i=0; i<ppar->N+2; i++){
        fprintf(fic, "%f %f \n", ppar->xi[i], ppar->un[i]);
    }

    fclose(fic);
    free(name_plot);
}

void par_create_execute_gnu(parameters * ppar){
    // Create and execute the gnuplot plotcom.gnu for parameters object ppar

    // Create of plotcom.gnu
    char * name_file = malloc(CHEMIN_MAX);
    strcpy(name_file, ppar->complete_path_output);
    strcat(name_file, "/");
    strcat(name_file, "plotcom.gnu");

    FILE *fic = fopen(name_file, "w");
    
    fprintf(fic, "set terminal pngcairo\n");
    fprintf(fic, "set output \'%s/graphe.png\'\n\n", ppar->complete_path_output);
    fprintf(fic, "set title \"Resolution de %s tmax=%f\"\n", ppar->option_godunov, ppar->tmax);
    fprintf(fic, "set xlabel \"x\"\n");
    fprintf(fic, "set ylabel \"u\"\n\n");
    fprintf(fic, "stats \'%s/plot.dat\' using 1:2 nooutput\n", ppar->complete_path_output);
    fprintf(fic, "set xrange [STATS_min_x:STATS_max_x]\n");
    fprintf(fic, "set yrange [STATS_min_y - %f * (STATS_max_y-STATS_min_y): STATS_max_y + %f * (STATS_max_y-STATS_min_y)]\n\n", BORDER, BORDER);
    fprintf(fic, "plot \'%s/plot.dat\' using 1:2 title \"solution numerique\" w lp pt 0", ppar->complete_path_output);
    if (ppar->keept_solexacte){
        fprintf(fic, ", \'%s/plot.dat\' using 1:3 title \"soluton exacte\" w lp pt 0", ppar->complete_path_output);
    }
    
    fclose(fic);
    free(name_file);

    // Execution de la commande gnuplot
    char * name_command = malloc(CHEMIN_MAX);
    strcpy(name_command, "gnuplot ");
    strcat(name_command, ppar->complete_path_output);
    strcat(name_command, "/");
    strcat(name_command, "plotcom.gnu");
    
    int status = system(name_command);
    assert(status == EXIT_SUCCESS);

    free(name_command);
}

void par_create_animation(parameters * ppar){
    // Create and execute the gnuplot plotcom.gnu for parameters object ppar

    // Create of folder animation
    char * path_animation = malloc(CHEMIN_MAX);
    strcpy(path_animation, ppar->complete_path_output);
    strcat(path_animation, "/animation");

    mkdir(path_animation, ACCESSPERMS);

    free(path_animation);

    // Create of plotcom.gnu
    char * name_file = malloc(CHEMIN_MAX);
    strcpy(name_file, ppar->complete_path_output);
    strcat(name_file, "/plotcom.gnu");

    FILE *fic = fopen(name_file, "w");
    
    fprintf(fic, "set terminal pngcairo\n\n");
    fprintf(fic, "stats \'%s/plots/plot1.dat\' using 1:2 nooutput\n\n", ppar->complete_path_output);
    fprintf(fic, "Xmin = STATS_min_x\n");
    fprintf(fic, "Xmax = STATS_max_x\n");
    fprintf(fic, "Ymin = STATS_min_y - %f * (STATS_max_y - STATS_min_y)\n", BORDER);
    fprintf(fic, "Ymax = STATS_max_y + %f * (STATS_max_y - STATS_min_y)\n\n", BORDER);

    fprintf(fic, "do for [i=0:%d] {\n", ppar->int_tnow);                        // <<--------- PEUT être une erreur ici (int_tnow-1 ou int_tnow ???)
    fprintf(fic, "\tset output sprintf(\'%s/animation/graphe%%d.png\',i) \n\n", ppar->complete_path_output);
    fprintf(fic, "\tset title sprintf(\"Animation de %s int_tnow=%%d\",i) \n", ppar->option_godunov);
    fprintf(fic, "\tset key off\n\n");
    fprintf(fic, "\tset xlabel \"x\" \n");
    fprintf(fic, "\tset ylabel \"u\" \n\n");
    fprintf(fic, "\tset xrange [Xmin:Xmax] \n");
    fprintf(fic, "\tset yrange [Ymin:Ymax] \n\n");
    fprintf(fic, "\tplot sprintf(\'%s/plots/plot%%d.dat\', i) using 1:2 w lp pt 0 \n", ppar->complete_path_output);
    fprintf(fic, "}");
    
    fclose(fic);
    free(name_file);
    
    // Execution de la commande gnuplot
    char * name_command = malloc(CHEMIN_MAX);
    strcpy(name_command, "gnuplot ");
    strcat(name_command, ppar->complete_path_output);
    strcat(name_command, "/plotcom.gnu");
    
    int status = system(name_command);
    assert(status == EXIT_SUCCESS);
    
    free(name_command);

    // Creation of video
    name_command = malloc(CHEMIN_MAX);
    sprintf(name_command,
                "mencoder mf://%s/animation/*.png -mf w=800:h=600:fps=25:type=png -ovc lavc -lavcopts vcodec=mpeg4 -oac copy -o %s/animation.avi",
                ppar->complete_path_output, ppar->complete_path_output);
    
    status = system(name_command);
    assert(status == EXIT_SUCCESS);
    
    free(name_command);

    printf("Fin create of animation\n");
}


//-----------------------------------------------------------------------------
// Fonction pour créer les fichiers outputs parameters_error

void parerr_create_parameters(parameters_error * pparerr, char * output_path){
    // Create parameters for parameters_error object pparerr

    char * name_output = malloc(CHEMIN_MAX);
    strcpy(name_output, output_path);
    strcat(name_output, "/");
    strcat(name_output, "parameters");

    FILE *fic = fopen(name_output, "w");

    fprintf(fic, "Parametres %s:\n\n", output_path);
    fprintf(fic, "len_liste_N %d\n", pparerr->len_liste_N);
    fprintf(fic, "liste_N ");
    for (int i=0; i<pparerr->len_liste_N; i++){
        fprintf(fic, "%d ", pparerr->liste_N[i]);
    }
    fprintf(fic, "\n");
    fprintf(fic, "m %d\n", pparerr->m);
    fprintf(fic, "xmin %f\n", pparerr->xmin);
    fprintf(fic, "xmax %f\n", pparerr->xmax);
    fprintf(fic, "cfl %f\n", pparerr->cfl);

    fclose(fic);
    free(name_output);
}

void parerr_create_plot(parameters_error * pparerr, char * output_path){
    // Create plot.dat for parameters_error object pparerr

    char * name_output = malloc(CHEMIN_MAX);
    strcpy(name_output, output_path);
    strcat(name_output, "/");
    strcat(name_output, "plot.dat");

    FILE *fic = fopen(name_output, "w");

    for (int i = 0; i < pparerr->len_liste_N; i++){

        fprintf(fic, "%d %f %ld\n", pparerr->liste_N[i], pparerr->liste_error[i], pparerr->liste_time[i]);

    }

    fclose(fic);
    free(name_output);
}

void parerr_create_execute_gnu(parameters_error * pparerr, char * output_path){
    // Create and execute the gnuplot script plotcom.gnu for parameters_error object pparerr

    char * name_file = malloc(CHEMIN_MAX);
    strcpy(name_file, output_path);
    strcat(name_file, "/");
    strcat(name_file, "plotcom.gnu");

    FILE *fic = fopen(name_file, "w");
    
    fprintf(fic, "set terminal pngcairo\n\n");
    fprintf(fic, "# Graphic of error\n");
    fprintf(fic, "set output \'%s/error.png\'\n\n", output_path);
    fprintf(fic, "set title \"Erreur du schéma de Godunov pour %s en %s\"\n", pparerr->option_godunov, pparerr->option_error);
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