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
        ppar->plambda_ma = lambda_ma_trans1;

        ppar->pfluxnum_gd = fluxnum_gd_trans1;

        ppar->pfluxnum_ru = fluxnum_ru_trans1;

        ppar->pboundary_spatial = boundary_spatial_trans1;
        ppar->pboundary_temporal_left = boundary_temporal_left_trans1;
        ppar->pboundary_temporal_right = boundary_temporal_right_trans1;

        if (ppar->option_solexacte){
            ppar->psolexacte = solexacte_trans1;
        }
    }
    else if (strcmp(option,"burgers_1d_1") == 0){
        ppar->plambda_ma = lambda_ma_burgers;

        ppar->pfluxnum_gd = fluxnum_gd_burgers;

        ppar->pfluxnum_ru = fluxnum_ru_burgers;

        ppar->pboundary_spatial = boundary_spatial_burgers1;
        ppar->pboundary_temporal_left = boundary_temporal_left_burgers1;
        ppar->pboundary_temporal_right = boundary_temporal_right_burgers1;

        if (ppar->option_solexacte){
            ppar->psolexacte = solexacte_burgers1;
        }
    }
    else if (strcmp(option,"burgers_1d_2") == 0){
        ppar->plambda_ma = lambda_ma_burgers;

        ppar->pfluxnum_gd = fluxnum_gd_burgers;
        
        ppar->pfluxnum_ru = fluxnum_ru_burgers;

        ppar->pboundary_spatial = boundary_spatial_burgers2;
        ppar->pboundary_temporal_left = boundary_temporal_left_burgers2;
        ppar->pboundary_temporal_right = boundary_temporal_right_burgers2;

        if (ppar->option_solexacte){
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
    strcat(name_file, "/parameters");
    
    FILE *fic = fopen(name_file, "w");
    
    fprintf(fic, "Parametres %s :\n\n", ppar->name_file);
    fprintf(fic, "option_solexacte %d\n", ppar->option_solexacte);
    fprintf(fic, "N %d\n", ppar->N);
    fprintf(fic, "dx %f\n", ppar->dx);
    fprintf(fic, "tmax %f\n", ppar->tmax);
    fprintf(fic, "xmin %f\n", ppar->xmin);
    fprintf(fic, "xmax %f\n", ppar->xmax);
    fprintf(fic, "cfl %f\n", ppar->cfl);
    if (ppar->option_godunov)
        fprintf(fic, "time godunov %ld\n", ppar->time_gd);
    if (ppar->option_rusanov)
        fprintf(fic, "time rusanov %ld\n", ppar->time_ru);
    if (ppar->option_muscl)
        fprintf(fic, "time MUSCL %ld\n", ppar->time_muscl);
    
    fclose(fic);
    free(name_file);
}

void par_create_plot(parameters * ppar){
    // Create file plot.dat for parameters object ppar
    // Give exact solution if ppar->keept_solution=1

    char * name_file = malloc(CHEMIN_MAX);
    strcpy(name_file, ppar->complete_path_output);
    strcat(name_file, "/plot.dat");

    FILE *fic = fopen(name_file, "w");

    fprintf(fic, "x ");
    if (ppar->option_godunov)
        fprintf(fic, "gd ");
    if (ppar->option_rusanov)
        fprintf(fic, "ru ");
    if (ppar->option_muscl)
        fprintf(fic, "muscl ");
    if (ppar->option_solexacte)
        fprintf(fic, "sol ");
    fprintf(fic, "\n");

    for (int i=0; i<ppar->N+2; i++){
        fprintf(fic, "%f ", ppar->xi[i]);
        if (ppar->option_godunov)
            fprintf(fic, "%f ", ppar->un[i]);
        if (ppar->option_rusanov)
            fprintf(fic, "%f ", ppar->vn[i]);
        if (ppar->option_muscl)
            fprintf(fic, "%f ", ppar->wn[i]);
        if (ppar->option_solexacte)
            fprintf(fic, "%f ", ppar->sol[i]);
        fprintf(fic, "\n");
    }

    fclose(fic);
    free(name_file);
}

void par_create_execute_gnu(parameters * ppar){
    // Create and execute the gnuplot plotcom.gnu for parameters object ppar

    int place_gd = 2;
    int place_ru = ppar->option_godunov + 2;
    int place_sol = ppar->option_godunov + ppar->option_rusanov + 2;

    // Create of plotcom.gnu
    char * name_file = malloc(CHEMIN_MAX);
    strcpy(name_file, ppar->complete_path_output);
    strcat(name_file, "/plotcom.gnu");

    FILE *fic = fopen(name_file, "w");
    
    fprintf(fic, "set terminal pngcairo\n");
    fprintf(fic, "set output \'%s/graphe.png\'\n\n", ppar->complete_path_output);
    fprintf(fic, "set title \"Resolution de %s tmax=%f\"\n", ppar->option_equation, ppar->tmax);
    fprintf(fic, "set xlabel \"x\"\n");
    fprintf(fic, "set ylabel \"u\"\n\n");
    fprintf(fic, "stats \'%s/plot.dat\' using 1:2 nooutput\n", ppar->complete_path_output);
    fprintf(fic, "set xrange [STATS_min_x:STATS_max_x]\n");
    fprintf(fic, "set yrange [STATS_min_y - %f * (STATS_max_y-STATS_min_y): STATS_max_y + %f * (STATS_max_y-STATS_min_y)]\n\n", BORDER, BORDER);

    fprintf(fic, "plot ");

    if (ppar->option_godunov){
        fprintf(fic, "\'%s/plot.dat\' using 1:2 title \"godunov\" w lp pt 0 lc rgb \"blue\"", ppar->complete_path_output);
        
        if ((ppar->option_rusanov) || (ppar->option_muscl) || (ppar->option_solexacte))
            fprintf(fic, ", ");
    }

    if (ppar->option_rusanov){
        fprintf(fic, "\'%s/plot.dat\' using 1:%d title \"rusanov\" w lp pt 0 lc rgb \"red\"",
                    ppar->complete_path_output,
                    ppar->option_godunov + 2);
    
        if ((ppar->option_muscl) || (ppar->option_solexacte))
            fprintf(fic, ", ");
    }

    if (ppar->option_muscl){
        fprintf(fic, "\'%s/plot.dat\' using 1:%d title \"MUSCL\" w lp pt 0 lc rgb \"green\"",
                    ppar->complete_path_output,
                    ppar->option_godunov + ppar->option_rusanov + 2);
    
        if (ppar->option_solexacte)
            fprintf(fic, ", ");
    }

    if (ppar->option_solexacte){
        fprintf(fic, "\'%s/plot.dat\' using 1:%d title \"exacte\" w lp pt 0 lc rgb \"black\"",
                    ppar->complete_path_output,
                    ppar->option_godunov + ppar->option_rusanov + ppar->option_muscl + 2);
    }

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
}

//-----------------------------------------------------------------------------
// Fonction pour créer les fichiers outputs parameters_error

void parerr_create_parameters(parameters_error * pparerr){
    // Create parameters for parameters_error object pparerr

    char * name_output = malloc(CHEMIN_MAX);
    strcpy(name_output, pparerr->complete_path_output);
    strcat(name_output, "/parameters");

    FILE *fic = fopen(name_output, "w");

    fprintf(fic, "Parametres %s:\n\n", pparerr->complete_path_output);
    fprintf(fic, "len_liste_N %d\n", pparerr->len_liste_N);
    fprintf(fic, "liste_N ");
    for (int i=0; i<pparerr->len_liste_N; i++){
        fprintf(fic, "%d ", pparerr->liste_N[i]);
    }
    fprintf(fic, "\n");
    fprintf(fic, "xmin %f\n", pparerr->xmin);
    fprintf(fic, "xmax %f\n", pparerr->xmax);
    fprintf(fic, "cfl %f\n", pparerr->cfl);

    fclose(fic);
    free(name_output);
}

void parerr_create_plot(parameters_error * pparerr){
    // Create plot.dat for parameters_error object pparerr

    char * name_output = malloc(CHEMIN_MAX);
    strcpy(name_output, pparerr->complete_path_output);
    strcat(name_output, "/plot.dat");

    FILE *fic = fopen(name_output, "w");

    fprintf(fic, "N ");

    if (pparerr->option_godunov)
        fprintf(fic, "err_g time_g ");
        
    if (pparerr->option_rusanov)
        fprintf(fic, "err_r time_r ");

    if (pparerr->option_muscl)
        fprintf(fic, "err_m time_m" );

    fprintf(fic, "\n");

    for (int i = 0; i < pparerr->len_liste_N; i++){
        fprintf(fic, "%d ", pparerr->liste_N[i]);
        if (pparerr->option_godunov)
            fprintf(fic, "%f %ld ", pparerr->liste_error_gd[i], pparerr->liste_time_gd[i]);
        
        if (pparerr->option_rusanov)
            fprintf(fic, "%f %ld ", pparerr->liste_error_ru[i], pparerr->liste_time_ru[i]);
        
        if (pparerr->option_muscl)
            fprintf(fic, "%f %ld ", pparerr->liste_error_muscl[i], pparerr->liste_time_muscl[i]);
        
        fprintf(fic, "\n");
    }

    fclose(fic);
    free(name_output);
}

void parerr_create_execute_gnu(parameters_error * pparerr){
    // Create and execute the gnuplot script plotcom.gnu for parameters_error object pparerr

    char * name_file = malloc(CHEMIN_MAX);
    strcpy(name_file, pparerr->complete_path_output);
    strcat(name_file, "/plotcom.gnu");

    FILE *fic = fopen(name_file, "w");
    
    fprintf(fic, "set terminal pngcairo\n\n");

    fprintf(fic, "# Graphic of error\n");
    fprintf(fic, "set output \'%s/error.png\'\n\n", pparerr->complete_path_output);
    fprintf(fic, "set title \"Erreur pour %s en %s\"\n", pparerr->option_equation, pparerr->option_error);
    fprintf(fic, "set xlabel \"N\"\n");
    fprintf(fic, "set ylabel \"error\"\n\n");
    fprintf(fic, "set logscale x 10\n");
    fprintf(fic, "stats \'%s/plot.dat\' using 1:2 nooutput\n", pparerr->complete_path_output);
    fprintf(fic, "set xrange [STATS_min_x:STATS_max_x]\n");
    fprintf(fic, "set yrange [0: STATS_max_y + %f * (STATS_max_y-STATS_min_y)]\n\n", BORDER);

    fprintf(fic, "plot ");
    if (pparerr->option_godunov){
        fprintf(fic, "\'%s/plot.dat\' using 1:2 title \"godunov\" w lp lc rgb \"blue\"", pparerr->complete_path_output);

        if ((pparerr->option_rusanov) || (pparerr->option_muscl))
            fprintf(fic, ", ");
    }
    if (pparerr->option_rusanov){
        fprintf(fic, "\'%s/plot.dat\' using 1:%d title \"rusanov\" w lp lc rgb \"red\"",
                pparerr->complete_path_output,
                2 + 2 * pparerr->option_godunov);
        
        if (pparerr->option_muscl)
            fprintf(fic, ", ");
    }
    if (pparerr->option_muscl)
        fprintf(fic, "\'%s/plot.dat\' using 1:%d title \"MUSCL\" w lp lc rgb \"green\"",
                pparerr->complete_path_output,
                2 + 2 * (pparerr->option_godunov + pparerr->option_rusanov));
    
    fprintf(fic, "\n\n");

    fprintf(fic, "reset\n\n");
    
    fprintf(fic, "# Graphic of time\n");
    fprintf(fic, "set output \'%s/time.png\'\n\n", pparerr->complete_path_output);
    fprintf(fic, "set title \"Duree\"\n");
    fprintf(fic, "set xlabel \"N\"\n");
    fprintf(fic, "set ylabel \"time (s)\"\n\n");
    //fprintf(fic, "set logscale x 10\n");
    fprintf(fic, "set autoscale y\n");
    fprintf(fic, "stats \'%s/plot.dat\' using 1:3 nooutput\n", pparerr->complete_path_output);
    fprintf(fic, "set xrange [STATS_min_x:STATS_max_x]\n");
    //fprintf(fic, "set yrange [0: 1 + STATS_max_y + %f * (STATS_max_y-STATS_min_y)]\n\n", BORDER);
    
    fprintf(fic, "plot ");
    if (pparerr->option_godunov){
        fprintf(fic, "\'%s/plot.dat\' using 1:3 title \"godunov\" w lp lc rgb \"blue\"", pparerr->complete_path_output);

        if ((pparerr->option_rusanov) || (pparerr->option_muscl))
            fprintf(fic, ", ");
    }
    if (pparerr->option_rusanov){
        fprintf(fic, "\'%s/plot.dat\' using 1:%d title \"rusanov\" w lp lc rgb \"red\"",
                pparerr->complete_path_output,
                3 + 2 * pparerr->option_godunov);
        
        if (pparerr->option_muscl)
            fprintf(fic, ", ");
    }
    if (pparerr->option_muscl)
        fprintf(fic, "\'%s/plot.dat\' using 1:%d title \"MUSCL\" w lp lc rgb \"green\"",
                pparerr->complete_path_output,
                3 + 2 * (pparerr->option_godunov + pparerr->option_rusanov));
    
    fclose(fic);
    free(name_file);

    // Execution de la commande gnuplot
    char * name_command = malloc(CHEMIN_MAX);
    strcpy(name_command, "gnuplot ");
    strcat(name_command, pparerr->complete_path_output);
    strcat(name_command, "/plotcom.gnu");
    
    int status = system(name_command);
    assert(status == EXIT_SUCCESS);

    free(name_command);
}