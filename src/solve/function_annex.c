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
// Fonction init
//-----------------------------------------------------------------------------

void give_parameters(parameters * ppar, char * option_equation){

    if (strcmp(option_equation, "saint_venant_1d") == 0){

        ppar->m = 2;

        ppar->plambda_ma = lambda_ma_sv;
        ppar->pfluxnum_gd = fluxnum_gd_sv;
        ppar->pfluxnum_ru = fluxnum_ru_sv;

        ppar->pboundary_spatial = boundary_spatial_1;
        ppar->pboundary_temporal_left = boundary_temporal_left_1;
        ppar->pboundary_temporal_right = boundary_temporal_right_1;

        //ppar->psolexacte = solexacte_sv;
    }
    else if (strcmp(option_equation, "bassin") == 0){

        ppar->m = 2;

        ppar->plambda_ma = lambda_ma_sv;
        ppar->pfluxnum_gd = fluxnum_gd_sv;
        ppar->pfluxnum_ru = fluxnum_ru_sv;

        ppar->pboundary_spatial = boundary_spatial_bassin;
        ppar->pboundary_temporal_left = boundary_temporal_left_bassin;
        ppar->pboundary_temporal_right = boundary_temporal_right_bassin;
    }
    else {
        raler(0,"The option \"%s\" dot not exist", option_equation);
    }
}

/*
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
*/

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
    fprintf(fic, "m %d\n", ppar->m);
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

    int m = ppar->m;

    char * name_file = malloc(CHEMIN_MAX);
    strcpy(name_file, ppar->complete_path_output);
    strcat(name_file, "/plot.dat");

    FILE *fic = fopen(name_file, "w");

    fprintf(fic, "x ");
    if (ppar->option_godunov){
        for (int iv=0; iv<ppar->m; iv++)
            fprintf(fic, "gd_%d ", iv);
    }
    if (ppar->option_rusanov){
        for (int iv=0; iv<ppar->m; iv++)
            fprintf(fic, "ru_%d ", iv);
    }
    if (ppar->option_muscl){
        for (int iv=0; iv<ppar->m; iv++)
            fprintf(fic, "muscl_%d ", iv);
    }
    if (ppar->option_solexacte){
        for (int iv=0; iv<ppar->m; iv++)
            fprintf(fic, "sol_%d ", iv);
    }
    fprintf(fic, "\n");

    for (int i=0; i<ppar->N+2; i++){

        fprintf(fic, "%f ", ppar->xi[i]);

        if (ppar->option_godunov){
            for (int iv=0; iv<ppar->m; iv++)
                fprintf(fic, "%f ", ppar->un[i*m + iv]);
        }
        if (ppar->option_rusanov){
            for (int iv=0; iv<ppar->m; iv++)
                fprintf(fic, "%f ", ppar->vn[i*m + iv]);
        }
        if (ppar->option_muscl){
            for (int iv=0; iv<ppar->m; iv++)
                fprintf(fic, "%f ", ppar->wn[i*m + iv]);
        }
        if (ppar->option_solexacte){
            for (int iv=0; iv<ppar->m; iv++)
                fprintf(fic, "%f ", ppar->sol[i*m + iv]);
        }
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

    fprintf(fic, "\n");
    fprintf(fic, "set multiplot layout %d, 1 title \"Resolution de Saint Venant tmax=%f\"", ppar->m, ppar->tmax);
    fprintf(fic, "\n");

    fprintf(fic, "set xlabel \"x\"\n");
    fprintf(fic, "set ylabel \"h\"\n\n");
    //fprintf(fic, "stats \'%s/plot.dat\' using 1:2 nooutput\n", ppar->complete_path_output);
    //fprintf(fic, "set xrange [STATS_min_x:STATS_max_x]\n");
    //fprintf(fic, "set yrange [STATS_min_y - %f * (STATS_max_y-STATS_min_y): STATS_max_y + %f * (STATS_max_y-STATS_min_y)]\n\n", BORDER, BORDER);

    fprintf(fic, "plot ");

    if (ppar->option_godunov){
        fprintf(fic, "\'%s/plot.dat\' using 1:2 title \"godunov\" w lp pt 0 lc rgb \"blue\"", ppar->complete_path_output);
        
        if ((ppar->option_rusanov) || (ppar->option_muscl) || (ppar->option_solexacte))
            fprintf(fic, ", ");
    }

    if (ppar->option_rusanov){
        fprintf(fic, "\'%s/plot.dat\' using 1:%d title \"rusanov\" w lp pt 0 lc rgb \"red\"",
                    ppar->complete_path_output,
                    2 * ppar->option_godunov + 2);
    
        if ((ppar->option_muscl) || (ppar->option_solexacte))
            fprintf(fic, ", ");
    }

    if (ppar->option_muscl){
        fprintf(fic, "\'%s/plot.dat\' using 1:%d title \"MUSCL\" w lp pt 0 lc rgb \"green\"",
                    ppar->complete_path_output,
                    2*(ppar->option_godunov+ppar->option_rusanov) + 2);
    
        if (ppar->option_solexacte)
            fprintf(fic, ", ");
    }

    if (ppar->option_solexacte){
        fprintf(fic, "\'%s/plot.dat\' using 1:%d title \"exacte\" w lp pt 0 lc rgb \"black\"",
                    ppar->complete_path_output,
                    2*(ppar->option_godunov + ppar->option_rusanov + ppar->option_muscl) + 2);
    }

    fprintf(fic, "\n\n");

    fprintf(fic, "set xlabel \"x\"\n");
    fprintf(fic, "set ylabel \"u\"\n\n");
    //fprintf(fic, "stats \'%s/plot.dat\' using 1:2 nooutput\n", ppar->complete_path_output);
    //fprintf(fic, "set xrange [STATS_min_x:STATS_max_x]\n");
    //fprintf(fic, "set yrange [STATS_min_y - %f * (STATS_max_y-STATS_min_y): STATS_max_y + %f * (STATS_max_y-STATS_min_y)]\n\n", BORDER, BORDER);

    fprintf(fic, "plot ");

    if (ppar->option_godunov){
        fprintf(fic, "\'%s/plot.dat\' using 1:3 title \"godunov\" w lp pt 0 lc rgb \"blue\"", ppar->complete_path_output);
        
        if ((ppar->option_rusanov) || (ppar->option_muscl) || (ppar->option_solexacte))
            fprintf(fic, ", ");
    }

    if (ppar->option_rusanov){
        fprintf(fic, "\'%s/plot.dat\' using 1:%d title \"rusanov\" w lp pt 0 lc rgb \"red\"",
                    ppar->complete_path_output,
                    2 * ppar->option_godunov + 3);
    
        if ((ppar->option_muscl) || (ppar->option_solexacte))
            fprintf(fic, ", ");
    }

    if (ppar->option_muscl){
        fprintf(fic, "\'%s/plot.dat\' using 1:%d title \"MUSCL\" w lp pt 0 lc rgb \"green\"",
                    ppar->complete_path_output,
                    2*(ppar->option_godunov+ppar->option_rusanov) + 3);
    
        if (ppar->option_solexacte)
            fprintf(fic, ", ");
    }

    if (ppar->option_solexacte){
        fprintf(fic, "\'%s/plot.dat\' using 1:%d title \"exacte\" w lp pt 0 lc rgb \"black\"",
                    ppar->complete_path_output,
                    2*(ppar->option_godunov + ppar->option_rusanov + ppar->option_muscl) + 3);
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
