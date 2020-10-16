#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>

#include "godunov.h"

// vitesse de transport
#define _C1 1

//-----------------------------------------------------------------------------
// Parameters of equation
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Example 1 of resolution of equation of transport

double speed_trans1(void){
     return _C1;
}

void fluxnum_trans1(double *a, double *b, double *flux){
    
     flux[0] = _C1 * a[0];
}

double lambda_ma_trans1(double *u){
     return _C1;
}

void solexacte_trans1(double x, double t, double *w){

     double uL = exp(x/_C1 - t);
     double uR = 0;

     if (x < _C1 * t){
         w[0] = uL;
     } else {
         w[0] = uR;
     }

}


//-----------------------------------------------------------------------------
// Calcul des normes
//-----------------------------------------------------------------------------

double error_L1(int I, int m, double * a, double * b){
    // Calcul l'erreur en norme L_1 entre le tableau a et le tableau b

     double norme_L1 = 0;
     
     for (int i=0; i<I; i++){
         for (int iv=0; iv<m; iv++){
            norme_L1 += abs(a[i*m+iv] - b[i*m+iv]);
         }
     }

     return sqrt(norme_L1);
}

double error_L2(int I, int m, double * a, double * b){
    // Calcul l'erreur en norme L_1 entre le tableau a et le tableau b

     double norme_L1 = 0;

     for (int i=0; i<I; i++){
         for (int iv=0; i<m; i++){
            norme_L1 += (a[i*m+iv] - b[i*m+iv]) * (a[i*m+iv] - b[i*m+iv]);
         }
     }

     return norme_L1;
}


//-----------------------------------------------------------------------------
// Manage of parameters
//-----------------------------------------------------------------------------

void godunov_parameters(godunov * pgd, char * option){

    if (option = "transport_1d_1"){
        pgd->pspeed = speed_trans1;
        pgd->pfluxnum = fluxnum_trans1;
        pgd->plambda_ma = lambda_ma_trans1;
        pgd->psolexacte = solexacte_trans1;
    }
}

void godunov_error_parameters(godunov_error * pgderr, char * option_error){

    if (option_error = "norm_L1"){
        pgderr->perror = error_L1;
    }
    else if (option_error = "norm_L2"){
        pgderr->perror = error_L2;
    }
}