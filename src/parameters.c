#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>

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

void boundary_spatial_trans1(double x, double *w){

    solexacte_trans1(x, 0, w);
}

void boundary_temporal_left_trans1(double xmin, double t, double *w){

    solexacte_trans1(xmin, t, w);
}

void boundary_temporal_right_trans1(double xmax, double t, double *w){

    solexacte_trans1(xmax, t, w);
}


//-----------------------------------------------------------------------------
// Calcul des normes
//-----------------------------------------------------------------------------

double norm_L1(int I, int m, double * a, double * b){
    // Calcul l'erreur en norme L_1 entre le tableau a et le tableau b

     double sum = 0;
     
     for (int i=0; i<I; i++){
         for (int iv=0; iv<m; iv++){
            sum += fabs(a[i*m+iv] - b[i*m+iv]);
         }
     }

     return sum / I / m;
}

double norm_L2(int I, int m, double * a, double * b){
    // Calcul l'erreur en norme L_1 entre le tableau a et le tableau b

     double sum = 0;

     for (int i=0; i<I; i++){
         for (int iv=0; i<m; i++){
            sum += (a[i*m+iv] - b[i*m+iv]) * (a[i*m+iv] - b[i*m+iv]);
         }
     }

     return sum / I / m;
}
