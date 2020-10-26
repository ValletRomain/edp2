#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>

//-----------------------------------------------------------------------------
// Parameters of equation
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Example 1 of resolution of equation of transport

// vitesse de transport
#define _C1 1

double lambda_ma_trans1(double *u){
     return _C1;
}

void fluxnum_trans1(double *a, double *b, double vmax, double *flux){
    
     flux[0] = _C1 * a[0];
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
// Example 1 of resolution of equation of burgers

#define u_L1 2
#define u_R1 1

double lambda_ma_burgers1(double *a){
     return a[0];
}

void fluxnum_burgers1(double *a, double *b, double vmax, double *flux){
    
     flux[0] = (b[0]-a[0]) / 2 * ((a[0]+b[0])/2 - vmax);
}

void solexacte_burgers1(double x, double t, double *w){

    double sigma = (u_L1+u_R1) / 2;

     if (x < sigma * t){
         w[0] = u_L1;
     } else {
         w[0] = u_R1;
     }

}

void boundary_spatial_burgers1(double x, double *w){

    solexacte_trans1(x, 0, w);
}

void boundary_temporal_left_burgers1(double xmin, double t, double *w){

    solexacte_burgers1(xmin, t, w);
}

void boundary_temporal_right_burgers1(double xmax, double t, double *w){

    solexacte_burgers1(xmax, t, w);
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
