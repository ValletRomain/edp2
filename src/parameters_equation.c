#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>


//-----------------------------------------------------------------------------
// Equation of transport

double Riemann_global_transport(double u_L, double u_R, double z, double c){

    double r;

    if (u_L < u_R){
        if (z < u_L)
            r = u_L;
        else if (z > u_R)
            r = u_R;
        else
            r = z;
    }
    else {        
        if (z < c)
            r = u_L;
        else 
            r = u_R;
    }

    return r;
}

double fluxnum_global_transport_ru(double a, double b, double c){

    if (c > 0)
        return c * a;
    else
        return c * b;
    
}

// Example 1 of resolution of equation of transport

// vitesse de transport
#define _C1 1

double Riemann_transport1(double u_L, double u_R, double z){

    return Riemann_global_transport(u_L, u_R, z, _C1);
}

double lambda_ma_trans1(double u){
     return _C1;
}

double fluxnum_gd_trans1(double a, double b){
    
    return _C1 * Riemann_transport1(a, b, 0);
}

double fluxnum_ru_trans1(double a, double b){

    return fluxnum_global_transport_ru(a, b, _C1);
}

double solexacte_trans1(const double x, const double t){

     double uL = exp(x/_C1 - t);
     double uR = 0;

     if (x < _C1 * t){
        return uL;
     } else {
        return uR;
     }

}

double boundary_spatial_trans1(double x){

    return solexacte_trans1(x, 0);
}

double boundary_temporal_left_trans1(double xmin, double t){

    return solexacte_trans1(xmin, t);
}

double boundary_temporal_right_trans1(double xmax, double t){

    return solexacte_trans1(xmax, t);
}

//-----------------------------------------------------------------------------
// Equation of Burgers

double lambda_ma_burgers(double a){
     return a;
}

double Riemann_burgers(double u_L, double u_R, double z){

    double sigma = (u_L + u_R) / 2;
    double r;

    if (u_L < u_R){
        if (z < u_L)
            r = u_L;
        else if (z > u_R)
            r = u_R;
        else 
            r = z;
    }
    else {
        if (z < sigma)
            r = u_L;
        else
            r = u_R;
    }

    return r;
}

double fluxnum_gd_burgers(double a, double b){
    
    double r = Riemann_burgers(a, b, 0);
    
    return r * r / 2;
}

double fluxnum_ru_burgers(double a, double b){

    return 1/2 * (b-a) * ( (a+b)/2 - fmax(fabs(a),fabs(b)) );
}

// Example 1 of resolution of equation of burgers avec u_L > u_R

#define u_L1 2
#define u_R1 1

double solexacte_burgers1(double x, double t){

    double sigma = (u_L1+u_R1) / 2;

    if (x < sigma * t){
        return u_L1;
    } else {
        return u_R1;
    }
}

double boundary_spatial_burgers1(double x){

    return solexacte_burgers1(x, 0);
}

double boundary_temporal_left_burgers1(double xmin, double t){

    return solexacte_burgers1(xmin, t);
}

double boundary_temporal_right_burgers1(double xmax, double t){

    return solexacte_burgers1(xmax, t);
}

// Example 1 of resolution of equation of burgers avec u_L < u_R

#define u_L2 2
#define u_R2 1

double solexacte_burgers2(double x, double t){

    double sigma = (u_L1+u_R1) / 2;

     if ((x < u_L2 * t) || ((x==0) && (t==0))){ // je force sole(0,0) = u_L
         return u_L1;
     }
     else if (x > u_R2 * t){
         return u_R1;
     }
     else {
         return x / t;
     }
}

double boundary_spatial_burgers2(double x){

    return solexacte_burgers2(x, 0);
}

double boundary_temporal_left_burgers2(double xmin, double t){

    return solexacte_burgers2(xmin, t);
}

double boundary_temporal_right_burgers2(double xmax, double t){

    return solexacte_burgers2(xmax, t);
}


//-----------------------------------------------------------------------------
// Pour la méthode MUSCL
//-----------------------------------------------------------------------------

double minmod(double a, double b, double c){

    if ((a < 0) && (b < 0) && (c < 0))
        return fmin(fmin(a, b), c);
    else if ((a > 0) && (b > 0) && (c > 0))
        return fmax(fmax(a, b), c);
    else
        return 0;
}

double w_half_l(double w_il1, double w_i, double w_ip1){
    // Remplie dest tel que dest[0] = a_{1+1/2},-

    double r = minmod(w_i-w_il1, w_ip1-w_i, (w_ip1-w_il1)/(double)2);

    return w_i + r/(double)2;
}

double w_half_p(double w_i, double w_ip1, double w_ip2){
    // Remplie dest tel que dest[0] = a_{1+1/2},-

    double r = minmod(w_ip1-w_i, w_ip2-w_ip1, (w_ip2-w_i)/(double)2);

    return w_ip1 - r/(double)2;
}

//-----------------------------------------------------------------------------
// Calcul des normes
//-----------------------------------------------------------------------------

double norm_L1(int I, double * a, double * b){
    // Calcul l'erreur en norme L_1 entre le tableau a et le tableau b

     double sum = 0;
     
     for (int i=0; i<I; i++)
            sum += fabs(a[i] - b[i]);

     return sum / I;
}

double norm_L2(int I, double * a, double * b){
    // Calcul l'erreur en norme L_1 entre le tableau a et le tableau b

     double sum = 0;

     for (int i=0; i<I; i++)
            sum += (a[i] - b[i]) * (a[i] - b[i]);

     return sum / I;
}

double norm_inf(int I, double * a, double * b){
    // Calcul l'erreur en norme L_1 entre le tableau a et le tableau b

    double norm = fabs(a[0] - b[0]);
    double inter;

    for (int i=0; i<I; i++){
        inter = fabs(a[i] - b[i]);
        if (inter > norm)
            norm = inter;
    }
    return norm;
}
