#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>



//-----------------------------------------------------------------------------
// Saint Venant 1d
//-----------------------------------------------------------------------------

void flux_riem_2d(double *wL, double *wR, double *vnorm, double *flux){

  double qnL = wL[1] * vnorm[0] +  wL[2] * vnorm[1]; 
  double qnR = wR[1] * vnorm[0] +  wR[2] * vnorm[1];

  double qtL = -wL[1] * vnorm[1] +  wL[2] * vnorm[0];
  double qtR = -wR[1] * vnorm[1] +  wR[2] * vnorm[0];

  double vL[2] = {wL[0], qnL};
  double vR[2] = {wR[0], qnR};

  double v[2];
  double xi = 0;

  riem_stvenant(vL, vR, xi, v);

  double un = v[1] / v[0];

  double ut;

  if (un > 0)
    ut = qtL / wL[0];
  else
    ut = qtR / wR[0];

  double qn = v[1];
  double qt = ut * v[0];

  double w[3];
  w[0] = v[0];
  w[1] = qn * vnorm[0] - qt * vnorm[1];
  w[2] = qn * vnorm[1] + qt * vnorm[0];

}

void riem_stvenant(double *wL,
		   double *wR,
		   double xi,
		   double *w){

  double hL = wL[0];
  double uL = wL[1]/wL[0];
  double hR = wR[0];
  double uR = wR[1]/wR[0];

  double hs = 1e-6;
  int itermax = 10;
  
  for(int it = 0; it < itermax; it++){
    double f = uL - (hs - hL) * Z(hs, hL) -
      uR - (hs - hR) * Z(hs, hR);
    double df = -(hs - hL) * dZ(hs, hL) -
      Z(hs, hL) -
      (hs - hR) * dZ(hs, hR) -
      Z(hs, hR);
    double dhs = -f / df;

    hs += dhs;

    printf("it=%d f=%e df=%e hs=%e dhs=%e\n",it,f,df,hs,dhs);
    
  }

  double us = uL - (hs - hL) * Z(hs, hL);

  double v1m, v1p, v2m, v2p;

  // 1-onde
  if (hs < hL){
    v1m = uL - sqrt(g * hL);
    v1p = us - sqrt(g * hs);
  } else {
    double a = sqrt(hs) / (sqrt(hs) + sqrt(hL));
    double u = a * us + (1 - a) * uL;
    double h = (hs + hL) / 2;
    v1m = u - sqrt(g * h);
    v1p = v1m;
  }

  // 2 onde
  if (hs < hR){
    v2m = us + sqrt(g * hs);
    v2p = uR + sqrt(g * hR);
  } else {
    double a = sqrt(hs) / (sqrt(hs) + sqrt(hR));
    double u = a * us + (1 - a) * uR;
    double h = (hs + hR) / 2;
    v2m = u + sqrt(g * h);
    v2p = v2m;
  }

  //printf("v=%f %f %f %f\n hs=%f us=%f\n", v1m,v1p,v2m,v2p, hs,us);
  
  if (xi < v1m) {
    w[0] = wL[0];
    w[1] = wL[1];
  } else if (xi < v1p){
    double u = (uL + 2 * xi + 2 *sqrt(g * hL)) / 3;
    double h = (u - xi) * (u - xi) / g;
    w[0] = h;
    w[1] = h * u;
  } else if (xi < v2m){
    w[0] = hs;
    w[1] = hs * us;
  } else if (xi < v2p){
    double u = (uR + 2 * xi - 2 *sqrt(g * hR)) / 3;
    double h = (u - xi) * (u - xi) / g;
    w[0] = h;
    w[1] = h * u;
  } else {
    w[0] = wR[0];
    w[1] = wR[1];
  }

}

void plot_riem(double *wL,
	       double *wR){

  FILE *plotfile;

  double xmin = -5;
  double xmax = 5;

  int n = 1000;
  double dx = (xmax - xmin)/n;
  plotfile = fopen("riem.dat", "w");
  for(int i = 0; i < n; i++){
    double xi = xmin+i*dx;
    double w[2];
    riem_stvenant(wL,wR,xi,w);
    double h = w[0];
    double u = w[1]/w[0];
    fprintf(plotfile, "%f %f %f\n",
	    xi, h, u);
  }
  fclose(plotfile);

  system("gnuplot riemcom");

}

double Heaviside(double x){
  if (x > 0)
    return 1;
  else
    return 0;
}

double Dirac(double x){
    return 0;
}


double Z(double hs,
	 double h){

double t0 = 2.0*sqrt(g)/(sqrt(hs) + sqrt(h))*Heaviside(h-hs) + sqrt(2.0)*sqrt(g)*sqrt(h+hs)/sqrt(h*hs)/2.0
            - sqrt(2.0)*sqrt(g)*sqrt(h+hs)/sqrt(h*hs)*Heaviside(h-hs)/2.0;

 return t0;
  
}

double dZ(double hs,
	  double h){

   double t0 = - sqrt(g)/pow(sqrt(hs) + sqrt(h),2.0)*Heaviside(h-hs)/sqrt(hs) - 2.0*sqrt(g)/(sqrt(hs) + sqrt(h))*Dirac(-h+hs) + sqrt(2.0)*sqrt(g)/sqrt(h+hs)/sqrt(h*hs)/4.0
                - sqrt(2.0)*sqrt(g)*sqrt(h+hs)/sqrt(h*h*h*hs*hs*hs)*h/4.0 - sqrt(2.0)*sqrt(g)/sqrt(h+hs)/sqrt(h*hs)*Heaviside(h-hs)/4.0
                + sqrt(2.0)*sqrt(g)*sqrt(h+hs)/sqrt(h*h*h*hs*hs*hs)*Heaviside(h-hs)*h/4.0+sqrt(2.0)*sqrt(g)*sqrt(h+hs)/sqrt(h*hs)*Dirac(-h+hs)/2.0;

   return t0;
}

//-----------------------------------------------------------------------------
// Equation of Saint-Venant

double lambda_ma_sv_1d(double a){
     return a;
}

double Riemann_sv_1d(double u_L, double u_R, double z){

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

double fluxnum_gd_sv_1d(double a, double b){
    
    double r = Riemann_sv_1d(a, b, 0);
    
    return r * r / 2;
}

double fluxnum_ru_sv_1d(double a, double b){

    return 1/2 * (b-a) * ( (a+b)/2 - fmax(fabs(a),fabs(b)) );
}

// Example 1 of resolution of equation of burgers avec u_L > u_R

#define u_L1 2
#define u_R1 1

double solexacte_sv(double x, double t){

    double sigma = (u_L1+u_R1) / 2;

    if (x < sigma * t){
        return u_L1;
    } else {
        return u_R1;
    }
}

double boundary_spatial_sv(double x){

    return solexacte_burgers1(x, 0);
}

double boundary_temporal_left_sv(double xmin, double t){

    return solexacte_sv(xmin, t);
}

double boundary_temporal_right_sv(double xmax, double t){

    return solexacte_sv(xmax, t);
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
