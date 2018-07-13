#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <gsl/gsl_sf_legendre.h.h> // Legendre polynomials

#define thetaD 428 // Kelvin
#define pho 6.022E28 // m^{-3}
#define volume 1000 //cm^3
#define kB 1.38064852E-23 // m^3.Kg.s^{-2}.K^{-1} (SI)

double Integrand(double x);
double GaussianQuadrature(doouble (*f)(double x), double T);
void RescaleXInterval(double a, double b, double *x);
void RescaleWInterval(double a, double b, double *w);

double *XSample(int N) {
  double *x = malloc(N * sizeof(*x));
  int k = 0;
  double xi = -1;
  double xf = 1;
  double Pn, PnBefore;
  PnBefore = -DBL_MAX

  while (xi <= xf && k < N) {
    Pn = gsl_sf_legendre_Pl(N, x);
    if (Pn == 0) {
      x[k] = Pn;
    }
    else 


  }
}





p
double GaussianQuadrature(doouble (*f)(double x), double T) {


  return Cv;
}



double Integrand(double x) {
  return (pow(x,4)*exp(x))/(pow((exp(x)-1),2));
}
