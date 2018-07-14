#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <gsl/gsl_sf_legendre.h.h> // Legendre polynomials

#define thetaD 428 // Kelvin
#define pho 6.022E28 // m^{-3}
#define V 1000 //cm^3
#define kB 1.38064852E-23 // m^3.Kg.s^{-2}.K^{-1} (SI)

double Integrand(double x);
double CVGaussianQuadrature(doouble (*f)(double x), double T);
void RescaleXInterval(double a, double b, double *x);
void RescaleWInterval(double a, double b, double *w);

// I have to count the number of steps of this func and see if it acceptable
double *XSample(int N) {
  double xi, xf, dx, Pn, *x;
  unsigned int k = 0; // index of x (xk)
  short *sign = malloc(2 * sizeof(*sign)); // compare the last 2 element's sign
  double PnBefore;

  // Variables initializations
  x = malloc(N * sizeof(*x)); // roots of Pn
  xi = -1;
  xf = 1;
  dx = (xf-xi)/N;
  PnBefore = -DBL_MAX // necessary? If not, exclude the <float.h> lib
  sign[0] = (short)(gsl_sf_legendre_Pl(N, xi) > 0);

  while (xi <= xf && k < N) {
    Pn = gsl_sf_legendre_Pl(N, xi);
    if (Pn == 0) {
      x[k++] = Pn;
    }
    else {
      sign[1] = (short)(Pn > 0);
      if (sign[1] != sign[0]) {
        // There exists a root between xi and xi-dx; find it
        // ...
        // Read GSL documentation for this part
        x[k++] = // the root
      }
      sign[0] = sign[1];
      // Problably the next line (also, PnBefore) will not be necessary
      // PnBefore = Pn;
    }
    xi += dx;
  }

  if (k < N-1) {
    printf("Alguma raíz não fora encontrada, k: %d,  N: %d\n", k, N);
  }

  // Plot the Pn function and see if the roots are alright as a check
  return x;
}



double CVGaussianQuadrature(doouble (*f)(double x), double T) {


  return Cv;
}



double Integrand(double x) {
  return (pow(x,4)*exp(x))/(pow((exp(x)-1),2));
}
