#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_sf_legendre.h> // Legendre polynomials
#include <gsl/gsl_roots.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>

#define thetaD 428 // Kelvin
#define pho 6.022E28 // m^{-3}
#define V 1000 //cm^3
#define kB 1.38064852E-23 // m^3.Kg.s^{-2}.K^{-1} (SI)

double Integrand(double x);
double CVGaussianQuadrature(double (*f)(double x), double T);
void RescaleXInterval(double a, double b, double *x);
void RescaleWInterval(double a, double b, double *w);
double PlDerivative(double x, void *par);
double Pl(double x, void *par);
void My_fdf(double x, void *par, double *f, double *df);

struct my_params {
  int N;
  double dx;
  };

double Pl(double x, void *par) {
  struct my_params * params = (struct my_params *)par;
  int N = (params->N);
  return gsl_sf_legendre_Pl(N, x);
}

double PlDerivative(double x, void *par) {
  struct my_params * params = (struct my_params *)par;
  int N = (params->N);
  double dx = (params->dx);

  return ( (gsl_sf_legendre_Pl(N, x+dx) - gsl_sf_legendre_Pl(N, x-dx))/(2*dx) );
}

void My_fdf(double x, void *par, double *f, double *df) {
  struct my_params * params = (struct my_params *)par;
  int N = (params->N);
  double t = gsl_sf_legendre_Pl(N, x);

  *f = t;
  *df = PlDerivative(x, params);
}

// I have to count the number of steps of this func and see if it acceptable
double *XSample(int N) {
  const gsl_root_fdfsolver_type *T;
  gsl_root_fdfsolver *s;
  int status, iter = 0, maxiter = 100;
  double xi, xf, dx, _x, Pn, *x;
  unsigned int k = 0; // index of x (xk)
  short *sign = malloc(2 * sizeof(*sign)); // compare the last 2 element's sign
  double PnBefore;
  gsl_function_fdf FDF;

  // Variables initializations
  x = malloc(N * sizeof(*x)); // roots of Pn
  xi = -1;
  xf = 1;
  dx = (xf-xi)/N;
  sign[0] = (short)(gsl_sf_legendre_Pl(N, xi) > 0);

  // GSL solver for root
  T = gsl_root_fdfsolver_newton;
  s = gsl_root_fdfsolver_alloc(T);
  // Legendre Pol derivative for gsl_root_fdfsolver
  struct my_params params = {N, dx};
  FDF.f = &Pl;
  FDF.df = &PlDerivative;
  FDF.fdf = &My_fdf;
  FDF.params = &params;

  while (xi <= xf && k < N) {
    Pn = gsl_sf_legendre_Pl(N, xi);

    if (Pn == 0) {
      x[k++] = Pn;
    }

    else {
      sign[1] = (short)(Pn > 0);

      if (sign[1] != sign[0]) {
        // There exists a root between xi and xi-dx; find it
        x[k] = xi;
        gsl_root_fdfsolver_set(s, &FDF, x[k]);
        do {
          iter++;
          status = gsl_root_fdfsolver_iterate(s);
          _x = x[k];
          x[k] = gsl_root_fdfsolver_root(s);
          status = gsl_root_test_delta(x[k], _x, 0, 1e-10);

        } while(status == GSL_CONTINUE && iter < maxiter);
        k++;
      }

      sign[0] = sign[1];
    }

    xi += dx;
  }

  if (k < N-1) {
    printf("Alguma raíz não fora encontrada, k: %d,  N: %d\n", k, N);
  }

  // Plot the Pn function and see if the roots are alright as a check
  gsl_root_fdfsolver_free(s);
  return x;
}


// Arrumar os CVGaussianQuadrature
#define CVGaussianQuadrature _CVGaussianQuadrature(Integrand, 273)
double _CVGaussianQuadrature(double (*f)(double x), double T) {
  int N = 50;
  double *x = malloc(N * sizeof(*x));
  x = XSample(N);

  return 0.0;
}


double Integrand(double x) {
  return (pow(x,4)*exp(x))/(pow((exp(x)-1),2));
}
