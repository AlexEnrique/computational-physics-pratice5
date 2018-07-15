#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>

#define KB 1.38064852E-23 // J/K (SI)
#define C 299 792 458 // m/s
#define HBAR 1.05457180013E-34 //
#define PI M_PI

double Integrand(double x) {
  if (x == 0) {
    return (double)(0);
  }
  return (pow(x,3)/(exp(x)-1));
}

double RescIntegrd(double t, void *par) {
  return (Integrand(t/(1-t))/(pow(1-t,2)));
}

// Int_0^1 RescIntegrd(t) dt
double EvalIntegral() {
  const unsigned int n = 30; // number of max. iterations
  const double a = 0, b = 1; // interval [a,b] of integration
  const double epsabs = 1E-9, epsrel = 1E-10; // max. abs. and rel. errors
  int status; // status of the integration
  double result; // the result of the integral

  size_t neval; // the number of evaluations of the integral
  gsl_integration_romberg_workspace *w;
  gsl_function F; // the integrand

  w = gsl_integration_romberg_alloc(n);
  F.function = &RescIntegrd;
  F.params = 0; // No parameters needed for the function used

  do {
    status = gsl_integration_romberg(&F, a, b, epsabs, epsrel, &result, &neval, w);
    if (status == GSL_EMAXITER) {
      printf("\n!! gsl_integration_romberg status: GSL_EMAXITER !!\n");
      printf("The number 'n' of max. iterations isn't sufficient.\nThe best current estimation was returned.\n");
      printf("If you want more acurated precision, increase the value of 'n'\n\n");
    }

  } while (status == GSL_CONTINUE);

  gsl_integration_romberg_free(w);
  return result;
}
