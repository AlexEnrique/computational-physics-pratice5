#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>

struct my_params {double a;};

double Integrand(double x, void *par) {
  if (x == 0 || fabs(x) <= GSL_DBL_EPSILON) {
    return (double)(0);
  }
  return (pow(x,3)/(exp(x)-1));
}

double RescIntegrd(double t, void *par) {
  struct my_params * params = (struct my_params *)par;
  double a = (params->a);
  return ( Integrand(t/(1-t), params)/pow(1-t,2) );
}

// This function doesn't work
double EvalIntegralRomberg() {
  const unsigned int n = 20; // number of max. iterations
  const double a = 0.0, b = .9999999999999999; // interval [a,b] of integration
  const double epsabs = 1E-5, epsrel = 1E-4; // max. abs. and rel. errors
  int status; // status of the integration
  double result; // the result of the integral

  size_t neval; // the number of evaluations of the integral
  gsl_integration_romberg_workspace *w;
  gsl_function F; // the integrand

  struct my_params params = {a};
  w = gsl_integration_romberg_alloc(n);
  F.function = &RescIntegrd;
  F.params = &params; // No parameters needed for the function used

  status = gsl_integration_romberg(&F, a, b, epsabs, epsrel, &result, &neval, w);

  if (status == GSL_EMAXITER) {
    printf("!! gsl_integration_romberg status: GSL_EMAXITER !!\n");
    printf("The number 'n' of max. iterations isn't sufficient.\nThe best current estimation was returned.\n");
    printf("If you want more acurated precision, increase the value of 'n'\n\n");
  }

  gsl_integration_romberg_free(w);
  return result;
}

double EvalIntegralQAGS() {
  const unsigned int n = 25;
  const double a = 0.0, epsabs = 1E-6, epsrel = epsabs;
  size_t limit = n;
  double result, abserr;
  int status;
  gsl_function F;
  gsl_integration_workspace *w;

  F.function = &Integrand;
  F.params = 0;
  w = gsl_integration_workspace_alloc(n);

  status = gsl_integration_qagiu(&F, a, epsabs, epsrel, limit, w, &result, &abserr);

  gsl_integration_workspace_free(w);
  return result;
}
