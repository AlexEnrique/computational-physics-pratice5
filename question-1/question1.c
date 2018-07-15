#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "g_quadrature.h"
#include <gsl/gsl_integration.h>


double f(double x, void *par) {
  return pow(sin(sqrt(100*x)),2);
}

int main() {
  // Teste
  double a = 0.0;
  double b = 1.0;
  double I;
  const gsl_integration_fixed_type *T;
  gsl_integration_fixed_workspace *w;
  T = gsl_integration_fixed_legendre;
  w = gsl_integration_fixed_alloc(T, 20, a ,b, 0, 0);

  gsl_function func;
  func.function = &f;
  func.params = 0;

  gsl_integration_fixed(&func, &I, w);
  printf("I; %1.20lf\n", I);


  gsl_integration_fixed_free(w);

  return 0;
}
