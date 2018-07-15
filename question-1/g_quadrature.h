#include <math.h>
#include <gsl/gsl_integration.h>

#define THETA_D 428 // Kelvin
#define PHO 6.022E28 // m^{-3}
#define V 1000/1E6 //cm^3
#define KB 1.38064852E-23 // J/K (SI)

double Integrand(double x, void *par) {
  return pow(x,4)*exp(x)/(pow((exp(x)-1),2));
}

double CvQuadrature(double T) {
  double Cv;
  double upperLim = (double)THETA_D/T;
  unsigned int N = 50;
  const gsl_integration_fixed_type *Itype;
  gsl_integration_fixed_workspace *w;
  gsl_function func;


  func.function = &Integrand;
  func.params = 0;
  Itype = gsl_integration_fixed_legendre;

  w = gsl_integration_fixed_alloc(Itype, N, 0 ,upperLim, 0, 0);

  gsl_integration_fixed(&func, &Cv, w);
  Cv *= 9*V*PHO*KB*pow(T/THETA_D,3);

  gsl_integration_fixed_free(w);
  return Cv;
}
