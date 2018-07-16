#include <math.h>
#include <gsl/gsl_integration.h>

#define L 10.0 // m
#define M 1.0 // Kg
#define G 6.674E-11 // SI units
#define SIGMA (double)1000/100

double Integrand (double x, double y, double z);

double DBGaussianQuad(int N, double z) {
  const double lowerLim = (double) -L/2, upperLim = (double) L/2; // Same limits for x and y
  double xi, yj, I = 0;
  double deltaXY = (upperLim - lowerLim)/N;
  const gsl_integration_fixed_type *Itype;
  gsl_integration_fixed_workspace *w;

  Itype = gsl_integration_fixed_legendre;
  w = gsl_integration_fixed_alloc(Itype, N, lowerLim ,upperLim, 0, 0);

  double *weights = gsl_integration_fixed_weights(w);

  for (unsigned int i = 0; i < N; i++) {
    xi = lowerLim + i*deltaXY;
    for (unsigned int j = 0; j < N; j++) {
      yj = lowerLim + j*deltaXY;
      I += weights[i] * weights[j] * Integrand(xi, yj, z);
    }
  }


  return G*SIGMA*I*z;
}

double Integrand (double x, double y, double z) {
  return pow( (pow(x,2) + pow(y,2) + pow(z,2)), (double)-3/2 );
}
