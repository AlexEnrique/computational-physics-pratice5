#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "g_quadrature.h"

int main() {
  double T = 5, TMax = 500, dT = 0.1; // Kelvin
  double Cv;

  do {
    Cv = CvQuadrature(T);
    printf("%lf\t%1.5e\n", T, Cv);

    T += dT; 
  } while (T <= TMax);


  return 0;
}
