#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "g_quadrature.h"

int main() {
  // Teste
  double I = CVGaussianQuadrature();
  printf("I: %lf\n", I);


  return 0;
}
