#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "g_quadrature.h"

int main() {
  // Teste 
  int N = 50;
  double *x = malloc(N * sizeof(*x));
  x = XSample(N);

  return 0;
}
