#include <stdio.h>
#include "double_quad.h"

int main () {
  int k, N = 100;
  const double zUpper = 10.0, zLower = 0.0, dz = .05;
  double I, z = zUpper;

  do {
    I = DBGaussianQuad(N, z);
    printf("%lf\t%lf\n", z, I);
    z -= dz;
  } while (z >= zLower);

  return 0;
}
