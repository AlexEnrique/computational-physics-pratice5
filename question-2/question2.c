#include <stdio.h>
#include "functions.h"
#include "trap_adaptative.h"

#define KB 1.38064852e-23 // J/K (SI)
#define C 299792458 // m/s
#define HBAR 1.05457180013e-34 //
#define PI M_PI

int main () {
  double I;

  I = EvalIntegralRomberg();
  printf("Integral via Romberg: %lf\n", I);

  I = EvalIntegralTrapezoidal();
  printf("Integral via Trapezoidal Method: %1.10lf\n", I);

  I = EvalIntegralQAGS();
  printf("Integral via adaptive integration with singularities: %1.10lf\n", I);

  double sigma = I * pow(KB, 4)/(4*pow(PI,2)*pow(C,2)*pow(HBAR,3));
  printf("\nStefan-Boltzmann constant (using the second value above): %1.5e\n", sigma);

  return 0;
}
