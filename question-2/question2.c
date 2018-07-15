#include <stdio.h>
#include "functions.h"
#include "trap_adaptative.h"

#define KB 1.38064852E-23 // J/K (SI)
#define C 299792458 // m/s
#define HBAR 1.05457180013E-34 //
#define PI M_PI

int main () {
  // double I = EvalIntegral1();
  double I = EvalIntegral2();
  printf("Integral: %lf\n", I);
  double sigma = I * pow(KB, 4)/(4*pow(PI,2)*pow(C,2)*pow(HBAR,3));
  printf("Stefan-Boltzmann constant: %1.5e\n", sigma);

  return 0;
}
