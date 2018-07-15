#include <stdio.h>
#include "functions.h"

int main () {
  double I = EvalIntegral();
  printf("Integral: %lf\n", I);
  
  return 0;
}
