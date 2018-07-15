#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

typedef struct {
  double a, b, error;
  int N;
  double *oldIPtr;
} f_args;

double var_f(double (*f)(double x), f_args in);
double f_base(double (*f)(double x), double a, double b, double error, int N, double *oldIPtr);
double EvalTrapezoidal(double (*f)(double x), double a, double b);

double var_f(double (*f)(double x), f_args in) {
  int N = in.N ? in.N : 1;
  int error = in.error ? in.error : (in.b-in.a)/N; // Error ~ h
  return f_base(f, in.a, in.b, error, N, in.oldIPtr);
}
#define TrapezoidalAdaptative(f, ...) var_f(f, (f_args){__VA_ARGS__})

double f_base(double (*f)(double x), double a, double b, double error, int N, double *oldIPtr) {
  double t, I;
  double h = (b-a)/N;

  I = 0;
  for (unsigned int i = 0; i < N; i++) {
    t = a + i*h;
    I += EvalTrapezoidal(*f, t, t+h);
  }

  if (fabs(I - *oldIPtr) > error || N == 1) {
    printf("I: %1.10lf,  error of I: %1.2e\n", I, fabs(I-*oldIPtr)/3);
    *oldIPtr = I;
    return TrapezoidalAdaptative(f, a, b, error, 2*N, oldIPtr);
  }
  return I;
}

double EvalTrapezoidal(double (*f)(double x), double a, double b) {
  return (b-a)/2*((*f)(a) + (*f)(b));
}
