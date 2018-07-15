#include <stdio.h>
#include <gsl/gsl_sf_legendre.h> // Legendre polynomials

// #define L 200 // Test
double *XSample (int N);
double Bisection(double (*f)(int N, double x), double a, double b, int N);
double PnFunc(int L, double x);
short Sign(double a);

double *XSample (int N) {
  double xi = -1;
  double xf = 1;
  double dx = (double)2/(10*N);
  double *r, y;
  int k = 0;
  short *sign = malloc(2 * sizeof(*sign));

  r = malloc(N * sizeof(*r));
  sign[0] = Sign(PnFunc(N, xi));

  while (xi <= xf) {
    y = PnFunc(N, xi);

    if (y == 0) {
      printf("%lf\t%lf\n", xi, 0.0);
    }

    else {
      sign[1] = Sign(PnFunc(N, xi));

      if (sign[0] != sign[1]) {
        r[k++] = Bisection(PnFunc, xi-dx, xi, N);
        // printf("%lf\t%lf\n", r[k-1], 0.0);
      }
      sign[0] = sign[1];
    }

    xi += dx;
  }

  return r;
}


double PnFunc(int L, double x) {
  return gsl_sf_legendre_Pl(L, x);
}


short Sign(double a) {
  return ( (short)(a / fabs(a)) );
}


double Bisection(double (*f)(int N, double x), double a, double b, int N) {
  int n = 0;
  double c, tol = 1e-8;

  do {
    c = a + (b-a)/2;
    if (Sign((*f)(N, a)) == Sign((*f)(N, c))) {
      a = c;
    }

    else {
      b = c;
    }
  } while ( pow(2, -(++n))*(b-a) < tol );


  return c;
}
