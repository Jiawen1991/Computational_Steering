#include "polifitgsl.h"
#include <gsl/gsl_complex.h>
#include <gsl/gsl_poly.h>
#include <stdio.h>
 
bool polynomialfit(int obs, int degree, 
		   double *dx, double *dy, double *store) /* n, p */
{
  gsl_multifit_linear_workspace *ws;
  gsl_matrix *cov, *X;
  gsl_vector *y, *c;
  double chisq;
 
  int i, j;
 
  X = gsl_matrix_alloc(obs, degree);
  y = gsl_vector_alloc(obs);
  c = gsl_vector_alloc(degree);
  cov = gsl_matrix_alloc(degree, degree);
 
  for(i=0; i < obs; i++) {
    for(j=0; j < degree; j++) {
      gsl_matrix_set(X, i, j, pow(dx[i], j));
    }
    gsl_vector_set(y, i, dy[i]);
  }
 
  ws = gsl_multifit_linear_alloc(obs, degree);
  gsl_multifit_linear(X, y, c, cov, &chisq, ws);
 
  /* store result ... */
  for(i=0; i < degree; i++)
  {
    store[i] = gsl_vector_get(c, i);
  }
 
  gsl_multifit_linear_free(ws);
  gsl_matrix_free(X);
  gsl_matrix_free(cov);
  gsl_vector_free(y);
  gsl_vector_free(c);
  return true; /* we do not "analyse" the result (cov matrix mainly)
		  to know if the fit is "good" */
}

 
 
double x[] = {0,  1,  2,  3,  4,  5,  6,   7,   8,   9,   10};
double y[] = {1,  6,  17, 34, 57, 86, 121, 162, 209, 262, 321};
int NP = sizeof(x)/sizeof(*x); 
int DEGREE = 3;

double coef[128], sol[128];
    gsl_complex csol[128];
 
int main()
{
  int i;
double coeff[DEGREE];
  polynomialfit(NP, DEGREE, x, y, coeff);
  printf("Coefficients:\n");
  for(i=0; i < DEGREE; i++) {
    printf("%lf\n", coeff[i]);
  }


/* Solve Linear and Quadratic Equations */
    printf("Solve Quadratic Equations\n");
 
    /* (1) 3 * x - 2 = 0 */
    gsl_poly_solve_quadratic(0.0, 3.0, -2.0, &sol[0], NULL);
    printf("x = %g\n", sol[0]);
 
    /* (2) 2 * x^2 - 3 * x - 5 = 0 */
    gsl_poly_solve_quadratic(2.0, -3.0, -5.0, &sol[0], &sol[1]);
    printf("x1, x2 = %g, %g\n", sol[0], sol[1]);
 
    /* (3) 3 * x^2 - 5 * x + 9 = 0 */
    gsl_poly_complex_solve_quadratic(3.0, -5.0, 9.0, &csol[0], &csol[1]);
    printf("re(x1), im(x1) = %g, %g\n", GSL_REAL(csol[0]), GSL_IMAG(csol[0]));
    printf("re(x2), im(x2) = %g, %g\n", GSL_REAL(csol[1]), GSL_IMAG(csol[1]));
 
/* Solve Cubic Equation */
    printf("\nSolve Monic Cubic Equation\n");
 
    /* x^3 + 2 * x^2 + 3 * x + 4 = 0 */
    gsl_poly_complex_solve_cubic(2.0, 3.0, 4.0, &csol[0], &csol[1], &csol[2]);
    for(i = 0; i < 3; i++)
        printf("re(x%d), im(x%d) = %g, %g\n", i, i, GSL_REAL(csol[i]), GSL_IMAG(csol[i]));


}
