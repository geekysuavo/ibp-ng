
/* include the required headers. */
#include "base.h"
#include "../src/enum.h"
#include "../src/enum-reduce.h"

/* the required function is not declared in the ibp-ng headers. */
void solve_3x3 (const double *A, const double *b, vector_t *p);

/* solve-linear.x: test-case for solving three-dimensional linear
 * systems using cramer's rule.
 */
int main (int argc, char **argv) {
  unsigned int n_fails = 0;

  /* define the solution vector and the left-hand matrix. */
  double x[] = { 1.5, 2.0, 3.1 };
  const double A[] = {
    0.8147, 0.9134, 0.2785,
    0.9058, 0.6324, 0.5469,
    0.1270, 0.0975, 0.9575
  };

  /* compute the right-hand vector. */
  double b[3] = { 0.0, 0.0, 0.0 };
  for (unsigned int i = 0; i < 3; i++)
    for (unsigned int j = 0; j < 3; j++)
      b[i] += A[3 * i + j] * x[j];

  /* solve the system and test the results. */
  vector_t p;
  solve_3x3(A, b, &p);
  n_fails += test_eq_double(p.x, 1.5, 1.0e-8);
  n_fails += test_eq_double(p.y, 2.0, 1.0e-8);
  n_fails += test_eq_double(p.z, 3.1, 1.0e-8);

  return (n_fails > 0);
}

