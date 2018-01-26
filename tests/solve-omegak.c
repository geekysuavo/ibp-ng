
/* include the required headers. */
#include "base.h"
#include "../src/enum.h"
#include "../src/enum-reduce.h"

/* the required function is not declared in the ibp-ng headers. */
int solve_iomega_k (vector_t *x1, vector_t *x2, vector_t *x3, vector_t *xk,
                    double d01, double d02, double lk, double uk,
                    intervals_t *iomega_k);

/* solve-omegak.x: test-case for computing the dihedral angles that
 * begin and end the solution arcs for a single interval reduction.
 */
int main (int argc, char **argv) {
  unsigned int n_fails = 0;

  /* set up the system of points. */
  vector_t x1 = { .x = 0.0, .y = 0.0, .z = 0.0 };
  vector_t x2 = { .x = 1.0, .y = 0.0, .z = 0.0 };
  vector_t x3 = { .x = 1.0, .y = 1.0, .z = 0.0 };
  vector_t xk = { .x = 1.0, .y = 2.0, .z = 3.0 };

  /* define the required distances. */
  const double d01 = 1.000000000000000;
  const double d02 = 1.414213562373095;
  double lk  = 2.906567063783021;
  double uk  = 3.114419338591080;

  /* solve for the two dihedral intervals. */
  intervals_t *I = intervals_new(4);
  int ret = solve_iomega_k(&x1, &x2, &x3, &xk, d01, d02, lk, uk, I);

  /* test the results. */
  const double ans1[] = { -1.727875959474387, -0.551870752379251 };
  const double ans2[] = { -1.413716694115407, -0.237711487020272 };
  n_fails += test_eq_uint(I->size, 2);
  n_fails += test_eq_array_double(2, I->start, ans1, 1.0e-8);
  n_fails += test_eq_array_double(2, I->end,   ans2, 1.0e-8);

  /* test the return value (valid solution). */
  n_fails += test_eq_int(ret, 1);

  /* use different bounds and solve again. */
  lk = 4.117052841406208;
  uk = 4.545143345531613;
  ret = solve_iomega_k(&x1, &x2, &x3, &xk, d01, d02, lk, uk, I);

  /* test the results. */
  const double ans3[] = { -M_PI, 0.861845941736156, 2.827433388230813 };
  const double ans4[] = { -2.827433388230814, 1.490164472454115, M_PI };
  n_fails += test_eq_uint(I->size, 3);
  n_fails += test_eq_array_double(3, I->start, ans3, 1.0e-8);
  n_fails += test_eq_array_double(3, I->end,   ans4, 1.0e-8);

  /* test the return value (valid solution). */
  n_fails += test_eq_int(ret, 1);

  /* free the interval set. */
  intervals_free(I);

  return (n_fails > 0);
}

