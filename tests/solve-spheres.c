
/* include the required headers. */
#include "base.h"
#include "../src/enum.h"
#include "../src/enum-reduce.h"

/* the required function is not declared in the ibp-ng headers. */
int solve_tsi (vector_t *a, vector_t *b, vector_t *c,
               double ra,   double rb,   double rc,
               vector_t *n, vector_t *p, double *alpha);

/* solve-spheres.x: test-case for finding the points lying in the
 * intersection of three spheres.
 */
int main (int argc, char **argv) {
  unsigned int n_fails = 0;

  /* set up a three-vertex system. */
  vector_t a = { .x = 0.0, .y = 1.0, .z = 1.0 };
  vector_t b = { .x = 1.0, .y = 1.0, .z = 1.0 };
  vector_t c = { .x = 1.0, .y = 1.0, .z = 0.0 };

  /* define the sphere radii. */
  const double ra = 1.000000000000000;
  const double rb = 1.414213562373095;
  const double rc = 1.732050807568877;

  /* solve the three-sphere intersection problem. */
  double alpha;
  vector_t p, n;
  int ret = solve_tsi(&a, &b, &c, ra, rb, rc, &n, &p, &alpha);

  /* test the normal vector elements. */
  n_fails += test_eq_double(n.x, 0.0, 1.0e-8);
  n_fails += test_eq_double(n.y, 1.0, 1.0e-8);
  n_fails += test_eq_double(n.z, 0.0, 1.0e-8);

  /* test the plane projection elements. */
  n_fails += test_eq_double(p.x, 0.0, 1.0e-8);
  n_fails += test_eq_double(p.y, 1.0, 1.0e-8);
  n_fails += test_eq_double(p.z, 1.0, 1.0e-8);

  /* test the plane-to-solution length. */
  n_fails += test_eq_double(alpha, 1.0, 1.8e-8);

  /* test the return value (valid solution). */
  n_fails += test_eq_int(ret, 1);

  return (n_fails > 0);
}

