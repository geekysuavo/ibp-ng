
/* include the required headers. */
#include "base.h"
#include "../src/intervals.h"

/* intervals-grid.x: test-case for interval set discretization operations.
 */
int main (int argc, char **argv) {
  unsigned int n_fails = 0;

  /* create three interval sets for gridding. */
  intervals_t *I = intervals_new(10);
  intervals_t *J = intervals_new(10);
  intervals_t *K = intervals_new(10);

  /* initialize the first set of intervals. */
  intervals_union(I, 0.5, 1.0);
  intervals_union(I, 2.0, 3.0);
  intervals_union(I, 4.5, 5.0);
  intervals_union(I, 6.0, 6.5);

  /* initialize the second set of intervals. */
  intervals_union(J, 0.8, 0.8);
  intervals_union(J, 2.1, 2.1);
  intervals_union(J, 6.5, 6.5);

  /* compute the set of grid samples. */
  double samples[8];
  unsigned int n_samples = 5;
  intervals_grid(I, samples, &n_samples);

  /* test the result. */
  const double ans1[] = { 0.5, 2.125, 2.75, 4.875, 6.5 };
  n_fails += test_eq_uint(n_samples, 5);
  n_fails += test_eq_array_double(5, samples, ans1, 1.0e-8);

  /* intersect the interval set with the triple of points, then grid. */
  intervals_intersect(I, J, K);
  intervals_grid(K, samples, &n_samples);

  /* test the result. */
  const double ans2[] = { 0.8, 2.1, 6.5 };
  n_fails += test_eq_uint(n_samples, 3);
  n_fails += test_eq_array_double(3, samples, ans2, 1.0e-8);

  /* free the interval sets. */
  intervals_free(I);
  intervals_free(J);
  intervals_free(K);

  return (n_fails > 0);
}

