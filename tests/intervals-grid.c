
/* include the required headers. */
#include "base.h"
#include "../src/intervals.h"

/* intervals-grid.x: test-case for interval set discretization operations.
 */
int main (int argc, char **argv) {
  unsigned int n_fails = 0;

  /* create an interval set for gridding. */
  intervals_t *I = intervals_new(10);
  intervals_union(I, 0.5, 1.0);
  intervals_union(I, 2.0, 3.0);
  intervals_union(I, 4.5, 5.0);
  intervals_union(I, 6.0, 6.5);

  /* compute the set of grid samples. */
  double samples[8];
  intervals_grid(I, samples, 5);

  /* test the result. */
  const double ans[] = { 0.5, 2.125, 2.75, 4.875, 6.5 };
  n_fails += test_eq_array_double(5, samples, ans, 1.0e-8);

  /* free the interval set. */
  intervals_free(I);

  return (n_fails > 0);
}

