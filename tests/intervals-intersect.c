
/* include the required headers. */
#include "base.h"
#include "../src/intervals.h"

/* intervals-intersect.x: test-case for interval set intersection operations.
 */
int main (int argc, char **argv) {
  unsigned int n_fails = 0;

  /* allocate three interval sets for the operation. */
  intervals_t *Ia, *Ib, *Ic;
  Ia = intervals_new(10);
  Ib = intervals_new(10);
  Ic = intervals_new(10);

  /* add intervals to the first set. */
  intervals_union(Ia, 0.0, 1.0);
  intervals_union(Ia, 2.0, 3.0);
  intervals_union(Ia, 4.0, 5.0);
  intervals_union(Ia, 6.0, 7.0);

  /* add intervals to the second set. */
  intervals_union(Ib, 0.5, 1.0);
  intervals_union(Ib, 1.8, 3.4);
  intervals_union(Ib, 4.5, 6.5);

  /* intersect the two sets. */
  intervals_intersect(Ia, Ib, Ic);

  /* test the results of the intersection. */
  const double ans1[] = { 0.5, 2.0, 4.5, 6.0 };
  const double ans2[] = { 1.0, 3.0, 5.0, 6.5 };
  n_fails += test_eq_uint(Ic->size, 4);
  n_fails += test_eq_array_double(4, Ic->start, ans1, 1.0e-8);
  n_fails += test_eq_array_double(4, Ic->end,   ans2, 1.0e-8);

  /* free the interval sets. */
  intervals_free(Ia);
  intervals_free(Ib);
  intervals_free(Ic);

  return (n_fails > 0);
}

