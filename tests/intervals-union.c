
/* include the required headers. */
#include "base.h"
#include "../src/intervals.h"

/* intervals-union.x: test-case for the interval set union operation.
 */
int main (int argc, char **argv) {
  unsigned int n_fails = 0;

  /* allocate a new interval set. */
  intervals_t *I = intervals_new(5);

  /* add the first interval. */
  intervals_union(I, 1.0, 2.0);
  n_fails += test_eq_uint(I->size, 1);
  n_fails += test_eq_double(I->start[0], 1.0, 1.0e-8);
  n_fails += test_eq_double(I->end[0],   2.0, 1.0e-8);

  /* add three more intervals, all mutually disjoint. */
  intervals_union(I, 2.1, 3.0);
  intervals_union(I, 4.0, 5.0);
  intervals_union(I, 6.0, 7.0);

  /* test the results of the unions. */
  const double ans1[] = { 1.0, 2.1, 4.0, 6.0 };
  const double ans2[] = { 2.0, 3.0, 5.0, 7.0 };
  n_fails += test_eq_uint(I->size, 4);
  n_fails += test_eq_array_double(4, I->start, ans1, 1.0e-8);
  n_fails += test_eq_array_double(4, I->end,   ans2, 1.0e-8);

  /* add an overlapping interval and test the results. */
  intervals_union(I, 2.4, 4.5);
  const double ans3[] = { 1.0, 2.1, 6.0 };
  const double ans4[] = { 2.0, 5.0, 7.0 };
  n_fails += test_eq_uint(I->size, 3);
  n_fails += test_eq_array_double(3, I->start, ans3, 1.0e-8);
  n_fails += test_eq_array_double(3, I->end,   ans4, 1.0e-8);

  /* free the interval set. */
  intervals_free(I);

  return (n_fails > 0);
}

