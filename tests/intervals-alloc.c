
/* include the required headers. */
#include "base.h"
#include "../src/intervals.h"

/* intervals-alloc.x: test-case for interval set allocation and freeing.
 */
int main (int argc, char **argv) {
  unsigned int n_fails = 0;

  /* test allocation. */
  intervals_t *I = intervals_new(5);
  if (I) {
    n_fails += test_eq_uint(I->size, 0);
    n_fails += test_eq_uint(I->capacity, 5);
  }
  else
    n_fails++;

  /* test freeing. */
  intervals_free(I);
  if (I)
    n_fails++;

  return (n_fails > 0);
}

