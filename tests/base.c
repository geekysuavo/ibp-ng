
/* include the base tests header. */
#include "base.h"

/* test_eq_array_double(): test if two arrays of doubles are
 * element-wise equal.
 *
 * arguments:
 *  @n: number of elements in both arrays.
 *  @u, @v: arrays of doubles to check.
 *  @tol: comparison tolerance.
 *
 * returns:
 *  integer indicating whether (1) or not (0) any pair of values differs
 *  more than the specified tolerance.
 */
unsigned int test_eq_array_double (unsigned int n,
                                   const double *u,
                                   const double *v,
                                   double tol) {
  for (unsigned int i = 0; i < n; i++) {
    if (fabs(u[i] - v[i]) > tol)
      return 1;
  }

  return 0;
}

/* test_eq_array_uint(): test if two arrays of uints are
 * element-wise equal.
 *
 * arguments:
 *  @n: number of elements in both arrays.
 *  @u, @v: arrays of uints to check.
 *
 * returns:
 *  integer indicating whether (1) or not (0) any pair of values differs.
 */
unsigned int test_eq_array_uint (unsigned int n,
                                 const unsigned int *u,
                                 const unsigned int *v) {
  for (unsigned int i = 0; i < n; i++) {
    if (u[i] != v[i])
      return 1;
  }

  return 0;
}

/* test_eq_uint(): test if two uints are equal.
 *
 * arguments:
 *  @u, @v: values to check for equality.
 *
 * returns:
 *  integer indicating whether (1) or not (0) the values differ.
 */
unsigned int test_eq_uint (unsigned int a, unsigned int b) {
  if (a != b)
    return 1;

  return 0;
}

/* test_eq_double(): test if two doubles are equal.
 *
 * arguments:
 *  @u, @v: values to check for equality.
 *  @tol: comparison tolerance.
 *
 * returns:
 *  integer indicating whether (1) or not (0) the values differ
 *  more than the specified tolerance.
 */
unsigned int test_eq_double (double a, double b, double tol) {
  if (fabs(a - b) > tol)
    return 1;

  return 0;
}

/* test_ptr_is_null_uint(): test if a pointer to a uint is null.
 *
 * arguments:
 *  @u: pointer to a uint to test against null.
 *
 * returns:
 *  integer indicating whether (1) or not (0) the pointer is non-null.
 */
unsigned int test_ptr_is_null_uint (const unsigned int *u) {
  if (u != NULL)
    return 1;

  return 0;
}

/* test_ptr_is_null_double(): test if a pointer to a double is null.
 *
 * arguments:
 *  @u: pointer to a double to test against null.
 *
 * returns:
 *  integer indicating whether (1) or not (0) the pointer is non-null.
 */
unsigned int test_ptr_is_null_double (const double *u) {
  if (u != NULL)
    return 1;

  return 0;
}

