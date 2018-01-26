
/* ensure once-only inclusion. */
#pragma once

/* include c library headers. */
#include <stdio.h>
#include <math.h>

/* function declarations (base.c): */

unsigned int test_eq_array_double (unsigned int n,
                                   const double *u,
                                   const double *v,
                                   double tol);

unsigned int test_eq_array_uint (unsigned int n,
                                 const unsigned int *u,
                                 const unsigned int *v);

unsigned int test_eq_int (int a, int b);

unsigned int test_eq_uint (unsigned int a, unsigned int b);

unsigned int test_eq_double (double a, double b, double tol);

unsigned int test_ptr_is_null_uint (const unsigned int *u);

unsigned int test_ptr_is_null_double (const double *u);

