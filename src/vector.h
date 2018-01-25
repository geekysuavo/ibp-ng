
/* ensure once-only inclusion. */
#pragma once

/* include the c math library header. */
#include <math.h>

/* include the traceback header. */
#include "trace.h"

/* vector_print(): macro to pretty-print vectors.
 */
#define vector_print(v)  vector_printfn(v, #v)

/* vector_free(): macro to free allocated vector structure pointers.
 */
#define vector_free(v)  { free(v); v = NULL; }

/* vector_t: structure for holding vertices in real(3).
 */
typedef struct {
  /* @x, @y, @z: elements of the vector.
   */
  double x, y, z;
}
vector_t;

/* function declarations (vector.c): */

vector_t *vector_new (void);

vector_t *vector_new_with_value (double x, double y, double z);

void vector_set (vector_t *v, double x, double y, double z);

void vector_dot (vector_t *a, vector_t *b, double *z);

void vector_axpy (vector_t *a, double alpha, vector_t *b);

double vector_sqdist (vector_t *a, vector_t *b);

double vector_dist (vector_t *a, vector_t *b);

void vector_normalize (vector_t *a);

void vector_cross (vector_t *a, vector_t *b, vector_t *z);

double vector_angle (vector_t *a, vector_t *b, vector_t *c);

double vector_dihedral (vector_t *a, vector_t *b,
                        vector_t *c, vector_t *d);

void vector_printfn (vector_t *v, const char *id);

