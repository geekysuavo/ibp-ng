
/* ensure once-only inclusion. */
#pragma once

/* include standard c library headers. */
#include <stdio.h>
#include <math.h>

/* define pi, if necessary. */
#ifndef M_PI
#define M_PI  3.14159265358979323846
#endif
#ifndef M_PI_2
#define M_PI_2  1.57079632679489661923
#endif

/* value_print(): macro to pretty-print values.
 */
#define value_print(v)  value_printfn(v, #v)

/* value_type_t: enumeration of all types of values that may
 * be stored within a value_t structure.
 */
typedef enum {
  VALUE_TYPE_UNDEFINED = 0,
  VALUE_TYPE_SCALAR,
  VALUE_TYPE_INTERVAL
}
value_type_t;

/* value_semantic_t: enumeration of all semantic meanings of values
 * that may be stored within a value_t.
 */
typedef enum {
  VALUE_IS_DISTANCE = 0,
  VALUE_IS_ANGLE,
  VALUE_IS_DIHEDRAL
}
value_semantic_t;

/* value_t: generalized parameter type structure that holds either
 * scalar values or interval values.
 */
typedef struct value value_t;
struct value {
  /* @type: type of value contained by the structure.
   */
  value_type_t type;

  /* @l: scalar exact value, interval lower bound.
   * @u: scalar exact value, interval upper bound.
   */
  double l, u;

  /* @sem: semantic meaning of the source value.
   * @src: source value from which this value was derived.
   */
  value_semantic_t sem;
  value_t *src;
};

/* function declarations: */

value_t value_undefined (void);

value_t value_scalar (double v);

value_t value_interval (double l, double u);

int value_is_undefined (value_t val);

int value_is_scalar (value_t val);

int value_is_interval (value_t val);

int value_is_distance (value_t val);

int value_is_angle (value_t val);

int value_is_dihedral (value_t val);

void value_set_source (value_t *v, value_t *vsrc,
                       const value_semantic_t sem);

value_t value_add (value_t va, value_t vb);

value_t value_sub (value_t va, value_t vb);

value_t value_mul (value_t va, value_t vb);

value_t value_div (value_t va, value_t vb);

value_t value_pow (value_t v, double p);

value_t value_scal (value_t v, double p);

value_t value_bound (value_t v, value_t b);

value_t value_intersect (value_t va, value_t vb);

value_t value_sin (value_t v);

value_t value_cos (value_t v);

void value_printfn (value_t v, const char *id);

value_t value_from_angle (value_t a, value_t b,
                          value_t theta);

value_t value_from_dihedral(value_t d01, value_t d02, value_t d12,
                            value_t d13, value_t d23, value_t omega);

value_t values_to_angle (value_t d01, value_t d02, value_t d12);

value_t values_to_dihedral (value_t d01, value_t d02, value_t d03,
                            value_t d12, value_t d13, value_t d23);

value_t values_to_chord (value_t d01, value_t d02, value_t d03,
                         value_t d12, value_t d13, value_t d23);

double distances_to_angle (double d01, double d02, double d12);

double distances_to_dihedral (double d01, double d02, double d03,
                              double d12, double d13, double d23);

