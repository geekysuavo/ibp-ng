
/* include the value header. */
#include "value.h"

/* value_undefined(): construct an undefined generalized value.
 *
 * returns:
 *  a value_t with content understood to be undefined.
 */
value_t value_undefined (void) {
  /* declare required variables:
   *  @val: the output value.
   */
  value_t val;

  /* set the value type. */
  val.type = VALUE_TYPE_UNDEFINED;

  /* set the value bounds. */
  val.l = val.u = NAN;

  /* return the value. */
  return val;
}

/* value_scalar(): construct a scalar generalized value.
 *
 * arguments:
 *  @v: exact numerical content of the value.
 *
 * returns:
 *  a value_t with content understood to be a scalar.
 */
value_t value_scalar (double v) {
  /* declare required variables:
   *  @val: the output value.
   */
  value_t val;

  /* set the value type. */
  val.type = VALUE_TYPE_SCALAR;

  /* set the value bounds. */
  val.l = val.u = v;

  /* return the value. */
  return val;
}

/* value_interval(): construct an interval generalized value.
 *
 * if the provided interval bounds are equal, then this function
 * will automatically convert the return value into a scalar.
 *
 * arguments:
 *  @l: lower bound of the value.
 *  @u: upper bound of the value.
 *
 * returns:
 *  a value_t with content understood to be an interval.
 */
value_t value_interval (double l, double u) {
  /* declare required variables:
   *  @val: the output value.
   */
  value_t val;

  /* return undefined if the interval bounds are invalid. */
  if (l > u)
    return value_undefined();

  /* return scalar if the interval bounds are equal. */
  if (l == u)
    return value_scalar(l);

  /* set the value type. */
  val.type = VALUE_TYPE_INTERVAL;

  /* set the value bounds. */
  val.l = l;
  val.u = u;

  /* return the value. */
  return val;
}

/* value_is_undefined(): return whether a generalized value is undefined.
 *
 * arguments:
 *  @val: the value to check.
 *
 * returns:
 *  integer indicating whether (1) or not (0) the value is undefined.
 */
int value_is_undefined (value_t val) {
  /* return whether the value is defined or not. */
  return (val.type == VALUE_TYPE_UNDEFINED);
}

/* value_is_scalar(): return whether a generalized value is a scalar.
 *
 * arguments:
 *  @val: the value to check.
 *
 * returns:
 *  integer indicating whether (1) or not (0) the value is a scalar.
 */
int value_is_scalar (value_t val) {
  /* return whether the value is a scalar or not. */
  return (val.type == VALUE_TYPE_SCALAR);
}

/* value_is_interval(): return whether a generalized value is an interval.
 *
 * arguments:
 *  @val: the value to check.
 *
 * returns:
 *  integer indicating whether (1) or not (0) the value is an interval.
 */
int value_is_interval (value_t val) {
  /* return whether the value is an interval or not. */
  return (val.type == VALUE_TYPE_INTERVAL);
}

/* value_add(): compute the sum of two generalized values.
 *
 * arguments:
 *  @va: first value in the operation.
 *  @vb: second value in the operation.
 *
 * returns:
 *  value containing the sum of the two input values.
 */
value_t value_add (value_t va, value_t vb) {
  /* check if either value is undefined. */
  if (va.type == VALUE_TYPE_UNDEFINED ||
      vb.type == VALUE_TYPE_UNDEFINED)
    return value_undefined();

  /* compute and return the result. */
  return value_interval(va.l + vb.l, va.u + vb.u);
}

/* value_sub(): compute the difference of two generalized values.
 *
 * arguments:
 *  @va: first value in the operation.
 *  @vb: second value in the operation.
 *
 * returns:
 *  value containing the difference of the two input values.
 */
value_t value_sub (value_t va, value_t vb) {
  /* check if either value is undefined. */
  if (va.type == VALUE_TYPE_UNDEFINED ||
      vb.type == VALUE_TYPE_UNDEFINED)
    return value_undefined();

  /* compute and return the result. */
  return value_interval(va.l - vb.u, va.u - vb.l);
}

/* value_mul(): compute the product of two generalized values.
 *
 * arguments:
 *  @va: first value in the operation.
 *  @vb: second value in the operation.
 *
 * returns:
 *  value containing the product of the two input values.
 */
value_t value_mul (value_t va, value_t vb) {
  /* declare required variables:
   *  @l, @u: final bounds of the computed value.
   *  @ll, @lu, @ul, @uu: intermediate bounds of the value.
   */
  double l, u, ll, lu, ul, uu;

  /* check if either value is undefined. */
  if (va.type == VALUE_TYPE_UNDEFINED ||
      vb.type == VALUE_TYPE_UNDEFINED)
    return value_undefined();

  /* compute all possible intermediate bounds on the value. */
  ll = va.l * vb.l;
  lu = va.l * vb.u;
  ul = va.u * vb.l;
  uu = va.u * vb.u;

  /* determine the lower bound of the new value. */
  l = ll;
  if (lu < l) l = lu;
  if (ul < l) l = ul;
  if (uu < l) l = uu;

  /* determine the upper bound of the new value. */
  u = ll;
  if (lu > u) u = lu;
  if (ul > u) u = ul;
  if (uu > u) u = uu;

  /* return the final value. */
  return value_interval(l, u);
}

/* value_div(): compute the quotient of two generalized values.
 *
 * arguments:
 *  @va: first value in the operation.
 *  @vb: second value in the operation.
 *
 * returns:
 *  value containing the quotient of the two input values.
 */
value_t value_div (value_t va, value_t vb) {
  /* declare required variables:
   *  @l, @u: final bounds of the computed value.
   *  @ll, @lu, @ul, @uu: intermediate bounds of the value.
   */
  double l, u, ll, lu, ul, uu;

  /* check if either value is undefined. */
  if (va.type == VALUE_TYPE_UNDEFINED ||
      vb.type == VALUE_TYPE_UNDEFINED)
    return value_undefined();

  /* compute all possible intermediate bounds on the value. */
  ll = va.l / vb.l;
  lu = va.l / vb.u;
  ul = va.u / vb.l;
  uu = va.u / vb.u;

  /* determine the lower bound of the new value. */
  l = ll;
  if (lu < l) l = lu;
  if (ul < l) l = ul;
  if (uu < l) l = uu;

  /* determine the upper bound of the new value. */
  u = ll;
  if (lu > u) u = lu;
  if (ul > u) u = ul;
  if (uu > u) u = uu;

  /* return the final value. */
  return value_interval(l, u);
}

/* value_pow(): compute the exponentiation of a generalized value.
 *
 * arguments:
 *  @v: base value in the operation.
 *  @p: scalar exponent in the operation.
 *
 * returns:
 *  value containing the exponentiation of the input values.
 */
value_t value_pow (value_t v, double p) {
  /* declare required variables:
   *  @l, @u: final bounds of the computed value.
   *  @s: intermediate bound swap variable.
   */
  double l, u, s;

  /* check if the input value is undefined. */
  if (v.type == VALUE_TYPE_UNDEFINED)
    return v;

  /* compute all possible intermediate bounds on the value. */
  l = pow(v.l, p);
  u = pow(v.u, p);

  /* check if the bounds are out of order. */
  if (u < l) {
    /* yes. swap the bounds. */
    s = l;
    l = u;
    u = s;
  }

  /* return the final value. */
  return value_interval(l, u);
}

/* value_scal(): compute the scalar multiple of a generalized value.
 *
 * arguments:
 *  @v: first value in the operation.
 *  @p: scalar factor in the operation.
 *
 * returns:
 *  value containing the product of the input values.
 */
value_t value_scal (value_t v, double p) {
  /* declare required variables:
   *  @l, @u: final bounds of the computed value.
   *  @s: intermediate bound swap variable.
   */
  double l, u, s;

  /* check if the input value is undefined. */
  if (v.type == VALUE_TYPE_UNDEFINED)
    return v;

  /* compute all possible intermediate bounds on the value. */
  l = v.l * p;
  u = v.u * p;

  /* check if the bounds are out of order. */
  if (u < l) {
    /* yes. swap the bounds. */
    s = l;
    l = u;
    u = s;
  }

  /* return the final value. */
  return value_interval(l, u);
}

/* value_bound(): bound a generalized value within another.
 *
 * arguments:
 *  @v: generalized value to bound.
 *  @b: bound to apply.
 *
 * returns:
 *  bounded generalized value.
 */
value_t value_bound (value_t v, value_t b) {
  /* declare required variables:
   *  @l, @u: final bounds of the computed value.
   */
  double l, u;

  /* check if either value is undefined. */
  if (v.type == VALUE_TYPE_UNDEFINED || b.type == VALUE_TYPE_UNDEFINED)
    return v;

  /* bound the bounds. :) */
  l = (v.l < b.l ? b.l : v.l > b.u ? b.u : v.l);
  u = (v.u < b.l ? b.l : v.u > b.u ? b.u : v.u);

  /* return the final value. */
  return value_interval(l, u);
}

/* value_sin(): compute the sine of a generalized value.
 *
 * arguments:
 *  @v: input value in the operation.
 *
 * returns:
 *  value containing the sine of the input value.
 */
value_t value_sin (value_t v) {
  /* declare required variables:
   *  @l, @u: final bounds of the computed value.
   *  @s: intermediate bound swap variable.
   *  @n, @nl, @nu: sine critical points.
   */
  double l, u, s;
  int n, nl, nu;

  /* check if the input value is undefined. */
  if (v.type == VALUE_TYPE_UNDEFINED)
    return v;

  /* compute all possible intermediate bounds on the value. */
  l = sin(v.l);
  u = sin(v.u);

  /* check if the bounds are out of order. */
  if (u < l) {
    /* yes. swap the bounds. */
    s = l;
    l = u;
    u = s;
  }

  /* bound the set of critical points to explicitly check. */
  nl = (int) ceil(v.l / M_PI - 0.5);
  nu = (int) floor(v.u / M_PI - 0.5);

  /* loop over the sine critical points. */
  for (n = nl; n <= nu; n++) {
    /* compute the value at the critical point. */
    s = sin(M_PI * (double) n + M_PI_2);

    /* update the bounds accordingly. */
    if (s < l) l = s;
    if (s > u) u = s;
  }

  /* return the final value. */
  return value_interval(l, u);
}

/* value_cos(): compute the cosine of a generalized value.
 *
 * arguments:
 *  @v: input value in the operation.
 *
 * returns:
 *  value containing the cosine of the input value.
 */
value_t value_cos (value_t v) {
  /* declare required variables:
   *  @l, @u: final bounds of the computed value.
   *  @s: intermediate bound swap variable.
   *  @n, @nl, @nu: cosine critical points.
   */
  double l, u, s;
  int n, nl, nu;

  /* check if the input value is undefined. */
  if (v.type == VALUE_TYPE_UNDEFINED)
    return v;

  /* compute all possible intermediate bounds on the value. */
  l = cos(v.l);
  u = cos(v.u);

  /* check if the bounds are out of order. */
  if (u < l) {
    /* yes. swap the bounds. */
    s = l;
    l = u;
    u = s;
  }

  /* bound the set of critical points to explicitly check. */
  nl = (int) ceil(v.l / M_PI);
  nu = (int) floor(v.u / M_PI);

  /* loop over the cosine critical points. */
  for (n = nl; n <= nu; n++) {
    /* compute the value at the critical point. */
    s = cos(M_PI * (double) n);

    /* update the bounds accordingly. */
    if (s < l) l = s;
    if (s > u) u = s;
  }

  /* return the final value. */
  return value_interval(l, u);
}

/* value_printfn(): core function used by the value_print() macro to
 * output generalized values to standard output.
 *
 * arguments:
 *  @v: value to print.
 *  @id: string identifier of the value.
 */
inline void value_printfn (value_t v, const char *id) {
  /* print the variable name. */
  printf("%s = ", id);

  /* check the variable type. */
  if (v.type == VALUE_TYPE_SCALAR) {
    /* print the scalar. */
    printf("%lf\n", v.l);
  }
  else if (v.type == VALUE_TYPE_INTERVAL) {
    /* print the interval. */
    printf("[%lf, %lf]\n", v.l, v.u);
  }
  else {
    /* print an undefined. */
    printf("undefined\n");
  }
}

/* value_from_angle(): compute the distance between the terminal points
 * of a chain of three points, where the distance between the first
 * and second point, the distance between the second and third point,
 * and the angle formed by the three points, are all known.
 *
 * equation:
 *  d02 = sqrt(d01^2 + d12^2 - 2*d01*d12*cos(theta))
 *
 * arguments:
 *  @d01: distance between the first and second point.
 *  @d12: distance between the second and third point.
 *  @theta: angle formed by the points, in radians.
 *
 * returns:
 *  distance between the first and third point.
 */
value_t value_from_angle (value_t d01, value_t d12,
                          value_t theta) {
  /* declare required variables:
   *  @x: value containing d01^2 + b^2.
   *  @y: value containing 2*d01*b.
   *  @z: value containing cos(theta).
   */
  value_t x, y, z;

  /* compute the intermediate values. */
  x = value_add(value_pow(d01, 2.0), value_pow(d12, 2.0));
  y = value_scal(value_mul(d01, d12), 2.0);
  z = value_cos(theta);

  /* compute the final value. */
  return value_pow(value_sub(x, value_mul(y, z)), 0.5);
}

/* value_from_dihedral(): compute the distance between the terminal
 * points of a chain of four points, where all other distances between
 * each pair of points is known.
 *
 * arguments:
 *  @d01: distance between the first and second point.
 *  @d02: distance between the first and third point.
 *  @d12: distance between the second and third point.
 *  @d13: distance between the second and fourth point.
 *  @d23: distance between the third and fourth point.
 *  @omega: dihedral angle formed by the points, in radians.
 *
 * returns:
 *  distance between the first and fourth point.
 */
value_t value_from_dihedral(value_t d01, value_t d02, value_t d12,
                            value_t d13, value_t d23, value_t omega) {
  /* declare required variables:
   *  @s??: squared values of each respective @d??.
   *  @x, @y, @z, @xy, @xn, @yn, @xyn: intermediates.
   */
  value_t s01, s02, s12, s13, s23;
  value_t x, y, z, xy, xn, yn, xyn;

  /* square the input distances. */
  s01 = value_pow(d01, 2.0); /* s01 = d01^2 */
  s02 = value_pow(d02, 2.0); /* s02 = d02^2 */
  s12 = value_pow(d12, 2.0); /* s12 = d12^2 */
  s13 = value_pow(d13, 2.0); /* s13 = d13^2 */
  s23 = value_pow(d23, 2.0); /* s23 = d23^2 */

  /* x = 0.5 * (s13 + s12 - s23) / (d13 * d12) */
  x = value_scal(value_div(value_sub(value_add(s13, s12), s23),
                           value_mul(d13, d12)), 0.5);

  /* y = 0.5 * (s01 + s12 - s02) / (d01 * d12) */
  y = value_scal(value_div(value_sub(value_add(s01, s12), s02),
                           value_mul(d01, d12)), 0.5);

  /* xn = 1.0 - x^2
   * yn = 1.0 - y^2
   */
  xn = value_sub(value_scalar(1.0), value_pow(x, 2.0));
  yn = value_sub(value_scalar(1.0), value_pow(y, 2.0));

  /* xy = x * y
   * xyn = xn * yn
   */
  xy = value_mul(x, y);
  xyn = value_mul(xn, yn);

  /* z = 2.0 * (sqrt(xyn) * cos(omega) + xy) * d01 * d13 */
  z = value_mul(value_scal(value_add(value_mul(value_pow(xyn, 0.5),
                                               value_cos(omega)),
                                     xy), 2.0),
                value_mul(d01, d13));

  /* sqrt(s01 + s13 - z) */
  return value_pow(value_sub(value_add(s01, s13), z), 0.5);
}

/* values_to_angle(): compute the cosine of an angle from the
 * relevant *generalized* distances between all three vertices.
 *
 * arguments:
 *  @dij: distances between vertices @i and @j.
 *
 * returns:
 *  computed cosine of the angle.
 */
value_t values_to_angle (value_t d01, value_t d02, value_t d12) {
  /* x = d12 * d12 + d01 * d01 - d02 * d02 */
  value_t x = value_sub(value_add(value_mul(d12, d12),
                                  value_mul(d01, d01)),
                        value_mul(d02, d02));

  /* y = 2.0 * d12 * d01 */
  value_t y = value_scal(value_mul(d12, d01), 2.0);

  /* cos(theta) = x / y */
  return value_div(x, y);
}

/* values_to_dihedral(): compute the cosine of a dihedral angle
 * from the relevant *generalized* distances between all four
 * vertices.
 *
 * arguments:
 *  @dij: distances between vertices @i and @j.
 *
 * returns:
 *  computed cosine of the dihedral angle. no bounds checking is performed.
 */
value_t values_to_dihedral (value_t d01, value_t d02, value_t d03,
                            value_t d12, value_t d13, value_t d23) {
  /* handle the special case of x0 == x3 */
  if (d03.l == 0.0 && d03.u == 0.0)
    return value_scalar(1.0);

  /* a = 0.5 * (d01*d01 + d13*d13 - d03*d03) / (d01*d13) */
  value_t a = value_scal(value_div(value_sub(value_add(value_mul(d01, d01),
                                                       value_mul(d13, d13)),
                                             value_mul(d03, d03)),
                                   value_mul(d01, d13)), 0.5);

  /* b = 0.5 * (d13*d13 + d12*d12 - d23*d23) / (d13*d12) */
  value_t b = value_scal(value_div(value_sub(value_add(value_mul(d13, d13),
                                                       value_mul(d12, d12)),
                                             value_mul(d23, d23)),
                                   value_mul(d13, d12)), 0.5);

  /* c = 0.5 * (d01*d01 + d12*d12 - d02*d02) / (d01*d12) */
  value_t c = value_scal(value_div(value_sub(value_add(value_mul(d01, d01),
                                                       value_mul(d12, d12)),
                                             value_mul(d02, d02)),
                                   value_mul(d01, d12)), 0.5);

  /* e = 1.0 - b * b
   * f = 1.0 - c * c
   */
  value_t e = value_sub(value_scalar(1.0), value_mul(b, b));
  value_t f = value_sub(value_scalar(1.0), value_mul(c, c));

  /* cos(omega) = (a - b * c) / sqrt(e * f) */
  return value_div(value_sub(a, value_mul(b, c)),
                   value_pow(value_mul(e, f), 0.5));
}

/* values_to_chord(): compute the chords of a given iBP solution set for
 * one node in a tree. used for RMSD pruning.
 *
 * arguments:
 *  @dij: distances between vertices @i and @j.
 *
 * returns:
 *  lengths of the min-chord and max-chord on the circle spanned by the
 *  dihedral angles of the set of distances.
 */
value_t values_to_chord (value_t d01, value_t d02, value_t d03,
                         value_t d12, value_t d13, value_t d23) {
  /* compute cos(theta) and cos(omega) */
  value_t ct = values_to_angle(d12, d13, d23);
  value_t cw = values_to_dihedral(d01, d02, d03, d12, d13, d23);

  /* bound the dihedral's cosine value. */
  cw.l = (cw.l > 1.0 ? 1.0 : cw.l < -1.0 ? -1.0 : cw.l);
  cw.u = (cw.u > 1.0 ? 1.0 : cw.u < -1.0 ? -1.0 : cw.u);

  /* sin(theta) = sqrt(1 - cos^2(theta))
   * sin(omega) = sqrt(1 - cos^2(omega))
   */
  value_t st = value_pow(value_sub(value_scalar(1.0),
                                   value_mul(ct, ct)), 0.5);
  value_t sw = value_pow(value_sub(value_scalar(1.0),
                                   value_mul(cw, cw)), 0.5);

  /* c = 2.0 * d23 * sin(theta) * sin(omega) */
  return value_scal(value_mul(d23, value_mul(st, sw)), 2.0);
}

/* distances_to_angle(): compute the cosine of an angle from the
 * relevant *scalar* distances between all three vertices.
 *
 * arguments:
 *  @dij: distances between vertices @i and @j.
 *
 * returns:
 *  computed cosine of the angle.
 */
inline double distances_to_angle (double d01, double d02, double d12) {
  /* compute and return the cosine of the angle. */
  return (d12 * d12 + d01 * d01 - d02 * d02) / (2.0 * d12 * d01);
}

/* distances_to_dihedral(): compute the cosine of a dihedral angle
 * from the relevant *scalar* distances between all four vertices.
 *
 * arguments:
 *  @dij: distances between vertices @i and @j.
 *
 * returns:
 *  computed cosine of the dihedral angle. no bounds checking is performed.
 */
double distances_to_dihedral (double d01, double d02, double d03,
                              double d12, double d13, double d23) {
  /* handle the special case of x0 == x3. */
  if (d03 == 0.0)
    return 1.0;

  /* compute a first set of temporary values. */
  const double a = 0.5 * (d01*d01 + d13*d13 - d03*d03) / (d01*d13);
  const double b = 0.5 * (d13*d13 + d12*d12 - d23*d23) / (d13*d12);
  const double c = 0.5 * (d01*d01 + d12*d12 - d02*d02) / (d01*d12);

  /* compute a second set of temporary values. */
  const double e = 1.0 - b*b;
  const double f = 1.0 - c*c;

  /* compute and return the cosine of the dihedral. */
  return (a - b * c) / sqrt(e * f);
}

