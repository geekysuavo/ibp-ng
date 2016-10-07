
/* include the vector header. */
#include "vector.h"

/* vector_new(): allocate a new zero vector.
 *
 * returns:
 *  pointer to a new vector with zero values.
 */
vector_t *vector_new (void) {
  /* allocate a new vector structure pointer. */
  vector_t *v = (vector_t*) malloc(sizeof(vector_t));

  /* check if allocation failed. */
  if (!v) {
    /* raise an exception and return null. */
    raise("unable to allocate vector structure pointer");
    return NULL;
  }

  /* initialize the values and return the new vector. */
  v->x = v->y = v->z = 0.0;
  return v;
}

/* vector_new_with_value(): allocate a new vector with specified values.
 *
 * arguments:
 *  @x, @y, @z: element values of the new vector.
 *
 * returns:
 *  pointer to the newly allocated and initialized vector.
 */
vector_t *vector_new_with_value (double x, double y, double z) {
  /* allocate a new vector structure pointer. */
  vector_t *v = vector_new();
  if (!v)
    return NULL;

  /* set the values. */
  v->x = x;
  v->y = y;
  v->z = z;

  /* return the new vector. */
  return v;
}

/* vector_set(): set the components of a vector to specified values.
 *
 * arguments:
 *  @v: pointer to the vector structure to modify.
 *  @x, @y, @z: values to set.
 */
inline void vector_set (vector_t *v, double x, double y, double z) {
  /* set the components of the vector. */
  v->x = x;
  v->y = y;
  v->z = z;
}

/* vector_dot(): compute the dot product between two vectors.
 *
 * arguments:
 *  @a, @b: vectors to use for the computation.
 *  @z: pointer to the output value.
 */
inline void vector_dot (vector_t *a, vector_t *b, double *z) {
  /* compute the value. */
  *z = a->x * b->x + a->y * b->y + a->z * b->z;
}

/* vector_sqdist(): compute the squared distance between two vectors.
 *
 * arguments:
 *  @a, @b: vectors to use for the computation.
 *
 * returns:
 *  computed squared distance value.
 */
inline double vector_sqdist (vector_t *a, vector_t *b) {
  /* compute the sum of squares. */
  return pow(a->x - b->x, 2.0) +
         pow(a->y - b->y, 2.0) +
         pow(a->z - b->z, 2.0);
}

/* vector_dist(): compute the distance between two vectors.
 *
 * arguments:
 *  @a, @b: vectors to use for the computation.
 *
 * returns:
 *  computed distance value.
 */
inline double vector_dist (vector_t *a, vector_t *b) {
  /* compute the square root of the sum of squares. */
  return sqrt(pow(a->x - b->x, 2.0) +
              pow(a->y - b->y, 2.0) +
              pow(a->z - b->z, 2.0));
}

/* vector_normalize(): normalize a vector to unit length.
 *
 * arguments:
 *  @a: pointer to the vector to modify.
 */
inline void vector_normalize (vector_t *a) {
  /* compute the length of the vector. */
  const double len = sqrt(a->x * a->x + a->y * a->y + a->z * a->z);

  /* scale the vector by its length. */
  a->x /= len;
  a->y /= len;
  a->z /= len;
}

/* vector_cross(): compute the cross product between two vectors.
 *
 * arguments:
 *  @a, @b: input vector pointers for the computation.
 *  @z: pointer to the output vector.
 */
inline void vector_cross (vector_t *a, vector_t *b, vector_t *z) {
  /* compute the result. */
  z->x = a->y * b->z - a->z * b->y;
  z->y = a->z * b->x - a->x * b->z;
  z->z = a->x * b->y - a->y * b->x;
}

/* vector_angle(): compute the angle among three vectors.
 *
 * arguments:
 *  @a, @b, @c: vectors to use for the computation.
 *
 * returns:
 *  computed angle value.
 */
inline double vector_angle (vector_t *a, vector_t *b, vector_t *c) {
  /* declare required variables:
   *  @ba, @bc: vectors from @b to @a and @c, respectively.
   *  @dot: dot product between the two centered vectors.
   */
  vector_t ba, bc;
  double dot;

  /* compute the first centered vector. */
  ba.x = a->x - b->x;
  ba.y = a->y - b->y;
  ba.z = a->z - b->z;

  /* compute the second centered vector. */
  bc.x = c->x - b->x;
  bc.y = c->y - b->y;
  bc.z = c->z - b->z;

  /* compute the dot product of the two vectors after normalization. */
  vector_normalize(&ba);
  vector_normalize(&bc);
  vector_dot(&ba, &bc, &dot);

  /* return the angle. */
  return acos(dot);
}

/* vector_dihedral(): compute the dihedral angle among four vectors.
 *
 * arguments:
 *  @a, @b, @c, @d: vectors to use for the computation.
 *
 * returns:
 *  computed angle value.
 */
inline double vector_dihedral (vector_t *a, vector_t *b,
                               vector_t *c, vector_t *d) {
  /* declare required variables. */
  vector_t b1, b2, b3;
  vector_t n1, n2, m;
  double x, y;

  /* compute the first vector in the three-vector system. */
  b1.x = a->x - b->x;
  b1.y = a->y - b->y;
  b1.z = a->z - b->z;

  /* compute the second vector in the three-vector system. */
  b2.x = b->x - c->x;
  b2.y = b->y - c->y;
  b2.z = b->z - c->z;

  /* compute the third vector in the three-vector system. */
  b3.x = c->x - d->x;
  b3.y = c->y - d->y;
  b3.z = c->z - d->z;

  /* compute the unit normal of the first plane. */
  vector_cross(&b1, &b2, &n1);
  vector_normalize(&n1);

  /* compute the unit normal of the second plane. */
  vector_cross(&b2, &b3, &n2);
  vector_normalize(&n2);

  /* compute the final vector of the orthonormal frame. */
  vector_normalize(&b2);
  vector_cross(&n1, &b2, &m);

  /* finally, compute the dihedral angle. */
  vector_dot(&n1, &n2, &x);
  vector_dot(&m, &n2, &y);

  /* return the angle. */
  return atan2(y, x);
}

