
/* include the interval set header. */
#include "intervals.h"

/* intervals_bytes(): compute the storage requirements of an interval set.
 *
 * arguments:
 *  @capacity: desired number of intervals in the set.
 *
 * returns:
 *  number of bytes required for an interval set having
 *  the specified capacity.
 */
size_t intervals_bytes (unsigned int capacity) {
  /* compute and return the byte count. */
  return sizeof(intervals_t) + 2 * capacity * sizeof(double);
}

/* intervals_new(): allocate a new empty interval set.
 *
 * arguments:
 *  @capacity: maximum number of intervals to be stored in the set.
 *
 * returns:
 *  pointer to a new interval set that contains no intervals.
 */
intervals_t *intervals_new (unsigned int capacity) {
  /* check that a valid capacity was requested. */
  if (capacity == 0) {
    /* nope. raise an exception and return null. */
    raise("attempted to allocate an empty interval set");
    return NULL;
  }

  /* allocate a new interval set pointer of the required size. */
  intervals_t *I = malloc(intervals_bytes(capacity));
  if (!I)
    return NULL;

  /* initialize the array of interval lower bound values. */
  char *ptr = ((char*) I) + sizeof(intervals_t);
  I->start = (double*) ptr;

  /* initialize the array of interval upper bound values. */
  ptr += capacity * sizeof(double);
  I->end = (double*) ptr;

  /* initialize the capacity and size of the interval set. */
  I->capacity = capacity;
  I->size = 0;

  /* return the new structure pointer. */
  return I;
}

/* intervals_compress(): combine all overlapping intervals contained
 * within an interval set.
 *
 * arguments:
 *  @I: pointer to the interval set to modify.
 */
static void intervals_compress (intervals_t *I) {
  /* declare required variables:
   *  @i: index of the intervals in the set.
   *  @j: index for storing combined intervals.
   */
  unsigned int i, j;

  for (i = 1, j = 0; i < I->size; i++) {
    /* check for non-empty intersections of successive intervals. */
    if (I->start[i] <= I->end[j]) {
      /* end[j] = max{ end[i], end[j] }. */
      if (I->end[j] < I->end[i])
        I->end[j] = I->end[i];
    }
    else {
      /* increment the storage index. */
      j++;

      /* copy the new interval into the last available position. */
      if (j < i) {
        I->start[j] = I->start[i];
        I->end[j] = I->end[i];
      }
    }
  }

  /* update the size of the compressed set. */
  I->size = j + 1;
}

/* intervals_union: compute the union of a set of intervals with
 * an additional interval.
 *
 * arguments:
 *  @I: pointer to the interval set to modify in-place.
 *  @start: lower bound of the additional interval.
 *  @end: upper bound of the additional interval.
 *
 * returns:
 *  integer indicating whether (1) or not (0) the operation succeeded.
 */
int intervals_union (intervals_t *I, double start, double end) {
  /* check if the interval set can admit a new interval. */
  if (I->size < I->capacity) {
    /* locate the array index where the new interval will be added. */
    unsigned int i = 0;
    while (i < I->size && I->start[i] < start) i++;

    /* shift all remaining intervals down one position in the array. */
    for (unsigned int j = I->size; j > i; j--) {
      I->start[j] = I->start[j - 1];
      I->end[j] = I->end[j - 1];
    }

    /* store the new interval bounds in the array. */
    I->start[i] = start;
    I->end[i] = end;

    /* if the interval overlaps with its neighbors in the array, then
     * the interval set needs to be compressed back into a disjoint
     * set of intervals.
     */
    const unsigned int last = I->size++;
    if ((i > 0 && I->start[i] <= I->end[i - 1]) ||
        (i < last && I->start[i + 1] <= I->end[i]))
      intervals_compress(I);
  }
  else {
    /* throw an exception. */
    throw("interval capacity (%u) exceeded in union", I->capacity);
  }

  /* return success. */
  return 1;
}

/* intervals_intersect(): compute the intersection of two interval sets,
 * which is itself an interval set.
 *
 * arguments:
 *  @Ia: interval set structure pointer, first input operand.
 *  @Ib: interval set structure pointer, second input operand.
 *  @Ic: interval set structure pointer, output.
 *
 * returns:
 *  integer indicating whether (1) or not (0) the operation succeeded.
 */
int intervals_intersect (intervals_t *Ia, intervals_t *Ib, intervals_t *Ic) {
  /* empty the output set. */
  Ic->size = 0;

  /* loop through the first input interval set. */
  for (unsigned int i = 0; i < Ia->size; i++) {
    /* get the lower and upper bounds of the current interval in @Ia. */
    const double la = Ia->start[i];
    const double ua = Ia->end[i];

    /* loop through the second input interval set. */
    for (unsigned int j = 0; j < Ib->size; j++) {
      /* get the lower and upper bounds of the current interval in @Ib. */
      double lb = Ib->start[j];
      double ub = Ib->end[j];

      /* uses the fact that both input sets are sorted. */
      if (lb > ua)
        break;

      /* compute the intersection of the current pair of intervals. */
      if (la > lb) lb = la;
      if (ua < ub) ub = ua;

      /* if the intersection is non-empty, it should be added to the
       * output interval set.
       */
      if (lb <= ub) {
        /* throw an exception if the output interval set is full. */
        if (Ic->size == Ic->capacity)
          throw("interval capacity (%u) exceeded in intersection",
                Ic->capacity);

        /* store the intersection result into the output interval set. */
        Ic->start[Ic->size] = lb;
        Ic->end[Ic->size] = ub;
        Ic->size++;
      }
    }
  }

  /* return success. */
  return 1;
}

/* intervals_grid(): uniformly discretize an interval set with a specified
 * number of points, where the points are equally spaced on the interval
 * formed by concatenating all intervals in the set together end-to-end.
 *
 * when the interval set consists of a union of degenerate intervals
 * (i.e. exact values), the resulting sample set will consist of
 * each of the interval values.
 *
 * note:
 *  this function does not handle the case when the interval set contains
 *  both proper and degenerate intervals, as this situation never arises
 *  within the mathematics of enum_reduce().
 *
 * arguments:
 *  @I: pointer to the interval set to access.
 *  @samp: output array of sampled values.
 *  @n_samp: pointer to both:
 *            - the desired number of values to sample.
 *            - the final number of sampled values.
 */
void intervals_grid (intervals_t *I, double *samp, unsigned int *n_samp) {
  /* the final number of sampled points will be the lesser of the
   * number of available and requested point counts.
   */
  const unsigned int ns = (I->size < *n_samp ? I->size : *n_samp);

  /* immediately finish if the interval set is empty. */
  if (ns == 0)
    goto final;

  /* compute the combined length of all intervals in the set. */
  double len = 0.0;
  for (unsigned int i = 0; i < I->size; i++)
    len += I->end[i] - I->start[i];

  /* check if the intervals are degenerate. */
  if (len == 0) {
    /* store the lower bounds (values) of each interval (point)
     * into the output samples array, until we have exhausted
     * either the interval set or the array size.
     */
    for (unsigned int i = 0; i < ns; i++)
      samp[i] = I->start[i];
  }
  else {
    /* compute the spacing between sampled points. */
    const double h = len / (ns - 1);

    /* initialize the first sampled point and zero the length. */
    samp[0] = I->start[0];
    len = 0.0;

    /* loop until all sampled points are computed. */
    for (unsigned int i = 0, j = 0; i < ns - 1; /**/) {
      /* compute the distance from the last sampled point to the end
       * of the currently indexed interval.
       */
      const double start = (I->start[j] > samp[i] ? I->start[j] : samp[i]);
      const double end = I->end[j];
      const double dl = end - start;

      /* check if the distance is sufficient. */
      if (len + dl >= (h - 1.0e-8)) {
        /* yes. add the new sampled point and zero the length. */
        samp[++i] = start + (h - len);
        len = 0.0;
      }
      else {
        /* not yet. accumulate the length and move to the next interval. */
        len += dl;
        j++;
      }
    }
  }

final:
  /* store the final sample count. */
  *n_samp = ns;
}

/* intervals_printfn(): core function used by the intervals_print() macro
 * to write interval sets to standard output.
 *
 * arguments:
 *  @I: pointer to the interval set to access.
 *  @id: string identifier of the interval.
 */
void intervals_printfn (intervals_t *I, const char *id) {
  /* print the variable name. */
  printf("%s (%u/%u) = ", id, I->size, I->capacity);

  /* loop over the intervals in the set. */
  for (unsigned int i = 0; i < I->size; i++)
    printf("[%lf, %lf]%s", I->start[i], I->end[i],
           i + 1 == I->size ? "" : ", ");

  /* write a newline. */
  printf("\n");
}

