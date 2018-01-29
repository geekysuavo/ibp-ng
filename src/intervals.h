
/* ensure once-only inclusion. */
#pragma once

/* include the traceback header. */
#include "trace.h"

/* intervals_print(): macro to pretty-print interval sets.
 */
#define intervals_print(I)  intervals_printfn(I, #I)

/* intervals_free(): macro to free allocated interval set structure pointers.
 */
#define intervals_free(I)  { free(I); I = NULL; }

/* intervals_t: structure for holding sets of intervals.
 */
typedef struct {
  /* interval data:
   *  @start: array of interval lower bounds.
   *  @end: array of interval upper bounds.
   */
  double *start, *end;

  /* interval set data:
   *  @size: current number of active intervals in the set.
   *  @capacity: total available number of intervals.
   */
  unsigned int size, capacity;
}
intervals_t;

/* function declarations (intervals.c): */

intervals_t *intervals_new (unsigned int capacity);

int intervals_union (intervals_t *I, double start, double end);

int intervals_intersect (intervals_t *Ia, intervals_t *Ib, intervals_t *Ic);

double intervals_len (intervals_t *I);

void intervals_grid (intervals_t *I, double *samp, unsigned int *n_samp);

void intervals_printfn (intervals_t *I, const char *id);

