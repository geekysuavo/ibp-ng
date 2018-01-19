
/* include the enumerator header. */
#include "enum.h"
#include "enum-thread.h"
#include "enum-prune.h"

/* include the c float header. */
#include <float.h>

/* enum_prune_future_t: structure for holding information required for
 * future distance feasibility pruning closures.
 */
typedef struct {
  /* @ntest, @nprune: test and prune counts of the current closure.
   * @i: graph level of the upstream (embedded) atom.
   * @j: graph level of the test (current) atom.
   * @k: graph level of the downstream (future) atom.
   * @limit: upper bound on d(xi,xj).
   */
  unsigned int ntest, nprune;
  unsigned int i, j, k;
  double limit;
}
enum_prune_future_t;

/* enum_prune_future_init(): initialize the future feasibility pruner.
 */
int enum_prune_future_init (enum_t *E, unsigned int lev) {
  /* declare required variables:
   *  @i, @j, @k: upstream, current and future graph levels.
   *  @klim: future graph level of the least upper bound.
   *  @dik: outer graph edge of the triangle.
   *  @djk: second inner graph edge of the triangle.
   *  @lim: least upper bound identified so far on d(xi,xj).
   */
  unsigned int i, j, k, klim;
  enum_prune_future_t *data;
  value_t dik, djk;
  double lim;

  /* do not test future feasibility of the initial clique. */
  if (lev < 3)
    return 1;

  /* locally store the re-order and originality arrays. */
  const unsigned int *order = E->G->order;
  const unsigned int *dup = E->G->orig;
  const unsigned int n = E->G->n_order;

  /* initialize the graph level constants. */
  j = lev;

  /* loop over upstream atoms in the order. */
  for (i = 0; i < j; i++) {
    /* initialize the limit values to correspond to the case
     * where no upper bound exists.
     */
    lim = DBL_MAX;
    klim = i;

    /* loop over downstream atoms in the order. */
    for (k = j + 1; k < n; k++) {
      /* skip duplicate atoms. */
      if (dup[i] || dup[k])
        continue;

      /* get the required graph edges used for testing feasibility. */
      dik = graph_get_edge(E->G, order[i], order[k]);
      djk = graph_get_edge(E->G, order[j], order[k]);

      /* skip atom sets without defined graph edges. */
      if (dik.type == VALUE_TYPE_UNDEFINED ||
          djk.type == VALUE_TYPE_UNDEFINED)
        continue;

      /* check if the new bound is less than the current bound. */
      if (dik.u + djk.u < lim) {
        /* yes, it is: replace the bound. */
        lim = dik.u + djk.u;
        klim = k;
      }
    }

    /* if no bound was identified, move to the next upstream atom. */
    if (klim == i)
      continue;

    /* allocate a new closure payload. */
    data = (enum_prune_future_t*) malloc(sizeof(enum_prune_future_t));
    if (!data)
      return 0;

    /* initialize the payload contents. */
    data->ntest = data->nprune = 0;
    data->i = order[i];
    data->j = order[j];
    data->k = order[klim];
    data->limit = lim;

    /* register the closure with the enumerator. */
    if (!enum_prune_add_closure(E, lev, enum_prune_future, data))
      return 0;
  }

  /* return success. */
  return 1;
}

/* enum_prune_future(): determine whether an enumerator tree may be pruned
 * at a given node based on future distance feasibility.
 */
int enum_prune_future (enum_t *E, enum_thread_t *th, void *data) {
  /* get the closure payload. */
  enum_prune_future_t *future_data = (enum_prune_future_t*) data;

  /* extract pretty handles to the atom positions. */
  vector_t xi = th->state[future_data->i].pos;
  vector_t xj = th->state[future_data->j].pos;

  /* compute the current distance. */
  const double dij = vector_dist(&xi, &xj);

  /* check if the computed distance obeys the bounds. */
  future_data->ntest++;
  if (dij > future_data->limit + E->ddf_tol) {
    /* out of bounds; prune. */
    future_data->nprune++;
    return 1;
  }

  /* do not prune. */
  return 0;
}

/* enum_prune_future_report(): output a report for the future distance
 * feasibility pruning closure.
 */
void enum_prune_future_report (enum_t *E, unsigned int lev, void *data) {
  /* get the closure payload. */
  enum_prune_future_t *future_data = (enum_prune_future_t*) data;

  /* return if no prunes were performed by the closure. */
  if (!future_data->nprune) return;

  /* get the atom indices. */
  unsigned int a0 = future_data->i;
  unsigned int a1 = future_data->j;
  unsigned int a2 = future_data->k;

  /* get the residue indices. */
  unsigned int r0 = E->P->atoms[a0].res_id;
  unsigned int r1 = E->P->atoms[a1].res_id;
  unsigned int r2 = E->P->atoms[a2].res_id;

  /* get the atom names. */
  const char *atom0 = E->P->atoms[a0].name;
  const char *atom1 = E->P->atoms[a1].name;
  const char *atom2 = E->P->atoms[a2].name;

  /* get the residue codes. */
  const char res0 = peptide_get_reschar(E->P, r0);
  const char res1 = peptide_get_reschar(E->P, r1);
  const char res2 = peptide_get_reschar(E->P, r2);

  /* compute the percentage. */
  double f = ((double) future_data->nprune) /
             ((double) future_data->ntest) * 100.0;

  /* output the statistics. */
  printf("  %c%-4u %-4s | %c%-4u %-4s | %c%-4u %-4s : "
         "%16u/%-16u  %6.2lf%%\n",
         res0, r0 + 1, atom0,
         res1, r1 + 1, atom1,
         res2, r2 + 1, atom2,
         future_data->nprune, future_data->ntest, f);
}

