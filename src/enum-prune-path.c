
/* include the enumerator header. */
#include "enum.h"
#include "enum-thread.h"
#include "enum-prune.h"

/* enum_prune_path_t: structure for holding information required for
 * shortest path feasibility pruning closures.
 */
typedef struct {
  /* @ntest: number of calls made or each neighbor pair.
   * @nprune: number of prunes performed for each neighbor pair.
   */
  unsigned int **ntest, **nprune;
}
enum_prune_path_t;

/* enum_prune_path_init(): initialize the path pruning device.
 */
int enum_prune_path_init (enum_t *E, unsigned int lev) {
  /* declare required variables:
   *  @data: closure data pointer.
   *  @sz: closure memory footprint.
   *  @na, @nb: row and column counts.
   */
  enum_prune_path_t *data;
  unsigned int sz, na, nb;
  char *rtest, *rprune;

  /* compute the number of table rows and columns. */
  na = lev;
  nb = E->G->n_order - lev - 1;

  /* compute the closure memory block size. */
  sz = sizeof(enum_prune_path_t)
     + 2 * na * sizeof(unsigned int*)
     + 2 * na * nb * sizeof(unsigned int);

  /* allocate the closure payload. */
  data = (enum_prune_path_t*) malloc(sz);
  if (!data)
    return 0;

  /* initialize the first array pointer. */
  sz = sizeof(enum_prune_path_t);
  data->ntest = (unsigned int**) (((char*) data) + sz);

  /* initialize the second array pointer. */
  sz += na * (nb * sizeof(unsigned int) + sizeof(unsigned int*));
  data->nprune = (unsigned int**) (((char*) data) + sz);

  /* initialize the array row pointers. */
  rtest  = ((char*) data->ntest)  + (na * sizeof(unsigned int*));
  rprune = ((char*) data->nprune) + (na * sizeof(unsigned int*));
  sz  = nb * sizeof(unsigned int);
  for (unsigned int i = 0; i < na; i++) {
    /* set the row pointers. */
    data->ntest[i]  = (unsigned int*) (rtest  + (i * sz));
    data->nprune[i] = (unsigned int*) (rprune + (i * sz));

    /* initialize the row contents. */
    for (unsigned int k = 0; k < nb; k++)
      data->ntest[i][k] = data->nprune[i][k] = 0;
  }

  /* register a closure. */
  if (!enum_prune_add_closure(E, lev, enum_prune_path, data))
    return 0;

  /* return success. */
  return 1;
}

/* enum_prune_path(): determine whether an enumerator tree may be pruned
 * at a given node based on shortest path feasibility.
 */
int enum_prune_path (enum_t *E, enum_thread_t *th, void *data) {
  /* declare required variables:
   *  @dik: distance bound from x(i) to x(k).
   *  @djk: distance bound from x(j) to x(k).
   *  @dij: current distance between x(i) and x(j).
   */
  value_t dik, djk;
  double dij;

  /* get the payload. */
  enum_prune_path_t *path_data = (enum_prune_path_t*) data;

  /* locally store the level, thread length, and originality array. */
  const unsigned int n = E->G->n_order;
  const unsigned int *dup = E->G->orig;
  const unsigned int j = th->level;

  /* locally store the end position. */
  vector_t *thpos = &th->state[j].pos;

  /* loop over all upstream embedded atoms. */
  for (unsigned int i = j - 1; i < n; i--) {
    /* skip duplicate atoms. */
    if (dup[i]) continue;

    /* compute the distance between the two nodes. */
    dij = vector_dist(&th->state[i].pos, thpos);

    /* loop over all future (un-embedded) atoms. */
    for (unsigned int k = j + 1; k < n; k++) {
      /* skip duplicate atoms. */
      if (dup[k]) continue;

      /* obtain the distance bounds to the future atom. */
      dik = graph_get_edge(E->G, E->G->order[i], E->G->order[k]);
      djk = graph_get_edge(E->G, E->G->order[j], E->G->order[k]);

      /* skip future atoms without defined bounds. */
      if (dik.type == VALUE_TYPE_UNDEFINED ||
          djk.type == VALUE_TYPE_UNDEFINED)
        continue;

      /* prune if the future atom is unreachable. */
      path_data->ntest[i][k - j - 1]++;
      if (dij - dik.u > djk.u) {
        path_data->nprune[i][k - j - 1]++;
        return 1;
      }
    }
  }

  /* do not prune. */
  return 0;
}

/* enum_prune_path_report(): output a report for the shortest path
 * feasbility (i.e. DSP) pruning closure.
 */
void enum_prune_path_report (enum_t *E, unsigned int lev, void *data) {
  /* get the closure payload. */
  enum_prune_path_t *path_data = (enum_prune_path_t*) data;

  /* get the inner atom/residue indices. */
  unsigned int j = lev;
  unsigned int aj = E->G->order[j];
  unsigned int rj = E->P->atoms[aj].res_id;

  /* get the inner atom/residue names. */
  const char *atomj = E->P->atoms[aj].name;
  const char *resj = peptide_get_resname(E->P, rj);

  /* loop over all tested prior neighbors. */
  for (unsigned int i = j - 1; i < E->G->n_order; i--) {
    /* skip duplicate prior atoms. */
    if (E->G->orig[i]) continue;

    /* get the first atom/residue indices. */
    unsigned int ai = E->G->order[i];
    unsigned int ri = E->P->atoms[ai].res_id;

    /* get the first atom/residue names. */
    const char *atomi = E->P->atoms[ai].name;
    const char *resi = peptide_get_resname(E->P, ri);

    /* loop over all tested posterior neighbors. */
    for (unsigned int k = j + 1; k < E->G->n_order; k++) {
      /* skip duplicates and non-pruned posteriors. */
      if (E->G->orig[k] || !path_data->nprune[i][k - j - 1]) continue;

      /* get the last atom/residue indices. */
      unsigned int ak = E->G->order[k];
      unsigned int rk = E->P->atoms[ak].res_id;

      /* get the last atom/residue names. */
      const char *atomk = E->P->atoms[ak].name;
      const char *resk = peptide_get_resname(E->P, rk);

      /* compute the percentage. */
      unsigned int np = path_data->nprune[i][k - j - 1];
      unsigned int nt = path_data->ntest[i][k - j - 1];
      double f = ((double) np) / ((double) nt) * 100.0;

      /* output the statistics. */
      printf("  %3s%-4u %-4s | %3s%-4u %-4s | %3s%-4u %-4s : "
             "%16u/%-16u  %6.2lf%%\n",
             resi, ri + 1, atomi,
             resj, rj + 1, atomj,
             resk, rk + 1, atomk,
             np, nt, f);
    }
  }
}

