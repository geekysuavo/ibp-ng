
/* include the enumerator header. */
#include "enum.h"
#include "enum-node.h"
#include "enum-prune.h"

/* enum_prune_ddf_t: structure for holding information required for
 * direct distance feasibility pruning closures.
 */
typedef struct {
  /* @ntest: number of calls made for each preceeding atom.
   * @nprune: number of prunes performed for each preceeding atom.
   */
  unsigned int *ntest, *nprune;
}
enum_prune_ddf_t;

/* enum_prune_ddf_init(): initialize the ddf enumerator.
 */
int enum_prune_ddf_init (enum_t *E, unsigned int lev) {
  /* declare required variables:
   *  @data: closure data pointer.
   *  @sz: closure memory footprint.
   */
  enum_prune_ddf_t *data;
  unsigned int sz;

  /* allocate the closure payload. */
  sz = sizeof(enum_prune_ddf_t) + 2 * lev * sizeof(unsigned int);
  data = (enum_prune_ddf_t*) malloc(sz);
  if (!data)
    return 0;

  /* initialize the payload array pointers. */
  data->ntest = (unsigned int*)
    (((char*) data) + sizeof(enum_prune_ddf_t));
  data->nprune = (unsigned int*)
    (((char*) data->ntest) + lev * sizeof(unsigned int));

  /* initialize the payload array contents. */
  memset(data->ntest, 0, lev * sizeof(unsigned int));
  memset(data->nprune, 0, lev * sizeof(unsigned int));

  /* register a closure. */
  if (!enum_prune_add_closure(E, lev, enum_prune_ddf, data))
    return 0;

  /* return success. */
  return 1;
}

/* enum_prune_ddf(): determine whether an enumerator tree may be pruned
 * at a given node based on direct distance feasibility (DDF).
 */
int enum_prune_ddf (enum_t *E, enum_node_t *end, void *data) {
  /* declare required variables:
   *  @node: upstream node pointer.
   *  @endpos: pointer to the position of the end.
   *  @bound: distance bound between upstream and end.
   *  @dist: computed distance between upstream and end.
   */
  enum_node_t *node;
  vector_t *endpos;
  value_t bound;
  double dist;

  /* get the payload. */
  enum_prune_ddf_t *ddf_data = (enum_prune_ddf_t*) data;

  /* initialize the local variables. */
  endpos = &end->pos;
  node = end->prev;

  /* loop over all upstream embedded nodes. */
  while (node) {
    /* skip duplicate embeddings. */
    if (E->G->orig[node->lev]) {
      node = node->prev;
      continue;
    }

    /* get the distance bound for the two nodes. */
    bound = graph_get_edge(E->G, E->G->order[node->lev],
                                 E->G->order[end->lev]);

    /* compute the distance between the two nodes. */
    dist = vector_dist(&node->pos, endpos);

    /* prune if the distance is outside of the bound. */
    ddf_data->ntest[node->lev]++;
    if (bound.l - dist > E->ddf_tol || dist - bound.u > E->ddf_tol) {
      ddf_data->nprune[node->lev]++;
      return 1;
    }
  }

  /* do not prune. */
  return 0;
}

/* enum_prune_ddf_report(): output a report for the direct distance
 * feasbility (DDF) pruning closure.
 */
void enum_prune_ddf_report (enum_t *E, unsigned int lev, void *data) {
  /* get the closure payload. */
  enum_prune_ddf_t *ddf_data = (enum_prune_ddf_t*) data;

  /* get the first atom/residue indices. */
  unsigned int ai = E->G->order[lev];
  unsigned int ri = E->P->atoms[ai].res_id;

  /* get the first atom/residue names. */
  const char *atomi = E->P->atoms[ai].name;
  const char resi = peptide_get_reschar(E->P, ri);

  /* loop over all tested predecessors. */
  for (unsigned int j = lev - 1; j < E->G->n_order; j--) {
    /* skip duplicates and non-pruned predecessors. */
    if (E->G->orig[j] || !ddf_data->nprune[j]) continue;

    /* get the second atom/residue indices. */
    unsigned int aj = E->G->order[j];
    unsigned int rj = E->P->atoms[aj].res_id;

    /* get the second atom/residue names. */
    const char *atomj = E->P->atoms[aj].name;
    const char resj = peptide_get_reschar(E->P, rj);

    /* compute the percentage. */
    double f = ((double) ddf_data->nprune[j]) /
               ((double) ddf_data->ntest[j]) * 100.0;

    /* output the statistics. */
    printf("  %c%-4u %-4s | %c%-4u %-4s : "
           "%16u/%-16u  %3.0lf%%\n",
           resi, ri + 1, atomi,
           resj, rj + 1, atomj,
           ddf_data->nprune[j], ddf_data->ntest[j], f);
  }
}

