
/* include the enumerator headers. */
#include "enum.h"
#include "enum-node.h"
#include "enum-write.h"

/* enum_node_init(): initialize an enumeration tree node.
 *
 * arguments:
 *  @node: pointer to the node to initialize.
 *  @prev: pointer to the parent node, or NULL for a root node.
 */
void enum_node_init (enum_node_t *node, enum_node_t *prev) {
  /* initialize the node parent. */
  node->prev = prev;

  /* initialize the node children. */
  node->next = NULL;
  node->nb = 0;

  /* initialize the node position and energy. */
  memset(&node->pos, 0, sizeof(vector_t));
  node->energy = 0.0;

  /* initialize the node feasibility flag. */
  node->feas = 0;

  /* determine the node level. */
  node->lev = (node->prev ? node->prev->lev + 1 : 0);
}

/* enum_node_free(): free all allocated memory associated with an
 * enumeration tree node.
 *
 * arguments:
 *  @E: pointer to the enumerator structure to access.
 *  @node: pointer to the tree node structure to free.
 */
void enum_node_free (enum_t *E, enum_node_t *node) {
  /* declare required variables:
   *  @i: general-purpose loop counter.
   */
  unsigned int i;

  /* return if the structure pointer is null. */
  if (!node) return;

  /* check if the node has allocated children. */
  if (node->next) {
    /* loop over the child nodes. */
    for (i = 0; i < node->nb; i++)
      enum_node_free(E, node->next[i]);

    /* free the child node array. */
    free(node->next);
  }
}

/* enum_node_set_branches(): set the number of child branches in an
 * enumerator tree node, and reallocate the array of child nodes.
 *
 * arguments:
 *  @node: pointer to the enumerator tree node to modify.
 *  @nb: number of branches to set for the node.
 *
 * returns:
 *  integer indicating whether (1) or not (0) the operation succeeded.
 */
int enum_node_set_branches (enum_node_t *node, unsigned int nb) {
  /* declare required variables:
   *  @i: node array index.
   */
  unsigned int i;

  /* fail if the node already has branches. */
  if (node->next)
    throw("node already contains branches");

  /* store the new number of branches. */
  node->nb = nb;

  /* reallocate the child node array. */
  unsigned long sz = node->nb * (sizeof(enum_node_t) + sizeof(enum_node_t*));
  node->next = malloc(sz);

  /* check if reallocated failed. */
  if (!node->next)
    throw("unable to reallocate node branches");

  /* allocate the child nodes. */
  for (i = 0; i < node->nb; i++) {
    /* initialize the current child node. */
    node->next[i] = ((enum_node_t*) (node->next + nb)) + i;
    enum_node_init(node->next[i], node);
  }

  /* return success. */
  return 1;
}

/* enum_node_feasible(): determine whether an enumeration tree node
 * is feasible or not.
 *
 * arguments:
 *  @E: pointer to the enumerator data structure to utilize.
 *  @end: pointer to the tree node to feasibility-check.
 *
 * returns:
 *  unsigned integer indicating whether (1) or not (0) the node is
 *  feasible.
 */
unsigned int enum_node_feasible (enum_t *E, enum_node_t *end) {
  /* declare required variables:
   *  @func: pruning function pointer array.
   *  @data: pruning data payload array.
   *  @i: pruning function array index.
   */
  enum_prune_test_fn *func;
  unsigned int i, n;
  void **data;

  /* get the arrays for the current position in the order. */
  data = E->prune_data[end->lev];
  func = E->prune[end->lev];
  n = E->prune_sz[end->lev];

  /* loop over the enumerator array of pruning function pointers. */
  for (i = 0; i < n; i++) {
    /* check if the current node is pruned by the function. */
    if ((func[i])(E, end, data[i]))
      return 0;
  }

  /* return true. */
  return 1;
}

/* enum_node_compute(): core interval branch and prune (iBP) recursive
 * computation routine. enumerates all solutions to an iDMDGP problem
 * instance.
 *
 * algorithm:
 *   Carlile Lavor, Leo Liberti, Antonio Mucherino, "The interval
 *   Branch-and-Prune algorithm for the discretizable molecular
 *   distance geometry problem with exact distances",
 *   J. Glob. Optim., 2013, 56: 855--871.
 *
 * computation:
 *   Carlile Lavor, Rafael Alves, Weber Figueiredo, Antonio Petraglia
 *   and Nelson Maculan, "Clifford Algebra and the Discretizable
 *   Molecular Distance Geometry Problem", Adv. Appl. Clifford
 *   Algebras, 2015, 25: 925--942.
 *
 * notation:
 *   x0: point x_{i-3}
 *   x1: point x_{i-2}
 *   x2: point x_{i-1}
 *   x3: point x_i
 *
 *   d01: dist (x_{i-3}, x_{i-2})    [exact]
 *   d02: dist (x_{i-3}, x_{i-1})    [exact]
 *   d03: dist (x_{i-3}, x_i)        [exact|interval]    [from graph]
 *   d12: dist (x_{i-2}, x_{i-1})    [exact]
 *   d13: dist (x_{i-2}, x_i)        [exact]             [from graph]
 *   d23: dist (x_{i-1}, x_i)        [exact]             [from graph]
 *
 *   theta: angle (x_{i-2}, x_{i-1}, x_i)
 *   omega: dihedral (x_{i-3}, x_{i-2}, x_{i-1}, x_i)
 *
 * arguments:
 *  @E: pointer to the enumerator structure to utilize.
 *  @node: pointer to the current computation node.
 *
 * returns:
 *  integer indicating whether (1) or not (0) computation succeeded.
 */
int enum_node_compute (enum_t *E, enum_node_t *node) {
  /* FIXME: implement enum_node_compute() */

  /* return success. */
  return 1;
}

