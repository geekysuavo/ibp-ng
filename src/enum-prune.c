
/* include the enumerator header. */
#include "enum.h"

/* all functions defined in the "enum-prune-*.[ch]" source files must
 * follow the function pointer specification outlined for enumerator
 * tree pruning.
 *
 * pruning init methods follow the enum_prune_init_fn specification,
 * and modify @prune, @prune_sz and @prune_data in an enumerator to
 * make it use that pruning method at specific points in the enumeration.
 *
 * pruning test methods follow the enum_prune_test_fn specification,
 * and return whether or not a tree may be pruned.
 *
 * for more details, consult:
 *  - enum_node_is_feasible() in enum-node.c
 *  - enum_prune_init_fn in enum.h
 *  - enum_prune_test_fn in enum.h
 */

/* enum_prune_add_closure(): register a single pruning method (with
 * accompanying data) with an enumerator structure.
 *
 * arguments:
 *  @E: pointer to the enumerator structure to modify.
 *  @lev: graph order index at which to add the closure.
 *  @func: pruning test function pointer to register.
 *  @data: pruning data payload to register.
 *
 * returns:
 *  integer indicating whether (1) or not (0) the operation succeeded.
 */
int enum_prune_add_closure (enum_t *E, unsigned int lev,
                            enum_prune_test_fn func,
                            void *data) {
  /* increment the size of the array. */
  E->prune_sz[lev]++;

  /* reallocate the function array. */
  E->prune[lev] = (enum_prune_test_fn*)
    realloc(E->prune[lev], E->prune_sz[lev] * sizeof(enum_prune_test_fn));

  /* reallocate the data array. */
  E->prune_data[lev] = (void**)
    realloc(E->prune_data[lev], E->prune_sz[lev] * sizeof(void*));

  /* check if either reallocation failed. */
  if (!E->prune[lev] || !E->prune_data[lev])
    throw("unable to reallocate pruning method arrays");

  /* store the closure. */
  E->prune[lev][E->prune_sz[lev] - 1] = func;
  E->prune_data[lev][E->prune_sz[lev] - 1] = data;

  /* return success. */
  return 1;
}

