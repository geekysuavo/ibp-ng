
/* include the enumerator header. */
#include "enum.h"
#include "enum-thread.h"
#include "enum-prune.h"

/* enum_prune_energy_t: structure for holding information required for
 * energetic feasibility pruning closures.
 */
typedef struct {
  /* FIXME: fill out the enum_prune_energy_t struct. */
}
enum_prune_energy_t;

/* enum_prune_energy_init(): initialize the energy enumerator.
 */
int enum_prune_energy_init (enum_t *E, unsigned int lev) {
  /* FIXME: implement enum_prune_energy_init() */

  /* return success. */
  return 1;
}

/* enum_prune_energy(): determine whether an enumerator tree may be pruned
 * at a given node based on energetic feasibility.
 */
int enum_prune_energy (enum_t *E, enum_thread_t *th, void *data) {
  /* FIXME: implement enum_prune_energy() */

  /* do not prune. */
  return 0;
}

/* enum_prune_ddf_report(): output a report for the energetic feasibility
 * pruning closure.
 */
void enum_prune_energy_report (enum_t *E, unsigned int lev, void *data) {
  /* FIXME: implement enum_prune_energy_report() */
}

