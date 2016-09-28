
/* include the enumerator header. */
#include "enum.h"
#include "enum-thread.h"
#include "enum-prune.h"

/* define constants to represent pruning closure @type.
 */
#define ENERGY_BOND     0
#define ENERGY_ANGLE    1
#define ENERGY_DIHEDRAL 2
#define ENERGY_DISTANCE 3

/* enum_prune_energy_t: structure for holding information
 * required for energetic pruning closures.
 */
typedef struct {
  /* @type: type of energy term to be feasibility checked.
   */
  unsigned int type;

  /* @ntest, @nprune: test and prune counts of the current closure.
   * @mu, @kappa: mean and precision force field parameters.
   * @n: backward step-counts for the prior atoms.
   */
  unsigned int ntest, nprune;
  unsigned int n[4];
  double mu, kappa;
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

