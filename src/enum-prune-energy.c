
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
#define ENERGY_CONTACT  4

/* enum_prune_energy_t: structure for holding information
 * required for energetic pruning closures.
 */
typedef struct enum_prune_energy enum_prune_energy_t;
struct enum_prune_energy {
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

  /* @next: pointer to the next energy term structure.
   */
  enum_prune_energy_t *next;
};

/* enum_prune_energy_init(): initialize the energy enumerator.
 */
int enum_prune_energy_init (enum_t *E, unsigned int lev) {
  /* declare required variables:
   *  @i: peptide bond/angle/torsion/improper array index.
   *  @k: atoms array index.
   *  @n_data: number of linked closure data structures.
   *  @id: atom index that has just been embedded.
   *  @ids: atom indices in each energetic term.
   *  @levs: reorder indices of each atom.
   *  @data: closure data pointer.
   */
  unsigned int i, k, n, n_data, id, *ids, levs[4];
  enum_prune_energy_t *data;

  /* define a contact force constant scale factor. */
  const double Z = 3.84147997;

  /* get the current atom index. */
  id = E->G->order[lev];

  /* initialize the closure data. */
  data = NULL;
  n_data = 0;

  /* loop over the bonds. */
  for (i = 0, n = 2; i < E->P->n_bonds; i++) {
    /* get the atom indices. */
    ids = E->P->bonds[i].atom_id;

    /* skip if the last-embedded atom is not in the bond. */
    if (id != ids[0] && id != ids[1])
      continue;

    /* get the graph level of each atom. */
    levs[0] = levs[1] = levs[2] = levs[3] = lev;
    for (k = 0; k < n; k++)
      levs[k] = enum_prune_get_level(E->G->order, lev, ids[k]);

    /* skip if all other atoms in the order have not been embedded. */
    if (levs[0] > lev ||
        levs[1] > lev)
      continue;

    /* reallocate the closure payload. */
    n_data++;
    data = (enum_prune_energy_t*)
      realloc(data, n_data * sizeof(enum_prune_energy_t));

    /* check for allocation failures. */
    if (!data)
      return 0;

    /* store the term type. */
    data[n_data - 1].type = (E->P->bonds[i].is_virtual
                              ? ENERGY_DISTANCE
                              : ENERGY_BOND);

    /* initialize the counters. */
    data[n_data - 1].ntest = data[n_data - 1].nprune = 0;

    /* store the offsets. */
    data[n_data - 1].n[0] = 0;
    data[n_data - 1].n[1] = 0;
    data[n_data - 1].n[2] = 0;
    data[n_data - 1].n[3] = 0;
    for (k = 0; k < n; k++)
      data[n_data - 1].n[k] = lev - levs[k];

    /* store the term parameters. */
    data[n_data - 1].mu = E->P->bonds[i].mu;
    data[n_data - 1].kappa = E->P->bonds[i].kappa;

    /* initialize the next-payload pointer. */
    data[n_data - 1].next = NULL;
  }

  /* loop over the angles. */
  for (i = 0, n = 3; i < E->P->n_angles; i++) {
    /* get the atom indices. */
    ids = E->P->angles[i].atom_id;

    /* skip if the last-embedded atom is not in the angle. */
    if (id != ids[0] && id != ids[1] && id != ids[2])
      continue;

    /* get the graph level of each atom. */
    levs[0] = levs[1] = levs[2] = levs[3] = lev;
    for (k = 0; k < n; k++)
      levs[k] = enum_prune_get_level(E->G->order, lev, ids[k]);

    /* skip if all other atoms in the order have not been embedded. */
    if (levs[0] > lev ||
        levs[1] > lev ||
        levs[2] > lev)
      continue;

    /* reallocate the closure payload. */
    n_data++;
    data = (enum_prune_energy_t*)
      realloc(data, n_data * sizeof(enum_prune_energy_t));

    /* check for allocation failures. */
    if (!data)
      return 0;

    /* store the term type and initialize the counters. */
    data[n_data - 1].type = ENERGY_ANGLE;
    data[n_data - 1].ntest = data[n_data - 1].nprune = 0;

    /* store the offsets. */
    data[n_data - 1].n[0] = 0;
    data[n_data - 1].n[1] = 0;
    data[n_data - 1].n[2] = 0;
    data[n_data - 1].n[3] = 0;
    for (k = 0; k < n; k++)
      data[n_data - 1].n[k] = lev - levs[k];

    /* store the term parameters. */
    data[n_data - 1].mu = E->P->angles[i].mu;
    data[n_data - 1].kappa = E->P->angles[i].kappa;

    /* initialize the next-payload pointer. */
    data[n_data - 1].next = NULL;
  }

  /* loop over the torsions. */
  for (i = 0, n = 4; i < E->P->n_torsions; i++) {
    /* get the atom indices. */
    ids = E->P->torsions[i].atom_id;

    /* skip if the last-embedded atom is not in the torsion. */
    if (id != ids[0] && id != ids[1] && id != ids[2] && id != ids[3])
      continue;

    /* get the graph level of each atom. */
    levs[0] = levs[1] = levs[2] = levs[3] = lev;
    for (k = 0; k < n; k++)
      levs[k] = enum_prune_get_level(E->G->order, lev, ids[k]);

    /* skip if all other atoms in the order have not been embedded. */
    if (levs[0] > lev ||
        levs[1] > lev ||
        levs[2] > lev ||
        levs[3] > lev)
      continue;

    /* reallocate the closure payload. */
    n_data++;
    data = (enum_prune_energy_t*)
      realloc(data, n_data * sizeof(enum_prune_energy_t));

    /* check for allocation failures. */
    if (!data)
      return 0;

    /* store the term type and initialize the counters. */
    data[n_data - 1].type = ENERGY_DIHEDRAL;
    data[n_data - 1].ntest = data[n_data - 1].nprune = 0;

    /* store the offsets. */
    data[n_data - 1].n[0] = 0;
    data[n_data - 1].n[1] = 0;
    data[n_data - 1].n[2] = 0;
    data[n_data - 1].n[3] = 0;
    for (k = 0; k < n; k++)
      data[n_data - 1].n[k] = lev - levs[k];

    /* store the term parameters. */
    data[n_data - 1].mu = E->P->torsions[i].mu;
    data[n_data - 1].kappa = E->P->torsions[i].kappa;

    /* initialize the next-payload pointer. */
    data[n_data - 1].next = NULL;
  }

  /* loop over the impropers. */
  for (i = 0, n = 4; i < E->P->n_impropers; i++) {
    /* get the atom indices. */
    ids = E->P->impropers[i].atom_id;

    /* skip if the last-embedded atom is not in the improper. */
    if (id != ids[0] && id != ids[1] && id != ids[2] && id != ids[3])
      continue;

    /* get the graph level of each atom. */
    levs[0] = levs[1] = levs[2] = levs[3] = lev;
    for (k = 0; k < n; k++)
      levs[k] = enum_prune_get_level(E->G->order, lev, ids[k]);

    /* skip if all other atoms in the order have not been embedded. */
    if (levs[0] > lev ||
        levs[1] > lev ||
        levs[2] > lev ||
        levs[3] > lev)
      continue;

    /* reallocate the closure payload. */
    n_data++;
    data = (enum_prune_energy_t*)
      realloc(data, n_data * sizeof(enum_prune_energy_t));

    /* check for allocation failures. */
    if (!data)
      return 0;

    /* store the term type and initialize the counters. */
    data[n_data - 1].type = ENERGY_DIHEDRAL;
    data[n_data - 1].ntest = data[n_data - 1].nprune = 0;

    /* store the offsets. */
    data[n_data - 1].n[0] = 0;
    data[n_data - 1].n[1] = 0;
    data[n_data - 1].n[2] = 0;
    data[n_data - 1].n[3] = 0;
    for (k = 0; k < n; k++)
      data[n_data - 1].n[k] = lev - levs[k];

    /* store the term parameters. */
    data[n_data - 1].mu = E->P->impropers[i].mu;
    data[n_data - 1].kappa = E->P->impropers[i].kappa;

    /* initialize the next-payload pointer. */
    data[n_data - 1].next = NULL;
  }

  /* loop over van der waals contacts. */
  for (i = 0, n = 4; i < lev; i++) {
    /* skip duplicate atoms. */
    if (E->G->orig[i]) continue;

    /* get the other atom index. */
    unsigned int jd = E->G->order[i];

    /* reallocate the closure payload. */
    n_data++;
    data = (enum_prune_energy_t*)
      realloc(data, n_data * sizeof(enum_prune_energy_t));

    /* check for allocation failures. */
    if (!data)
      return 0;

    /* store the term type. */
    data[n_data - 1].type = ENERGY_CONTACT;

    /* initialize the counters. */
    data[n_data - 1].ntest = data[n_data - 1].nprune = 0;

    /* store the offsets. */
    data[n_data - 1].n[0] = 0;
    data[n_data - 1].n[1] = lev - i;
    data[n_data - 1].n[2] = 0;
    data[n_data - 1].n[3] = 0;

    /* store the term parameters. */
    data[n_data - 1].kappa = Z * pow(E->ddf_tol, -2.0);
    data[n_data - 1].mu =
      E->P->atoms[id].radius + E->P->atoms[jd].radius;

    /* initialize the next-payload pointer. */
    data[n_data - 1].next = NULL;
  }

  /* return if no closures were created. */
  if (!n_data)
    return 1;

  /* link up the chain of energy structures. */
  for (i = 0; i < n_data - 1; i++)
    data[i].next = data + (i + 1);

  /* register the closure with the enumerator. */
  if (!enum_prune_add_closure(E, lev, enum_prune_energy, data))
    return 0;

  /* return success. */
  return 1;
}

/* enum_prune_energy(): determine whether an enumerator tree may
 * be pruned at a given node based on energetic feasibility.
 */
int enum_prune_energy (enum_t *E, enum_thread_t *th, void *data) {
  /* declare required variables:
   *  @energy_data: payload for energetic pruning terms.
   *  @Eterm: energy contribution of the current term.
   *  @Enew: total energy contribution of the embedded atom.
   *  @mu: mean parameter of the current term.
   *  @kappa: precision parameter of the current term.
   *  @obs: observed value of the current term.
   *  @x: array of atoms involved in the current term.
   */
  enum_prune_energy_t *energy_data;
  double Eterm, Enew, mu, kappa, obs;
  unsigned int k;
  vector_t x[4];

  /* get the payload and initialize the energy contribution. */
  energy_data = (enum_prune_energy_t*) data;
  Enew = 0.0;

  /* loop over each term in the new energy contribution. */
  while (energy_data) {
    /* extract the term parameters. */
    Eterm = 0.0;
    mu = energy_data->mu;
    kappa = energy_data->kappa;

    /* extract all atom positions. */
    for (k = 0; k < 4; k++)
      x[k] = th->state[th->level - energy_data->n[k]].pos;

    /* compute the current term. */
    switch (energy_data->type) {
      /* bonded distance. */
      case ENERGY_BOND:
        obs = vector_dist(&x[0], &x[1]);
        Eterm = 0.5 * kappa * pow(obs - mu, 2.0);
        break;

      /* two-bond angle. */
      case ENERGY_ANGLE:
        obs = vector_angle(&x[0], &x[1], &x[2]);
        Eterm = 0.5 * kappa * pow(obs - mu, 2.0);
        break;

      /* dihedral angle. */
      case ENERGY_DIHEDRAL:
        obs = vector_dihedral(&x[0], &x[1], &x[2], &x[3]);
        Eterm = 0.5 * kappa * pow(obs - mu, 2.0);
        break;

      /* non-bonded distance. */
      case ENERGY_DISTANCE:
        obs = vector_dist(&x[0], &x[1]);
        Eterm = 0.5 * kappa * pow(log(obs / mu), 2.0);
        break;

      /* close contacts (vdw repulsion). */
      case ENERGY_CONTACT:
        obs = vector_dist(&x[0], &x[1]);
        Eterm = 0.5 * kappa * pow(mu / obs, 6.0);
        break;

      /* otherwise, do nothing. */
      default: break;
    }

    /* sum the energy term into the new contribution. */
    Enew += Eterm;

    /* move to the next energy term. */
    energy_data = energy_data->next;
  }

  /* update the energy at the current node. */
  th->state[th->level].energy = th->state[th->level - 1].energy + Enew;

  /* check if the node should be pruned. */
  if (th->state[th->level].energy > E->energy_tol)
    return 1;

  /* do not prune. */
  return 0;
}

/* enum_prune_energy_report(): output a report for the energetic feasibility
 * pruning closure.
 */
void enum_prune_energy_report (enum_t *E, unsigned int lev, void *data) {
  /* FIXME: implement enum_prune_energy_report() */
}

