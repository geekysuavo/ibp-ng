
/* ensure once-only inclusion. */
#pragma once

/* include the traceback, graph and residue headers. */
#include "trace.h"
#include "graph.h"

/* include the repetition ordering header. */
#include "reorder.h"

/* peptide_atom_t: structure for holding per-atom information.
 */
typedef struct {
  /* @res_id: index of the residue to which the atom belongs.
   */
  unsigned int res_id;

  /* @name: name of atom, specific to the residue.
   * @type: type of atom, specific to the residue.
   */
  char *name, *type;

  /* @mass: mass of the atom.
   * @charge: charge on the atom.
   * @radius: atomic radius.
   */
  double mass, charge, radius;
}
peptide_atom_t;

/* peptide_bond_t: structure for holding per-bond information.
 */
typedef struct {
  /* @atom_id: indices of the bonded atoms.
   * @len: distance between the bonded atoms.
   */
  unsigned int atom_id[2];
  value_t len;

  /* @is_virtual: whether the bond is a non-physical bond
   * (i.e. from an interval restraint) or not.
   */
  unsigned int is_virtual;

  /* @mu: mode of the normal/log-normal probability.
   * @kappa: precision of the normal/log-normal probability.
   */
  double mu, kappa;
}
peptide_bond_t;

/* peptide_angle_t: structure for holding per-angle information.
 */
typedef struct {
  /* @atom_id: indices of the involved atoms.
   * @ang: angle formed by the atoms.
   */
  unsigned int atom_id[3];
  value_t ang;

  /* @mu: mode of the von mises probability.
   * @kappa: precision of the von mises probability.
   */
  double mu, kappa;
}
peptide_angle_t;

/* peptide_dihed_t: structure for holding per-dihedral information.
 */
typedef struct {
  /* @atom_id: indices of the involved atoms.
   * @ang: angle formed by the planes.
   */
  unsigned int atom_id[4];
  value_t ang;

  /* @mu: mode of the von mises probability.
   * @kappa: precision of the von mises probability.
   */
  double mu, kappa;
}
peptide_dihed_t;

/* peptide_t: structure for holding a single peptide chain.
 *
 * the information assimilated into this data structure is used to
 * construct an iDMDGP instance graph.
 */
typedef struct {
  /* @n_res: number of residues in the sequence.
   * @res: array of residue name strings in the sequence.
   */
  unsigned int n_res;
  char **res;

  /* @sc: array of residue indices having explicit sidechains.
   * @n_sc: number of residues having explicit sidechains.
   */
  unsigned int *sc, n_sc;

  /* @atoms: array of atoms in the peptide.
   * @n_atoms: number of atoms in the peptide.
   */
  peptide_atom_t *atoms;
  unsigned int n_atoms;

  /* @bonds: array of bonds in the peptide.
   * @n_bonds: number of bonds in the peptide.
   */
  peptide_bond_t *bonds;
  unsigned int n_bonds;

  /* @angles: array of known angles in the peptide.
   * @n_angles: number of known angles in the peptide.
   */
  peptide_angle_t *angles;
  unsigned int n_angles;

  /* @torsions: array of torsions in the peptide.
   * @n_torsions: number of torsions in the peptide.
   */
  peptide_dihed_t *torsions;
  unsigned int n_torsions;

  /* @impropers: array of impropers in the peptide.
   * @n_impropers: number of impropers in the peptide.
   */
  peptide_dihed_t *impropers;
  unsigned int n_impropers;
}
peptide_t;

/* function declarations (peptide-alloc.c): */

peptide_t *peptide_new (void);

peptide_t *peptide_new_from_file (const char *fname, const char *opts);

void peptide_free (peptide_t *P);

/* function declarations (peptide-residues.c): */

int peptide_add_residue (peptide_t *P, const char *res);

int peptide_add_sidechain (peptide_t *P, unsigned int res);

int peptide_has_sidechain (peptide_t *P, unsigned int res);

const char *peptide_get_resname (peptide_t *P, unsigned int idx);

const char *peptide_get_restype (peptide_t *P, unsigned int idx);

/* function declarations (peptide-field.c): */

int peptide_field (peptide_t *P, double tol);

/* function declarations (peptide-graph.c): */

graph_t *peptide_graph (peptide_t *P, reorder_t *ord,
                        unsigned int refine,
                        unsigned int complete);

