
/* ensure once-only inclusion. */
#pragma once

/* include the c math library header. */
#include <math.h>

/* include the peptide header. */
#include "peptide.h"

/* include the traceback and string function headers. */
#include "trace.h"
#include "str.h"

/* topol_mode_t: enumeration of all types of patch statements that may
 * be set for a topology structure.
 */
typedef enum {
  TOPOL_MODE_ADD = 0,
  TOPOL_MODE_MODIFY,
  TOPOL_MODE_DELETE
}
topol_mode_t;

/* topol_atom_t: structure for holding atomic information.
 */
typedef struct {
  /* @name: name of the atom within a given residue.
   * @type: global type string of the atom.
   */
  char *name, *type;

  /* @mass: mass of the atom.
   * @charge: partial charge on the atom.
   */
  double mass, charge;

  /* @mode: atom topology mode.
   * @off: residue offset.
   */
  topol_mode_t mode;
  int off;
}
topol_atom_t;

/* topol_bond_t: structure for holding molecular bonding connectivities
 * between pairs of atoms.
 */
typedef struct {
  /* @atoms: atoms involved in the connection.
   * @off: residue offsets of the connection.
   */
  char *atoms[2];
  int off[2];

  /* @mode: bond topology mode.
   */
  topol_mode_t mode;
}
topol_bond_t;

/* topol_angle_t: structure for holding angle connectivities between sets
 * of three atoms.
 */
typedef struct {
  /* @atoms: atoms involved in the connection.
   * @off: residue offsets of the connection.
   */
  char *atoms[3];
  int off[3];

  /* @mode: angle topology mode.
   */
  topol_mode_t mode;
}
topol_angle_t;

/* topol_dihedral_t: structure for holding dihedral connectivities between
 * sets of four atoms.
 */
typedef struct {
  /* @atoms: atoms involved in the connection.
   * @off: residue offsets of the connection.
   */
  char *atoms[4];
  int off[4];

  /* @mode: dihedral topology mode.
   */
  topol_mode_t mode;
}
topol_dihedral_t;

/* topol_residue_t: structure for holding the topology of a single residue.
 */
typedef struct {
  /* @name: residue name string.
   */
  char *name;

  /* @atoms: array of atoms.
   * @n_atoms: number of atoms.
   */
  topol_atom_t *atoms;
  unsigned int n_atoms;

  /* @bonds: array of two-atom connectivities.
   * @n_bonds: number of two-atom connectivities.
   */
  topol_bond_t *bonds;
  unsigned int n_bonds;

  /* @angles: array of three-atom connectivities.
   * @n_angles: number of three-atom connectivities.
   */
  topol_angle_t *angles;
  unsigned int n_angles;

  /* @torsions: array of torsional dihedral connectivities.
   * @n_torsions: number of torsional dihedral connectivities.
   */
  topol_dihedral_t *torsions;
  unsigned int n_torsions;

  /* @impropers: array of improper dihedral connectivities.
   * @n_impropers: number of improper dihedral connectivities.
   */
  topol_dihedral_t *impropers;
  unsigned int n_impropers;

  /* @patch: whether or not the residue is a patch residue.
   */
  unsigned int patch;
}
topol_residue_t;

/* topol_mass_t: structure for holding default masses of all atoms.
 */
typedef struct {
  /* @type: atom type string.
   * @mass: default mass.
   */
  char *type;
  double mass;
}
topol_mass_t;

/* topol_t: structure for holding peptide molecular topology information.
 */
typedef struct {
  /* @res: array of residue topologies.
   * @n_res: number of residue topologies.
   */
  topol_residue_t *res;
  unsigned int n_res;

  /* @mass: array of atom masses.
   * @n_mass: number of atom masses.
   */
  topol_mass_t *mass;
  unsigned int n_mass;

  /* @auto_angles: whether or not to autogenerate angles.
   * @auto_dihedrals: whether or not to autogenerate dihedrals.
   */
  unsigned int auto_angles, auto_dihedrals;
}
topol_t;

/* function declarations (topol-alloc.c): */

topol_t *topol_new (void);

topol_t *topol_new_from_file (const char *fname);

void topol_free (topol_t *top);

/* function declarations (topol.c): */

int topol_apply_all (topol_t *top, peptide_t *P);

