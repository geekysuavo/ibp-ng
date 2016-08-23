
/* include the molecular parameters header. */
#include "param.h"

/* param_add_bond(): add a bond entry into a parameter structure.
 *
 * arguments:
 *  @par: pointer to the parameter structure to modify.
 *  @a: string name of the first atom in the entry.
 *  @b: string name of the second atom in the entry.
 *  @length: numerical value of the entry.
 *
 * returns:
 *  integer indicating whether (1) or not (0) the operation succeeded.
 */
int param_add_bond (param_t *par,
                    const char *a,
                    const char *b,
                    value_t length) {
  /* declare required variables:
   *  @i: index of the new array entry.
   */
  unsigned int i;

  /* increment the array size. */
  i = par->n_bonds;
  par->n_bonds++;

  /* reallocate the array. */
  par->bonds = (param_bond_t*)
    realloc(par->bonds, par->n_bonds * sizeof(param_bond_t));

  /* check if reallocation failed. */
  if (!par->bonds)
    throw("unable to reallocate bond array");

  /* store the strings of the new array entry. */
  par->bonds[i].a = (char*) a;
  par->bonds[i].b = (char*) b;

  /* store the numerical value of the new array entry. */
  par->bonds[i].v = length;

  /* return success. */
  return 1;
}

/* param_add_angle(): add an angle entry into a parameter structure.
 *
 * arguments:
 *  @par: pointer to the parameter structure to modify.
 *  @a: string name of the first atom in the entry.
 *  @b: string name of the second atom in the entry.
 *  @c: string name of the third atom in the entry.
 *  @angle: numerical value of the entry.
 *
 * returns:
 *  integer indicating whether (1) or not (0) the operation succeeded.
 */
int param_add_angle (param_t *par,
                     const char *a,
                     const char *b,
                     const char *c,
                     value_t angle) {
  /* declare required variables:
   *  @i: index of the new array entry.
   */
  unsigned int i;

  /* increment the array size. */
  i = par->n_angles;
  par->n_angles++;

  /* reallocate the array. */
  par->angles = (param_angle_t*)
    realloc(par->angles, par->n_angles * sizeof(param_angle_t));

  /* check if reallocation failed. */
  if (!par->angles)
    throw("unable to reallocate angle array");

  /* store the strings of the new array entry. */
  par->angles[i].a = (char*) a;
  par->angles[i].b = (char*) b;
  par->angles[i].c = (char*) c;

  /* store the numerical value of the new array entry. */
  par->angles[i].v = angle;

  /* return success. */
  return 1;
}

/* param_add_torsion(): add a torsion angle entry into a parameter structure.
 *
 * arguments:
 *  @par: pointer to the parameter structure to modify.
 *  @a: string name of the first atom in the entry.
 *  @b: string name of the second atom in the entry.
 *  @c: string name of the third atom in the entry.
 *  @d: string name of the fourth atom in the entry.
 *  @angle: numerical value of the entry.
 *
 * returns:
 *  integer indicating whether (1) or not (0) the operation succeeded.
 */
int param_add_torsion (param_t *par,
                       const char *a,
                       const char *b,
                       const char *c,
                       const char *d,
                       value_t angle) {
  /* declare required variables:
   *  @i: index of the new array entry.
   */
  unsigned int i;

  /* increment the array size. */
  i = par->n_torsions;
  par->n_torsions++;

  /* reallocate the array. */
  par->torsions = (param_dihedral_t*)
    realloc(par->torsions, par->n_torsions * sizeof(param_dihedral_t));

  /* check if reallocation failed. */
  if (!par->torsions)
    throw("unable to reallocate torsion array");

  /* store the strings of the new array entry. */
  par->torsions[i].a = (char*) a;
  par->torsions[i].b = (char*) b;
  par->torsions[i].c = (char*) c;
  par->torsions[i].d = (char*) d;

  /* store the numerical value of the new array entry. */
  par->torsions[i].v = angle;

  /* return success. */
  return 1;
}

/* param_add_improper(): add an improper dihedral angle entry into
 * a parameter structure.
 *
 * arguments:
 *  @par: pointer to the parameter structure to modify.
 *  @a: string name of the first atom in the entry.
 *  @b: string name of the second atom in the entry.
 *  @c: string name of the third atom in the entry.
 *  @d: string name of the fourth atom in the entry.
 *  @angle: numerical value of the entry.
 *
 * returns:
 *  integer indicating whether (1) or not (0) the operation succeeded.
 */
int param_add_improper (param_t *par,
                        const char *a,
                        const char *b,
                        const char *c,
                        const char *d,
                        value_t angle) {
  /* declare required variables:
   *  @i: index of the new array entry.
   */
  unsigned int i;

  /* increment the array size. */
  i = par->n_impropers;
  par->n_impropers++;

  /* reallocate the array. */
  par->impropers = (param_dihedral_t*)
    realloc(par->impropers, par->n_impropers * sizeof(param_dihedral_t));

  /* check if reallocation failed. */
  if (!par->impropers)
    throw("unable to reallocate improper array");

  /* store the strings of the new array entry. */
  par->impropers[i].a = (char*) a;
  par->impropers[i].b = (char*) b;
  par->impropers[i].c = (char*) c;
  par->impropers[i].d = (char*) d;

  /* store the numerical value of the new array entry. */
  par->impropers[i].v = angle;

  /* return success. */
  return 1;
}

/* param_add_radius(): add a nonbonding van der waals radius entry into
 * a parameter structure.
 *
 * arguments:
 *  @par: pointer to the parameter structure to modify.
 *  @type: string type of the atom in the entry.
 *  @sigma: scaled radius of the entry.
 *
 * returns:
 *  integer indicating whether (1) or not (0) the operation succeeded.
 */
int param_add_radius (param_t *par, const char *type, double sigma) {
  /* declare required variables:
   *  @i: index of the new array entry.
   */
  unsigned int i;

  /* increment the array size. */
  i = par->n_radii;
  par->n_radii++;

  /* reallocate the array. */
  par->radii = (param_radius_t*)
    realloc(par->radii, par->n_radii * sizeof(param_radius_t));

  /* check if reallocation failed. */
  if (!par->radii)
    throw("unable to reallocate radius array");

  /* store the atom type string of the new array entry. */
  par->radii[i].a = (char*) type;

  /* store the radius of the array entry:
   *  equation: R_vdw = (sigma/2) * rt^6(2)
   */
  par->radii[i].v = value_scalar(par->vdw_scale * sigma *
                                 pow(2.0, -5.0 / 6.0));

  /* return success. */
  return 1;
}

