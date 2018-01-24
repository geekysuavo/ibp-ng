
/* ensure once-only inclusion. */
#pragma once

/* include the c math library header. */
#include <math.h>

/* include the peptide header. */
#include "peptide.h"

/* param_bond_t: structure for holding bond distances.
 */
typedef struct {
  /* @a, @b: atoms involved in the bond.
   * @v: length of the bond.
   */
  char *a, *b;
  value_t v;
}
param_bond_t;

/* param_angle_t: structure for holding two-bond angles.
 */
typedef struct {
  /* @a, @b, @c: atoms involved in the bond.
   * @v: angle of the bond.
   */
  char *a, *b, *c;
  value_t v;
}
param_angle_t;

/* param_dihedral_t: structure for holding three-bond angles.
 */
typedef struct {
  /* @a, @b, @c, @d: atoms involved in the dihedral unit.
   * @v: angle of the unit.
   */
  char *a, *b, *c, *d;
  value_t v;
}
param_dihedral_t;

/* param_radius_t: structure for holding van der waals radii.
 */
typedef struct {
  /* @a: atom type for the radius entry.
   * @v: radius numerical value.
   */
  char *a;
  value_t v;
}
param_radius_t;

/* param_t: structure for holding molecular parameter information.
 */
typedef struct {
  /* @bonds: array of bond distances.
   * @n_bonds: number of bond distances.
   */
  param_bond_t *bonds;
  unsigned int n_bonds;

  /* @angles: array of two-bond angles.
   * @n_angles: number of two-bond angles.
   */
  param_angle_t *angles;
  unsigned int n_angles;

  /* @torsions: array of torsional dihedral angles.
   * @n_torsions: number of torsional dihedral angles.
   */
  param_dihedral_t *torsions;
  unsigned int n_torsions;

  /* @impropers: array of improper dihedral angles.
   * @n_impropers: number of improper dihedral angles.
   */
  param_dihedral_t *impropers;
  unsigned int n_impropers;

  /* @radii: array of van der waals radii.
   * @n_radii: number of van der waals radii.
   * @vdw_scale: radius scaling factor.
   */
  param_radius_t *radii;
  unsigned int n_radii;
  double vdw_scale;
}
param_t;

/* function declarations (param-alloc.c): */

param_t *param_new_from_file (const char *fname, double vdw_scale);

void param_free (param_t *par);

/* function declarations (param-add.c): */

int param_add_bond (param_t *par,
                    const char *a,
                    const char *b,
                    value_t length);

int param_add_angle (param_t *par,
                     const char *a,
                     const char *b,
                     const char *c,
                     value_t angle);

int param_add_torsion (param_t *par,
                       const char *a,
                       const char *b,
                       const char *c,
                       const char *d,
                       value_t angle);

int param_add_improper (param_t *par,
                        const char *a,
                        const char *b,
                        const char *c,
                        const char *d,
                        value_t angle);

int param_add_radius (param_t *par, const char *type, double sigma);

/* function declarations (param-get.c): */

value_t param_get_bond (param_t *par,
                        const char *a,
                        const char *b);

value_t param_get_angle (param_t *par,
                         const char *a,
                         const char *b,
                         const char *c);

value_t param_get_torsion (param_t *par,
                           const char *a,
                           const char *b,
                           const char *c,
                           const char *d);

value_t param_get_improper (param_t *par,
                            const char *a,
                            const char *b,
                            const char *c,
                            const char *d);

value_t param_get_radius (param_t *par, const char *type);

/* function declarations (param.c): */

int param_apply_all (param_t *par, peptide_t *P);

