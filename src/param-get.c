
/* include the molecular parameters header. */
#include "param.h"

/* param_get_bond(): lookup the bond length between two particular
 * atom types.
 *
 * arguments:
 *  @par: pointer to the parameter structure to access.
 *  @a: string name of the first atom.
 *  @b: string name of the second atom.
 *
 * returns:
 *  length of the bond, or an undefined value if no such bond was found.
 */
value_t param_get_bond (param_t *par,
                        const char *a,
                        const char *b) {
  /* declare required variables:
   *  @i: loop counter.
   */
  unsigned int i;

  /* search the array for the query atom strings. */
  for (i = 0; par && i < par->n_bonds; i++) {
    /* check if the current array entry is a match. */
    if (strcmp(par->bonds[i].a, a) == 0 &&
        strcmp(par->bonds[i].b, b) == 0) {
      /* yes. return success. */
      return par->bonds[i].v;
    }
  }

  /* return failure. */
  return value_undefined();
}

/* param_get_angle(): lookup the two-bond angle between three particular
 * atom types.
 *
 * arguments:
 *  @par: pointer to the parameter structure to access.
 *  @a: string name of the first atom.
 *  @b: string name of the second atom.
 *  @c: string name of the third atom.
 *
 * returns:
 *  value of the angle, or an undefined value if no such angle was found.
 */
value_t param_get_angle (param_t *par,
                         const char *a,
                         const char *b,
                         const char *c) {
  /* declare required variables:
   *  @i: loop counter.
   */
  unsigned int i;

  /* search the array for the query atom strings. */
  for (i = 0; par && i < par->n_angles; i++) {
    /* check if the current array entry is a match. */
    if (strcmp(par->angles[i].a, a) == 0 &&
        strcmp(par->angles[i].b, b) == 0 &&
        strcmp(par->angles[i].c, c) == 0) {
      /* yes. return success. */
      return par->angles[i].v;
    }
  }

  /* return failure. */
  return value_undefined();
}

/* param_get_torsion(): lookup the three-bond angle between four particular
 * atom types.
 *
 * arguments:
 *  @par: pointer to the parameter structure to access.
 *  @a: string name of the first atom.
 *  @b: string name of the second atom.
 *  @c: string name of the third atom.
 *  @d: string name of the fourth atom.
 *
 * returns:
 *  value of the angle, or an undefined value if no such angle was found.
 */
value_t param_get_torsion (param_t *par,
                           const char *a,
                           const char *b,
                           const char *c,
                           const char *d) {
  /* declare required variables:
   *  @i: loop counter.
   */
  unsigned int i;

  /* search the array for the query atom strings. */
  for (i = 0; par && i < par->n_torsions; i++) {
    /* check if the current array entry is a match. */
    if (strcmp(par->torsions[i].a, a) == 0 &&
        strcmp(par->torsions[i].b, b) == 0 &&
        strcmp(par->torsions[i].c, c) == 0 &&
        strcmp(par->torsions[i].d, d) == 0) {
      /* yes. return success. */
      return par->torsions[i].v;
    }
  }

  /* return failure. */
  return value_undefined();
}

/* param_get_improper(): lookup the three-bond angle between four particular
 * atom types.
 *
 * arguments:
 *  @par: pointer to the parameter structure to access.
 *  @a: string name of the first atom.
 *  @b: string name of the second atom.
 *  @c: string name of the third atom.
 *  @d: string name of the fourth atom.
 *
 * returns:
 *  value of the angle, or an undefined value if no such angle was found.
 */
value_t param_get_improper (param_t *par,
                            const char *a,
                            const char *b,
                            const char *c,
                            const char *d) {
  /* declare required variables:
   *  @i: loop counter.
   */
  unsigned int i;

  /* search the array for the query atom strings. */
  for (i = 0; par && i < par->n_impropers; i++) {
    /* check if the current array entry is a match. */
    if (strcmp(par->impropers[i].a, a) == 0 &&
        strcmp(par->impropers[i].b, b) == 0 &&
        strcmp(par->impropers[i].c, c) == 0 &&
        strcmp(par->impropers[i].d, d) == 0) {
      /* yes. return success. */
      return par->impropers[i].v;
    }
  }

  /* return failure. */
  return value_undefined();
}

/* param_get_radius(): lookup the van der waals radius of a particular
 * atom type.
 *
 * arguments:
 *  @par: pointer to the parameter structure to access.
 *  @type: string type of the atom in the entry.
 *
 * returns:
 *  value of the radius, or an undefined value if no such atom type was found.
 */
value_t param_get_radius (param_t *par, const char *type) {
  /* declare required variables:
   *  @i: loop counter.
   */
  unsigned int i;

  /* search the array for the query atom type. */
  for (i = 0; par && i < par->n_radii; i++) {
    /* check if the current array entry is a match. */
    if (strcmp(par->radii[i].a, type) == 0) {
      /* yes. return success. */
      return par->radii[i].v;
    }
  }

  /* return failure. */
  return value_undefined();
}

