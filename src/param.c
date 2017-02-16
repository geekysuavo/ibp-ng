
/* include the molecular parameters header. */
#include "param.h"

/* param_check_angles(): ensure that all mutually adjacent triples of
 * angles (i.e. angles forming three planes that intersect at a
 * common vertex) sum to no more than 360 degrees.
 *
 * this resolves the issues encountered in prior versions of ibp, where
 * computed values of cosine(omega) were out of bounds. however, the
 * user is still presented with a warning message in order to provide
 * more acceptable means of permanently fixing the parameter file.
 *
 * arguments:
 *  @P: pointer to the peptide structure to modify.
 */
void param_check_angles (peptide_t *P) {
  /* declare required variables:
   *  @A, @B, @C, @D: atoms defining the trihedron/tetrahedron.
   *  @a1, @a2, @a3: peptide angle pointers.
   *  @i1, @i2, @i3: peptide angle array indices.
   *  @idx: index of the fourth angle in the set.
   *  @sum: sum of the three angles.
   */
  peptide_angle_t *a1, *a2, *a3;
  unsigned int i1, i2, i3;
  value_t sum;

  /* loop over the set of angles. */
  for (i1 = 0; i1 < P->n_angles; i1++) {
    /* obtain a pointer to the first angle. */
    a1 = P->angles + i1;

    /* loop again over the angles. */
    for (i2 = 0; i2 < P->n_angles; i2++) {
      /* obtain a pointer to the second angle. */
      a2 = P->angles + i2;

      /* skip the first angle. */
      if (i2 == i1)
        continue;

      /* skip angles with non-common central vertices. */
      if (a2->atom_id[1] != a1->atom_id[1])
        continue;

      /* skip angles that cannot yield a trihedron/tetrahedron. */
      if (a2->atom_id[0] != a1->atom_id[0] &&
          a2->atom_id[2] == a1->atom_id[0])
        continue;

      /* loop one last time over the angles. */
      for (i3 = 0; i3 < P->n_angles; i3++) {
        /* obtain a pointer to the third angle. */
        a3 = P->angles + i3;

        /* skip the first and second angles. */
        if (i3 == i2 || i3 == i1)
          continue;

        /* skip angles with non-common central vertices. */
        if (a3->atom_id[1] != a1->atom_id[1])
          continue;

        /* skip angles that cannot complete the trihedron/tetrahedron. */
        if (a3->atom_id[0] != a1->atom_id[2] &&
            a3->atom_id[2] != a1->atom_id[2])
          continue;

        /* compute the sum of the three angles. */
        sum = value_add(value_add(a1->ang, a2->ang), a3->ang);

        /* check if the sum is invalid. */
        if (sum.u >= 360.0) {
          /* rescale the angles. */
          a1->ang = value_scal(a1->ang, 360.0 / sum.u);
          a2->ang = value_scal(a2->ang, 360.0 / sum.u);
          a3->ang = value_scal(a3->ang, 360.0 / sum.u);
        }
      }
    }
  }
}

/* param_apply_all(): add parameter information to a peptide structure.
 *
 * arguments:
 *  @par: pointer to the parameter structure to access.
 *  @P: pointer to the peptide structure to modify.
 *
 * returns:
 *  integer indicating whether (1) or not (0) the operation succeeded.
 */
int param_apply_all (param_t *par, peptide_t *P) {
  /* declare required variables:
   *  @A, @B, @C, @D: peptide atom structure pointers.
   *  @i: general-purpose array index.
   *  @val: generalized value storage.
   */
  peptide_atom_t *A, *B, *C, *D;
  unsigned int i;
  value_t val;

  /* loop over the atoms in the peptide. */
  for (i = 0; i < P->n_atoms; i++) {
    /* lookup the atomic radius. */
    val = param_get_radius(par, P->atoms[i].type);

    /* validate the parameter. */
    if (value_is_undefined(val))
      throw("atom type %s has no nonbonded parameter entry",
            P->atoms[i].type);

    /* store the parameter. */
    P->atoms[i].radius = val.l;
  }

  /* loop over the bonds in the peptide. */
  for (i = 0; i < P->n_bonds; i++) {
    /* get the atoms. */
    A = P->atoms + P->bonds[i].atom_id[0];
    B = P->atoms + P->bonds[i].atom_id[1];

    /* lookup the parameter. */
    val = param_get_bond(par, A->type, B->type);

    /* reverse-lookup the parameter. */
    if (value_is_undefined(val))
      val = param_get_bond(par, B->type, A->type);

    /* validate the parameter. */
    if (value_is_undefined(val))
      throw("atom types %s,%s have no bond parameter entry",
            A->type, B->type);

    /* store the parameter. */
    P->bonds[i].len = val;
  }

  /* loop over the angles in the peptide. */
  for (i = 0; i < P->n_angles; i++) {
    /* get the atoms. */
    A = P->atoms + P->angles[i].atom_id[0];
    B = P->atoms + P->angles[i].atom_id[1];
    C = P->atoms + P->angles[i].atom_id[2];

    /* lookup the parameter. */
    val = param_get_angle(par, A->type, B->type, C->type);

    /* reverse-lookup the parameter. */
    if (value_is_undefined(val))
      val = param_get_angle(par, C->type, B->type, A->type);

    /* validate the parameter. */
    if (value_is_undefined(val))
      throw("atom types %s,%s,%s have no angle parameter entry",
            A->type, B->type, C->type);

    /* store the parameter. */
    P->angles[i].ang = val;
  }

  /* loop over the torsions in the peptide. */
  for (i = 0; i < P->n_torsions; i++) {
    /* get the atoms. */
    A = P->atoms + P->torsions[i].atom_id[0];
    B = P->atoms + P->torsions[i].atom_id[1];
    C = P->atoms + P->torsions[i].atom_id[2];
    D = P->atoms + P->torsions[i].atom_id[3];

    /* lookup the parameter. */
    val = param_get_torsion(par, A->type, B->type, C->type, D->type);

    /* reverse-lookup the parameter. */
    if (value_is_undefined(val))
      val = param_get_torsion(par, D->type, C->type, B->type, A->type);

    /* validate the parameter. */
    if (value_is_undefined(val))
      throw("atom types %s,%s,%s,%s have no dihedral parameter entry",
            A->type, B->type, C->type, D->type);

    /* store the parameter. */
    P->torsions[i].ang = val;
  }

  /* loop over the impropers in the peptide. */
  for (i = 0; i < P->n_impropers; i++) {
    /* get the atoms. */
    A = P->atoms + P->impropers[i].atom_id[0];
    B = P->atoms + P->impropers[i].atom_id[1];
    C = P->atoms + P->impropers[i].atom_id[2];
    D = P->atoms + P->impropers[i].atom_id[3];

    /* lookup the parameter. */
    val = param_get_improper(par, A->type, B->type, C->type, D->type);

    /* reverse-lookup the parameter. */
    if (value_is_undefined(val))
      val = param_get_improper(par, D->type, C->type, B->type, A->type);

    /* validate the parameter. */
    if (value_is_undefined(val))
      throw("atom types %s,%s,%s,%s have no improper parameter entry",
            A->type, B->type, C->type, D->type);

    /* store the parameter. */
    P->impropers[i].ang = val;
  }

  /* perform an angle check of the constructed peptide. */
  param_check_angles(P);

  /* return success. */
  return 1;
}

