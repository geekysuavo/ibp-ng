
/* include the molecular topology header. */
#include "topol.h"
#include "topol-add.h"

/* topol_is_connected(): determine whether a set of bonded atoms form a
 * connected tuple.
 *
 * this function modifies the contents of the input array of atoms in order
 * to determine connectivity. in the case of a connected tuple, the array
 * will be in-order, i.e.: (A-B, C-B, D-C) => (A-B, B-C, C-D).
 *
 * arguments:
 *  @ids: array of atom name strings to check for connectivity.
 *  @n: number of index pairs in the array.
 *
 * returns:
 *  integer indicating whether (1) or not (0) the indices can be rearranged
 *  (within each bonded pair) to form a connected tuple.
 */
int topol_is_connected (char **ids, unsigned int n) {
  /* declare required variables:
   *  @i: main loop counter.
   *  @idtemp: swap location.
   */
  unsigned int i;
  char *idtemp;

  /* return false in the trivial case of a single bond. */
  if (n <= 1)
    return 0;

  /* check the first pair explicitly. */
  if (strcmp(ids[1], ids[2]) && strcmp(ids[1], ids[3])) {
    /* swap the first pair. */
    idtemp = ids[0];
    ids[0] = ids[1];
    ids[1] = idtemp;
  }

  /* loop over the remaining pairs. */
  for (i = 1; i < n; i++) {
    /* check if the current pair is connected to the previous pair. */
    if (strcmp(ids[2 * i], ids[2 * i - 1]) == 0) {
      /* in-order connectivity. */
      continue;
    }
    else if (strcmp(ids[2 * i + 1], ids[2 * i - 1]) == 0) {
      /* out-of-order connectivity. */
      idtemp = ids[2 * i];
      ids[2 * i] = ids[2 * i + 1];
      ids[2 * i + 1] = idtemp;
    }
    else {
      /* no connectivity. */
      return 0;
    }
  }

  /* return true. */
  return 1;
}

/* topol_autogen_angles(): autogenerate all possible angle topology
 * entries based on available connectivities.
 *
 * arguments:
 *  @top: pointer to the main topology structure.
 *  @res: pointer to the residue topology structure to modify.
 *
 * returns:
 *  integer indicating whether (1) or not (0) the operation succeeded.
 */
int topol_autogen_angles (topol_t *top, topol_residue_t *res) {
  /* declare required variables:
   *  @i1, @i2: bond indices.
   *  @ids: atom name strings.
   */
  unsigned int i1, i2;
  char *ids[4];

  /* loop over all bonds. */
  for (i1 = 0; i1 < res->n_bonds; i1++) {
    /* set the first pair of atom names. */
    ids[0] = res->bonds[i1].atoms[0];
    ids[1] = res->bonds[i1].atoms[1];

    /* loop over all other bonds. */
    for (i2 = i1 + 1; i2 < res->n_bonds; i2++) {
      /* set the second pair of atom names. */
      ids[2] = res->bonds[i2].atoms[0];
      ids[3] = res->bonds[i2].atoms[1];

      /* check if the bonds form a connected triple of atoms. */
      if (topol_is_connected(ids, 2)) {
        /* attempt to add the angle to the topology structure. */
        if (!topol_add_angle(top, res->name,
                             ids[0], 0, ids[1], 0, ids[3], 0,
                             TOPOL_MODE_ADD))
          throw("unable to autogenerate angle %s-%s-%s on %s",
                ids[0], ids[1], ids[3], res->name);
      }
    }
  }

  /* return success. */
  return 1;
}

/* topol_autogen_torsions(): autogenerate all possible torsion topology
 * entries based on available connectivities.
 *
 * arguments:
 *  @top: pointer to the main topology structure.
 *  @res: pointer to the residue topology structure to modify.
 *
 * returns:
 *  integer indicating whether (1) or not (0) the operation succeeded.
 */
int topol_autogen_torsions (topol_t *top, topol_residue_t *res) {
  /* declare required variables:
   *  @i1, @i2, @i3: bond indices.
   */
  unsigned int i1, i2, i3;
  char *ids[6];

  /* loop over all bonds. */
  for (i1 = 0; i1 < res->n_bonds; i1++) {
    /* set the first pair of atom pointers. */
    ids[0] = res->bonds[i1].atoms[0];
    ids[1] = res->bonds[i1].atoms[1];

    /* loop over all other bonds. */
    for (i2 = i1 + 1; i2 < res->n_bonds; i2++) {
      /* set the second pair of atom pointers. */
      ids[2] = res->bonds[i2].atoms[0];
      ids[3] = res->bonds[i2].atoms[1];

      /* loop once more over all remaining bonds. */
      for (i3 = i2 + 1; i3 < res->n_bonds; i3++) {
        /* set the third pair of atom pointers. */
        ids[4] = res->bonds[i3].atoms[0];
        ids[5] = res->bonds[i3].atoms[1];

        /* check if the bonds form a connected quadruple of atoms. */
        if (topol_is_connected(ids, 3)) {
          /* attempt to add the torsion to the topology structure. */
          if (!topol_add_torsion(top, res->name,
                                 ids[0], 0, ids[1], 0,
                                 ids[3], 0, ids[5], 0,
                                 TOPOL_MODE_ADD))
            throw("unable to autogenerate dihedral %s-%s-%s-%s on %s",
                  ids[0], ids[1], ids[3], ids[5], res->name);
        }
      }
    }
  }

  /* return success. */
  return 1;
}

/* topol_autogen(): complete the last residue entry of a topology
 * structure by autogenerating any required angles and torsions.
 *
 * arguments:
 *  @top: pointer to the topology structure to modify.
 *
 * returns:
 *  integer indicating whether (1) or not (0) the operation succeeded.
 */
int topol_autogen (topol_t *top) {
  /* declare required variables:
   *  @res: residue topology structure pointer.
   */
  topol_residue_t *res;

  /* check that at least one residue exists. */
  if (top->n_res == 0)
    throw("no residue topology exists");

  /* get the residue structure pointer. */
  res = top->res + (top->n_res - 1);

  /* skip patch residues. */
  if (res->patch)
    return 1;

  /* if enabled, autogenerate angles. */
  if (top->auto_angles && !topol_autogen_angles(top, res))
    throw("unable to autogenerate angles for topology '%s'",
          res->name);

  /* if enabled, autogenerate torsions. */
  if (top->auto_dihedrals && !topol_autogen_torsions(top, res))
    throw("unable to autogenerate dihedrals for topology '%s'",
          res->name);

  /* return success. */
  return 1;
}

