
/* include the peptide header. */
#include "peptide.h"

/* peptide_atom_find(): lookup an atom in a peptide by its residue number
 * and atom name string.
 *
 * argument:
 *  @P: pointer to the peptide structure to access.
 *  @resid: residue sequence index of the atom.
 *  @name: name of the atom.
 *
 * returns: 
 *  index of the atom in the peptide's atom array, or -1 if no match.
 */
int peptide_atom_find (peptide_t *P, unsigned int resid, const char *name) {
  /* declare required variables:
   *  @i: atom loop counter.
   */
  unsigned int i;

  /* loop over the atoms of the peptide. */
  for (i = 0; i < P->n_atoms; i++) {
    /* check if the current atom is a match. */
    if (P->atoms[i].res_id == resid && strcmp(P->atoms[i].name, name) == 0)
      return (int) i;
  }

  /* return failure. */
  return -1;
}

/* peptide_atom_add(): add an atom to a peptide structure, given all
 * required atom information.
 *
 * arguments:
 *  @P: pointer to the peptide structure to modify.
 *  @resid: sequence index of the residue to which the atom belongs.
 *  @name: name of the atom, within the context of the residue.
 *  @type: type of the atom, within the global context of atom types.
 *  @mass: mass of the atom.
 *  @charge: partial charge of the atom.
 *  @radius: radius of the atom.
 *
 * returns:
 *  integer indicating whether (1) or not (0) the atom was successfully
 *  added to the peptide.
 */
int peptide_atom_add (peptide_t *P,
                      unsigned int resid,
                      const char *name,
                      const char *type,
                      double mass,
                      double charge,
                      double radius) {
  /* declare required variables:
   *  @i: atom index.
   */
  unsigned int i;

  /* check if the atom already exists. */
  if (peptide_atom_find(P, resid, name) >= 0)
    throw("atom %u.%s already exists", resid, name);

  /* increment the atom array length. */
  i = P->n_atoms;
  P->n_atoms++;

  /* reallocate the atom array. */
  P->atoms = (peptide_atom_t*)
    realloc(P->atoms, P->n_atoms * sizeof(peptide_atom_t));

  /* check if reallocation failed. */
  if (!P->atoms)
    throw("unable to reallocate atom array");

  /* store the residue index. */
  P->atoms[i].res_id = resid;

  /* store the string values (as duplicates). */
  P->atoms[i].name = strdup(name);
  P->atoms[i].type = strdup(type);

  /* store the numeric values. */
  P->atoms[i].mass = mass;
  P->atoms[i].charge = charge;
  P->atoms[i].radius = radius;

  /* return success. */
  return 1;
}

/* peptide_atom_modify(): modify an existing atom in a peptide structure,
 * or add it if none exists with the given name.
 *
 * arguments:
 *  see peptide_atom_add().
 *
 * returns:
 *  integer indicating whether (1) or not (0) the operation succeeded.
 */
int peptide_atom_modify (peptide_t *P,
                         unsigned int resid,
                         const char *name,
                         const char *type,
                         double mass,
                         double charge,
                         double radius) {
  /* declare required variables:
   *  @i: atom index.
   */
  int i;

  /* lookup the atom, and add it if none exists. */
  i = peptide_atom_find(P, resid, name);
  if (i < 0)
    return peptide_atom_add(P, resid, name, type, mass, charge, radius);

  /* check if a new type string was provided. */
  if (type && strlen(type)) {
    /* yes. free the existing type string. */
    if (P->atoms[i].type)
      free(P->atoms[i].type);

    /* store the new type string. */
    P->atoms[i].type = strdup(type);
  }

  /* store the new mass, if provided. */
  if (!isnan(mass))
    P->atoms[i].mass = mass;

  /* store the new charge, if provided. */
  if (!isnan(charge))
    P->atoms[i].charge = charge;

  /* store the new radius, if provided. */
  if (!isnan(radius))
    P->atoms[i].radius = radius;

  /* return success. */
  return 1;
}

/* peptide_atom_delete(): delete an existing atom from a peptide structure.
 *
 * arguments:
 *  @P: pointer to the peptide structure to modify.
 *  @resid: sequence index of the residue to which the atom belongs.
 *  @name: name of the atom, within the context of the residue.
 *
 * returns:
 *  integer indicating whether (1) or not (0) the operation succeeded.
 */
int peptide_atom_delete (peptide_t *P,
                         unsigned int resid,
                         const char *name) {
  /* declare required variables:
   *  @k, @ki: general-purpose indices.
   *  @i, @j: atom indices.
   */
  unsigned int j, k, ki;
  int i;

  /* lookup the atom, and return if none exists. */
  i = peptide_atom_find(P, resid, name);
  if (i < 0)
    return 1;

  /* free the atom name string. */
  if (P->atoms[i].name)
    free(P->atoms[i].name);

  /* free the atom type string. */
  if (P->atoms[i].type)
    free(P->atoms[i].type);

  /* get the index of the last atom. */
  j = P->n_atoms - 1;

  /* check if the atom swapping is required. */
  if ((unsigned int) i != j)
    memcpy(P->atoms + i, P->atoms + j, sizeof(peptide_atom_t));

  /* update atom indices in the bond array. */
  for (k = 0; k < P->n_bonds; k++) {
    for (ki = 0; ki < 2; ki++) {
      if (P->bonds[k].atom_id[ki] == j)
        P->bonds[k].atom_id[ki] = i;
    }
  }

  /* update atom indices in the angle array. */
  for (k = 0; k < P->n_angles; k++) {
    for (ki = 0; ki < 3; ki++) {
      if (P->angles[k].atom_id[ki] == j)
        P->angles[k].atom_id[ki] = i;
    }
  }

  /* update atom indices in the torsion array. */
  for (k = 0; k < P->n_torsions; k++) {
    for (ki = 0; ki < 4; ki++) {
      if (P->torsions[k].atom_id[ki] == j)
        P->torsions[k].atom_id[ki] = i;
    }
  }

  /* update atom indices in the improper array. */
  for (k = 0; k < P->n_impropers; k++) {
    for (ki = 0; ki < 4; ki++) {
      if (P->impropers[k].atom_id[ki] == j)
        P->impropers[k].atom_id[ki] = i;
    }
  }

  /* decrement the atom count and return success. */
  P->n_atoms--;
  return 1;
}

