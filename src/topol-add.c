
/* include the molecular topology header. */
#include "topol.h"

/* topol_find_residue(): lookup the index of a named residue topology entry.
 *
 * arguments:
 *  @top: pointer to the topology structure to access.
 *  @resname: string name of the residue to find.
 *
 * returns:
 *  index of the residue, or -1 if no match was found.
 */
int topol_find_residue (topol_t *top, const char *resname) {
  /* declare required variables:
   *  @i: loop counter.
   */
  unsigned int i;

  /* search the array for the query name string. */
  for (i = 0; top && i < top->n_res; i++) {
    /* return if the current array entry is a match. */
    if (strcmp(top->res[i].name, resname) == 0)
      return (int) i;
  }

  /* return failure. */
  return -1;
}

/* topol_get_mass(): lookup the mass of a particular atom type.
 *
 * arguments:
 *  @top: pointer to the topology structure to access.
 *  @type: string type of the query atom.
 *
 * returns:
 *  mass of the atom type, or NAN if no such atom type was found.
 */
double topol_get_mass (topol_t *top, const char *type) {
  /* declare required variables:
   *  @i: loop counter.
   */
  unsigned int i;

  /* search the array for the query type string. */
  for (i = 0; top && i < top->n_mass; i++) {
    /* return if the current array entry is a match. */
    if (strcmp(top->mass[i].type, type) == 0)
      return top->mass[i].mass;
  }

  /* return failure. */
  return NAN;
}

/* * * * functions declared in topol.h begin here: * * * */

/* topol_add_mass(): add a mass entry into a topology structure.
 *
 * arguments:
 *  @top: pointer to the topology structure to modify.
 *  @type: string type of the new atom mass entry.
 *  @mass: numerical value of the entry.
 *
 * returns:
 *  integer indicating whether (1) or not (0) the operation succeeded.
 */
int topol_add_mass (topol_t *top, const char *type, double mass) {
  /* declare required variables:
   *  @i: index of the new array entry.
   */
  unsigned int i;

  /* increment the array size. */
  i = top->n_mass;
  top->n_mass++;

  /* reallocate the array. */
  top->mass = (topol_mass_t*)
    realloc(top->mass, top->n_mass * sizeof(topol_mass_t));

  /* check if reallocation failed. */
  if (!top->mass)
    throw("unable to reallocate mass array");

  /* store the values of the new entry. */
  top->mass[i].type = strdup(type);
  top->mass[i].mass = mass;

  /* return success. */
  return 1;
}

/* topol_add_residue(): add a residue entry into a topology structure.
 *
 * arguments:
 *  @top: pointer to the topology structure to modify.
 *  @name: string type of the new residue entry.
 *  @patch: patch-status of the residue.
 *
 * returns:
 *  integer indicating whether (1) or not (0) the operation succeeded.
 */
int topol_add_residue (topol_t *top, const char *name, unsigned int patch) {
  /* declare required variables:
   *  @i: index of the new array entry.
   */
  unsigned int i;

  /* check that the residue does not already exist. */
  if (topol_find_residue(top, name) >= 0)
    throw("topology already contains residue '%s'", name);

  /* increment the array size. */
  i = top->n_res;
  top->n_res++;

  /* reallocate the array. */
  top->res = (topol_residue_t*)
    realloc(top->res, top->n_res * sizeof(topol_residue_t));

  /* check if reallocation failed. */
  if (!top->res)
    throw("unable to reallocate residue array");

  /* store the name of the new entry. */
  top->res[i].name = strdup(name);

  /* initialize the atom array. */
  top->res[i].atoms = NULL;
  top->res[i].n_atoms = 0;

  /* initialize the bond array. */
  top->res[i].bonds = NULL;
  top->res[i].n_bonds = 0;

  /* initialize the angle array. */
  top->res[i].angles = NULL;
  top->res[i].n_angles = 0;

  /* initialize the torsion array. */
  top->res[i].torsions = NULL;
  top->res[i].n_torsions = 0;

  /* initialize the improper array. */
  top->res[i].impropers = NULL;
  top->res[i].n_impropers = 0;

  /* store the patch status. */
  top->res[i].patch = patch;

  /* return success. */
  return 1;
}

/* topol_add_atom(): add an atom entry into the last residue of a
 * topology structure.
 *
 * if @resname is passed as NULL, then the atom entry will be added
 * to the last residue in the structure.
 *
 * arguments:
 *  @top: pointer to the topology structure to modify.
 *  @resname: string name of the residue to modify.
 *  @name: string name of the atom in the residue.
 *  @type: string type of the atom to add.
 *  @charge: partial charge on the atom.
 *  @mode: topology mode of the atom.
 *  @off: residue offset of the atom.
 *
 * returns:
 *  integer indicating whether (1) or not (0) the operation succeeded.
 */
int topol_add_atom (topol_t *top,
                    const char *resname,
                    const char *name,
                    const char *type,
                    double charge,
                    topol_mode_t mode, int off) {
  /* declare required variables:
   *  @res: pointer to the residue topology to modify.
   *  @i: index of the new array entry.
   *  @mass: mass of the new atom.
   */
  topol_residue_t *res;
  double mass;
  int i;

  /* look up the mass of the atom. */
  mass = topol_get_mass(top, type);

  /* check that a valid mass was identified. */
  if (isnan(mass))
    throw("topology contains no mass definition for '%s'", type);

  /* check that the residue array is allocated. */
  if (!top->n_res)
    throw("topology contains no residues");

  /* determine the residue index to be modified. */
  if (resname)
    i = topol_find_residue(top, resname);
  else
    i = top->n_res - 1;

  /* check that the residue index is valid. */
  if (i < 0 || (unsigned int) i >= top->n_res)
    throw("residue index out of bounds [0,%u]", top->n_res - 1);

  /* gain access to the residue pointer. */
  res = top->res + i;

  /* increment the atom array size. */
  i = res->n_atoms;
  res->n_atoms++;

  /* reallocate the array. */
  res->atoms = (topol_atom_t*)
    realloc(res->atoms, res->n_atoms * sizeof(topol_atom_t));

  /* check if reallocation failed. */
  if (!res->atoms)
    throw("unable to reallocate atom array");

  /* store the strings of the new atom entry. */
  res->atoms[i].name = strdup(name);
  res->atoms[i].type = strdup(type);

  /* store the mass and charge of the new atom entry. */
  res->atoms[i].mass = mass;
  res->atoms[i].charge = charge;

  /* store the mode and offset of the new atom entry. */
  res->atoms[i].mode = mode;
  res->atoms[i].off = off;

  /* return success. */
  return 1;
}

/* topol_add_bond(): add a bond entry into a particular residue of a
 * topology structure.
 *
 * if @resname is passed as NULL, then the bond entry will be added
 * to the last residue in the structure.
 *
 * arguments:
 *  @top: pointer to the topology structure to modify.
 *  @resname: string name of the residue to modify.
 *  @a: string name of the first atom in the bond.
 *  @b: string name of the second atom in the bond.
 *  @aoff: residue offset of the first atom.
 *  @boff: residue offset of the second atom.
 *  @mode:  topology mode of the bond.
 *
 * returns:
 *  integer indicating whether (1) or not (0) the operation succeeded.
 */
int topol_add_bond (topol_t *top,
                    const char *resname,
                    const char *a, int aoff,
                    const char *b, int boff,
                    topol_mode_t mode) {
  /* declare required variables:
   *  @res: pointer to the residue topology to modify.
   *  @i: index of the new array entry.
   */
  topol_residue_t *res;
  int i;

  /* check that the residue array is allocated. */
  if (!top->n_res)
    throw("topology contains no residues");

  /* determine the residue index to be modified. */
  if (resname)
    i = topol_find_residue(top, resname);
  else
    i = top->n_res - 1;

  /* check that the residue index is valid. */
  if (i < 0 || (unsigned int) i >= top->n_res)
    throw("residue index out of bounds [0,%u]", top->n_res - 1);

  /* gain access to the residue pointer. */
  res = top->res + i;

  /* increment the bond array size. */
  i = res->n_bonds;
  res->n_bonds++;

  /* reallocate the bond array. */
  res->bonds = (topol_bond_t*)
    realloc(res->bonds, res->n_bonds * sizeof(topol_bond_t));

  /* check if reallocation failed. */
  if (!res->bonds)
    throw("unable to reallocate bond array");

  /* store the atom names in the new bond entry. */
  res->bonds[i].atoms[0] = strdup(a);
  res->bonds[i].atoms[1] = strdup(b);

  /* store the residue offsets in the new bond entry. */
  res->bonds[i].off[0] = aoff;
  res->bonds[i].off[1] = boff;

  /* store the topology of the new bond entry. */
  res->bonds[i].mode = mode;

  /* return success. */
  return 1;
}

/* topol_add_angle(): add an angle entry into a particular residue of a
 * topology structure.
 *
 * if @resname is passed as NULL, then the angle entry will be added
 * to the last residue in the structure.
 *
 * arguments:
 *  @top: pointer to the topology structure to modify.
 *  @resname: string name of the residue to modify.
 *  @a: string name of the first atom in the angle.
 *  @b: string name of the second atom in the angle.
 *  @c: string name of the third atom in the angle.
 *  @aoff: residue offset of the first atom.
 *  @boff: residue offset of the second atom.
 *  @coff: residue offset of the third atom.
 *  @mode: topology mode of the angle.
 *
 * returns:
 *  integer indicating whether (1) or not (0) the operation succeeded.
 */
int topol_add_angle (topol_t *top,
                     const char *resname,
                     const char *a, int aoff,
                     const char *b, int boff,
                     const char *c, int coff,
                     topol_mode_t mode) {
  /* declare required variables:
   *  @res: pointer to the residue topology to modify.
   *  @i: index of the new array entry.
   */
  topol_residue_t *res;
  int i;

  /* check that the residue array is allocated. */
  if (!top->n_res)
    throw("topology contains no residues");

  /* determine the residue index to be modified. */
  if (resname)
    i = topol_find_residue(top, resname);
  else
    i = top->n_res - 1;

  /* check that the residue index is valid. */
  if (i < 0 || (unsigned int) i >= top->n_res)
    throw("residue index out of bounds [0,%u]", top->n_res - 1);

  /* gain access to the residue pointer. */
  res = top->res + i;

  /* increment the angle array size. */
  i = res->n_angles;
  res->n_angles++;

  /* reallocate the angle array. */
  res->angles = (topol_angle_t*)
    realloc(res->angles, res->n_angles * sizeof(topol_angle_t));

  /* check if reallocation failed. */
  if (!res->angles)
    throw("unable to reallocate angle array");

  /* store the atom names in the new angle entry. */
  res->angles[i].atoms[0] = strdup(a);
  res->angles[i].atoms[1] = strdup(b);
  res->angles[i].atoms[2] = strdup(c);

  /* store the residue offsets in the new angle entry. */
  res->angles[i].off[0] = aoff;
  res->angles[i].off[1] = boff;
  res->angles[i].off[2] = coff;

  /* store the topology mode of the new angle entry. */
  res->angles[i].mode = mode;

  /* return success. */
  return 1;
}

/* topol_add_torsion(): add a torsion entry into a particular residue of a
 * topology structure.
 *
 * if @resname is passed as NULL, then the torsion entry will be added
 * to the last residue in the structure.
 *
 * arguments:
 *  @top: pointer to the topology structure to modify.
 *  @resname: string name of the residue to modify.
 *  @a: string name of the first atom in the entry.
 *  @b: string name of the second atom in the entry.
 *  @c: string name of the third atom in the entry.
 *  @d: string name of the fourth atom in the entry.
 *  @aoff: residue offset of the first atom.
 *  @boff: residue offset of the second atom.
 *  @coff: residue offset of the third atom.
 *  @doff: residue offset of the fourth atom.
 *  @mode: topology mode of the angle.
 *
 * returns:
 *  integer indicating whether (1) or not (0) the operation succeeded.
 */
int topol_add_torsion (topol_t *top,
                       const char *resname,
                       const char *a, int aoff,
                       const char *b, int boff,
                       const char *c, int coff,
                       const char *d, int doff,
                       topol_mode_t mode) {
  /* declare required variables:
   *  @res: pointer to the residue topology to modify.
   *  @i: index of the new array entry.
   */
  topol_residue_t *res;
  int i;

  /* check that the residue array is allocated. */
  if (!top->n_res)
    throw("topology contains no residues");

  /* determine the residue index to be modified. */
  if (resname)
    i = topol_find_residue(top, resname);
  else
    i = top->n_res - 1;

  /* check that the residue index is valid. */
  if (i < 0 || (unsigned int) i >= top->n_res)
    throw("residue index out of bounds [0,%u]", top->n_res - 1);

  /* gain access to the residue pointer. */
  res = top->res + i;

  /* increment the torsion array size. */
  i = res->n_torsions;
  res->n_torsions++;

  /* reallocate the torsion array. */
  res->torsions = (topol_dihedral_t*)
    realloc(res->torsions, res->n_torsions * sizeof(topol_dihedral_t));

  /* check if reallocation failed. */
  if (!res->torsions)
    throw("unable to reallocate torsion array");

  /* store the atom names in the new array entry. */
  res->torsions[i].atoms[0] = strdup(a);
  res->torsions[i].atoms[1] = strdup(b);
  res->torsions[i].atoms[2] = strdup(c);
  res->torsions[i].atoms[3] = strdup(d);

  /* store the residue offsets in the new array entry. */
  res->torsions[i].off[0] = aoff;
  res->torsions[i].off[1] = boff;
  res->torsions[i].off[2] = coff;
  res->torsions[i].off[3] = doff;

  /* store the topology mode of the new array entry. */
  res->torsions[i].mode = mode;

  /* return success. */
  return 1;
}

/* topol_add_improper(): add an improper entry into a particular residue
 * of a topology structure.
 *
 * if @resname is passed as NULL, then the improper entry will be added
 * to the last residue in the structure.
 *
 * arguments:
 *  @top: pointer to the topology structure to modify.
 *  @resname: string name of the residue to modify.
 *  @a: string name of the first atom in the entry.
 *  @b: string name of the second atom in the entry.
 *  @c: string name of the third atom in the entry.
 *  @d: string name of the fourth atom in the entry.
 *  @d: string name of the fourth atom in the entry.
 *  @aoff: residue offset of the first atom.
 *  @boff: residue offset of the second atom.
 *  @coff: residue offset of the third atom.
 *  @doff: residue offset of the fourth atom.
 *
 * returns:
 *  integer indicating whether (1) or not (0) the operation succeeded.
 */
int topol_add_improper (topol_t *top,
                        const char *resname,
                        const char *a, int aoff,
                        const char *b, int boff,
                        const char *c, int coff,
                        const char *d, int doff,
                        topol_mode_t mode) {
  /* declare required variables:
   *  @res: pointer to the residue topology to modify.
   *  @i: index of the new array entry.
   */
  topol_residue_t *res;
  int i;

  /* check that the residue array is allocated. */
  if (!top->n_res)
    throw("topology contains no residues");

  /* determine the residue index to be modified. */
  if (resname)
    i = topol_find_residue(top, resname);
  else
    i = top->n_res - 1;

  /* check that the residue index is valid. */
  if (i < 0 || (unsigned int) i >= top->n_res)
    throw("residue index out of bounds [0,%u]", top->n_res - 1);

  /* gain access to the residue pointer. */
  res = top->res + i;

  /* increment the improper array size. */
  i = res->n_impropers;
  res->n_impropers++;

  /* reallocate the improper array. */
  res->impropers = (topol_dihedral_t*)
    realloc(res->impropers, res->n_impropers * sizeof(topol_dihedral_t));

  /* check if reallocation failed. */
  if (!res->impropers)
    throw("unable to reallocate improper array");

  /* store the atom names in the new array entry. */
  res->impropers[i].atoms[0] = strdup(a);
  res->impropers[i].atoms[1] = strdup(b);
  res->impropers[i].atoms[2] = strdup(c);
  res->impropers[i].atoms[3] = strdup(d);

  /* store the residue offsets in the new array entry. */
  res->impropers[i].off[0] = aoff;
  res->impropers[i].off[1] = boff;
  res->impropers[i].off[2] = coff;
  res->impropers[i].off[3] = doff;

  /* store the topology mode of the new array entry. */
  res->impropers[i].mode = mode;

  /* return success. */
  return 1;
}

