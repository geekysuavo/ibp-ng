
/* include the assignment header. */
#include "assign.h"

/* include the relevant peptide headers. */
#include "peptide-bonds.h"
#include "peptide-impropers.h"

/* function declarations from parsing code: */
int assign_parse(FILE *fh, peptide_t *P);

/* assign_free(): free allocated memory associated with an assignment set.
 *
 * arguments:
 *  @set: pointer to the atom set to free.
 */
void assign_free (assign_set_t *set) {
  /* return if the set pointer is null. */
  if (!set) return;

  /* free the index array of the set. */
  if (set->v)
    free(set->v);

  /* free the set pointer. */
  free(set);
}

/* assign_is_wild(): check if an assignment name string contains any
 * wildcard characters.
 *
 * arguments:
 *  @name: string to check for wildcards.
 *
 * returns:
 *  integer indicating whether (1) or not (0) wildcards were found.
 */
int assign_is_wild (const char *name) {
  /* return true if any wildcard is found. */
  if (strchr(name, '*') ||
      strchr(name, '%') ||
      strchr(name, '#') ||
      strchr(name, '+'))
    return 1;

  /* return false. */
  return 0;
}

/* assign_is_empty(): check if an assignment set contains no atoms.
 *
 * arguments:
 *  @set: pointer to the assignment set to check.
 *  @idx: index of the assignment set in the expression.
 *
 * returns:
 *  integer indicating whether (1) or not (0) the input set is empty.
 */
int assign_is_empty (assign_set_t *set, unsigned int idx) {
  /* fail if the set is empty. */
  if (set->n == 0) {
    /* raise an exception and return true. */
    raise("atom selector #%u is empty", idx);
    return 1;
  }

  /* return false. */
  return 0;
}

/* assign_is_ambiguous(): check if an assignment set contains multiple atoms.
 *
 * arguments:
 *  @set: pointer to the assignment set to check.
 *  @idx: index of the assignment set in the expression.
 *
 * returns:
 *  integer indicating whether (1) or not (0) the input set is ambiguous.
 */
int assign_is_ambiguous (assign_set_t *set, unsigned int idx) {
  /* declare required variables:
   *  @msg, @pmsg: verbose message string and pointer.
   *  @atom: atom structure pointer.
   *  @i: atom loop counter.
   */
  char *msg, *pmsg;
  peptide_atom_t *atom;
  unsigned int i;

  /* fail if the set contains multiple atoms. */
  if (set->n > 1) {
    /* raise an exception. */
    raise("atom selector #%u is ambiguous", idx);

    /* allocate a message string. */
    msg = (char*) malloc(128 * set->n * sizeof(char));
    if (!msg)
      return 1;

    /* initialize the message string. */
    pmsg = msg;

    /* loop over the matching atoms. */
    for (i = 0; i < set->n; i++) {
      /* get the atom pointer. */
      atom = set->P->atoms + set->v[i];

      /* add the atom information to the message string. */
      pmsg += sprintf(pmsg, "\n   | residue %s%u atom %s (%s)",
                      resid_get_code3(set->P->res[atom->res_id]),
                      atom->res_id + 1,
                      atom->name,
                      atom->type);
    }

    /* output and free the message string. */
    info("selector #%u matches multiple atoms:%s", idx, msg);
    free(msg);

    /* return true. */
    return 1;
  }

  /* return false. */
  return 0;
}

/* assign_none(): construct an empty assignment set.
 *
 * arguments:
 *  @P: pointer to the peptide to reference set operations against.
 *
 * returns:
 *  pointer to a newly allocated assignment set, or NULL on failure.
 */
assign_set_t *assign_none (peptide_t *P) {
  /* declare required variables:
   *  @set: pointer to the output set.
   */
  assign_set_t *set;

  /* check that the peptide structure pointer is valid. */
  if (!P) {
    /* raise an exception and return null. */
    raise("peptide structure pointer is null");
    return NULL;
  }

  /* check that the peptide structure contains atoms. */
  if (P->n_atoms == 0) {
    /* raise an exception and return null. */
    raise("peptide structure contains no atoms");
    return NULL;
  }

  /* allocate a set structure pointer. */
  set = (assign_set_t*) malloc(sizeof(assign_set_t));

  /* check if allocation failed. */
  if (!set) {
    /* raise an exception and return null. */
    raise("unable to allocate assignment structure pointer");
    return NULL;
  }

  /* set the maximum and current atom counts. */
  set->n_max = P->n_atoms;
  set->n = 0;

  /* store the peptide structure pointer. */
  set->P = P;

  /* allocate the atom index array. */
  set->v = (unsigned int*) malloc(set->n_max * sizeof(unsigned int));

  /* check if allocation failed. */
  if (!set->v) {
    /* raise an exception and return null. */
    raise("unable to allocate assignment index array");
    assign_free(set);
    return NULL;
  }

  /* return the new set pointer. */
  return set;
}

/* assign_all(): construct a complete assignment set.
 *
 * arguments:
 *  @P: pointer to the peptide to reference set operations against.
 *
 * returns:
 *  pointer to a newly allocated assignment set, or NULL on failure.
 */
assign_set_t *assign_all (peptide_t *P) {
  /* declare required variables:
   *  @set: pointer to the output set.
   *  @i: atom index loop counter.
   */
  assign_set_t *set;
  unsigned int i;

  /* allocate a new assignment set. */
  set = assign_none(P);
  if (!set)
    return NULL;

  /* update the atom count. */
  set->n = set->n_max;

  /* store all atom indices of the peptide in the assignment set. */
  for (i = 0; i < set->n; i++)
    set->v[i] = i;

  /* return the new set pointer. */
  return set;
}

/* assign_atomid(): construct an assignment set for a single indexed atom.
 *
 * arguments:
 *  @P: pointer to the peptide to reference set operations against.
 *  @id: atom index to include in the set.
 *
 * returns:
 *  pointer to the newly allocated assignment set, or NULL on failure.
 */
assign_set_t *assign_atomid (peptide_t *P, unsigned int id) {
  /* declare required variables:
   *  @set: pointer to the output set.
   */
  assign_set_t *set;

  /* check the atom index. */
  if (id < 1 || id > P->n_atoms) {
    /* raise an exception and return null. */
    raise("atom index %d out of bounds [1,%u]", id, P->n_atoms);
    return NULL;
  }

  /* allocate a new assignment set. */
  set = assign_none(P);
  if (!set)
    return NULL;

  /* store the single set element. */
  set->v[0] = (unsigned int) id - 1;
  set->n = 1;

  /* return the new set pointer. */
  return set;
}

/* assign_resid(): construct an assignment set for all atoms having a
 * specified residue index.
 *
 * arguments:
 *  @P: pointer to the peptide to reference set operations against.
 *  @id: residue index of atoms to include in the set.
 *
 * returns:
 *  pointer to a newly allocated assignment set, or NULL on failure.
 */
assign_set_t *assign_resid (peptide_t *P, unsigned int id) {
  /* declare required variables:
   *  @set: pointer to the output set.
   *  @i: atom index loop counter.
   */
  assign_set_t *set;
  unsigned int i;

  /* check the residue index. */
  if (id < 1 || id > P->n_res) {
    /* raise an exception and return null. */
    raise("residue index %d out of bounds [1,%u]", id, P->n_res);
    return NULL;
  }

  /* allocate a new assignment set. */
  set = assign_none(P);
  if (!set)
    return NULL;

  /* loop over all atoms in the peptide. */
  for (i = 0; i < set->n_max; i++) {
    /* if the atom has the correct residue index, add it to the set. */
    if (P->atoms[i].res_id == id - 1)
      set->v[set->n++] = i;
  }

  /* return the new set pointer. */
  return set;
}

/* assign_resname(): construct an assignment set for all atoms having a
 * specified residue name (i.e. single- or three-letter code).
 *
 * arguments:
 *  @P: pointer to the peptide to reference set operations against.
 *  @name: residue name of atoms to include in the set.
 *
 * returns:
 *  pointer to a newly allocated assignment set, or NULL on failure.
 */
assign_set_t *assign_resname (peptide_t *P, const char *name) {
  /* declare required variables:
   *  @set: pointer to the output set.
   *  @i: atom index loop counter.
   *  @id: residue lookup index.
   */
  assign_set_t *set;
  unsigned int i, id;

  /* check if the residue name string contains wildcards. */
  if (assign_is_wild(name)) {
    /* raise an exception and return null. */
    raise("illegal wildards in residue name '%s'", name);
    return NULL;
  }

  /* determine the residue table index from its name. */
  if (strlen(name) == 1)
    id = resid_lookup1(name[0]);
  else if (strlen(name) == 3)
    id = resid_lookup3(name);
  else
    id = 0;

  /* check that a valid residue table index was found. */
  if (!id) {
    /* raise an exception and return null. */
    raise("invalid residue name '%s'", name);
    return NULL;
  }

  /* allocate a new assignment set. */
  set = assign_none(P);
  if (!set)
    return NULL;

  /* loop over all atoms in the peptide. */
  for (i = 0; i < set->n_max; i++) {
    /* if the atom has the correct residue name, add it to the set. */
    if (P->res[P->atoms[i].res_id] == id)
      set->v[set->n++] = i;
  }

  /* return the new set pointer. */
  return set;
}

/* assign_name(): construct an assignment set for all atoms having a
 * specified atom name.
 *
 * arguments:
 *  @P: pointer to the peptide to reference set operations against.
 *  @name: atom name of atoms to include in the set.
 *
 * returns:
 *  pointer to a newly allocated assignment set, or NULL on failure.
 */
assign_set_t *assign_name (peptide_t *P, const char *name) {
  /* declare required variables:
   *  @set: pointer to the output set.
   *  @i: atom index loop counter.
   *  @uname: uppercase name string.
   */
  assign_set_t *set;
  unsigned int i;
  char *uname;

  /* check if the atom name string contains wildcards. */
  if (assign_is_wild(name)) {
    /* raise an exception and return null. */
    raise("illegal wildards in atom name '%s'", name);
    return NULL;
  }

  /* allocate a new assignment set. */
  set = assign_none(P);
  if (!set)
    return NULL;

  /* build the uppercase name string. */
  uname = strtoupper(name);

  /* check that the name string was duplicated successfully. */
  if (!uname) {
    /* raise an exception and return null. */
    raise("unable to duplicate atom name '%s'", name);
    assign_free(set);
    return NULL;
  }

  /* loop over all atoms in the peptide. */
  for (i = 0; i < set->n_max; i++) {
    /* if the atom has the correct name, add it to the set. */
    if (strcmp(P->atoms[i].name, uname) == 0)
      set->v[set->n++] = i;
  }

  /* return the new set pointer. */
  free(uname);
  return set;
}

/* assign_type(): construct an assignment set for all atoms having a
 * specified atom type (i.e. H, C, N, O).
 *
 * arguments:
 *  @P: pointer to the peptide to reference set operations against.
 *  @type: type character to use for atom selection. (uppercase only)
 *
 * returns:
 *  pointer to a newly allocated assignment set, or NULL on failure.
 */
assign_set_t *assign_type (peptide_t *P, char type) {
  /* declare required variables:
   *  @set: pointer to the output set.
   *  @i: atom index loop counter.
   */
  assign_set_t *set;
  unsigned int i;

  /* allocate a new assignment set. */
  set = assign_none(P);
  if (!set)
    return NULL;

  /* loop over all atoms in the peptide. */
  for (i = 0; i < set->n_max; i++) {
    /* if the atom has the specified type, add it to the set. */
    if (P->atoms[i].type[0] == type)
      set->v[set->n++] = i;
  }

  /* return the new set pointer. */
  return set;
}

/* assign_not(): construct an assignment set by negating an input set
 * against its complete set.
 *
 * arguments:
 *  @set1: input set to use for the operation.
 *
 * returns:
 *  pointer to a newly allocated assignment set, or NULL on failure.
 */
assign_set_t *assign_not (assign_set_t *set1) {
  /* declare required variables:
   *  @set: pointer to the output set.
   *  @i, @i1: atom index loop counters.
   */
  assign_set_t *set;
  unsigned int i, i1;

  /* check if the input set is null. */
  if (!set1) {
    /* raise an exception and return null. */
    raise("input set to not-op is null");
    return NULL;
  }

  /* allocate a new assignment set. */
  set = assign_none(set1->P);
  if (!set)
    return NULL;

  /* loop over all atoms in the peptide. */
  for (i = i1 = 0; i < set->n_max; i++) {
    /* skip all indices in the input set that have already been visited
     * by the main loop over @i.
     */
    while (set1->n && set1->v[i1] < i && i1 + 1 < set1->n) i1++;

    /* do not add indices that are in the input set. */
    if (set1->n && set1->v[i1] == i)
      continue;

    /* add the index to the output set. */
    set->v[set->n++] = i;
  }

  /* return the new set pointer. */
  return set;
}

/* assign_and(): construct an assignment set by intersecting two input sets
 * against each other.
 *
 * arguments:
 *  @set1: first input set to use for the operation.
 *  @set2: second input set to use for the operation.
 *
 * returns:
 *  pointer to a newly allocated assignment set, or NULL on failure.
 */
assign_set_t *assign_and (assign_set_t *set1, assign_set_t *set2) {
  /* declare required variables:
   *  @set: pointer to the output set.
   *  @i1, @i1: atom index loop counters.
   */
  assign_set_t *set;
  unsigned int i1, i2;

  /* check if the input set is null. */
  if (!set1 || !set2) {
    /* raise an exception and return null. */
    raise("input set to and-op is null");
    return NULL;
  }

  /* allocate a new assignment set. */
  set = assign_none(set1->P);
  if (!set)
    return NULL;

  /* return if either set contains no atoms. */
  if (set1->n == 0 || set2->n == 0)
    return set;

  /* loop over all atoms in the first set. */
  for (i1 = i2 = 0; i1 < set1->n; i1++) {
    /* skip all indices in the second set that have already
     * been visited by the main loop over the first set.
     */
    while (set2->v[i2] < set1->v[i1] && i2 + 1 < set2->n) i2++;

    /* add indices that are in both sets. */
    if (set1->v[i1] == set2->v[i2])
      set->v[set->n++] = set1->v[i1];
  }

  /* return the new set pointer. */
  return set;
}

/* assign_or(): construct an assignment set by unioning two input sets
 * with each other.
 *
 * arguments:
 *  @set1: first input set to use for the operation.
 *  @set2: second input set to use for the operation.
 *
 * returns:
 *  pointer to a newly allocated assignment set, or NULL on failure.
 */
assign_set_t *assign_or (assign_set_t *set1, assign_set_t *set2) {
  /* declare required variables:
   *  @set: pointer to the output set.
   *  @i, @i1, @i1: atom index loop counters.
   */
  assign_set_t *set;
  unsigned int i, i1, i2;

  /* check if the input set is null. */
  if (!set1 || !set2) {
    /* raise an exception and return null. */
    raise("input set to or-op is null");
    return NULL;
  }

  /* allocate a new assignment set. */
  set = assign_none(set1->P);
  if (!set)
    return NULL;

  /* loop over all atoms in the peptide. */
  for (i = i1 = i2 = 0; i < set->n_max; i++) {
    /* skip all indices in the first set that have already been visited
     * by the main loop over @i.
     */
    while (set1->n && set1->v[i1] < i && i1 + 1 < set1->n) i1++;

    /* skip all indices in the second set that have already been visited
     * by the main loop over @i.
     */
    while (set2->n && set2->v[i2] < i && i2 + 1 < set2->n) i2++;

    /* add indices that are in either input set. */
    if ((set1->n && set1->v[i1] == i) ||
        (set2->n && set2->v[i2] == i))
      set->v[set->n++] = i;
  }

  /* return the new set pointer. */
  return set;
}

/* assign_set_distance(): add a distance restraint to a peptide structure.
 *
 * arguments:
 *  @P: pointer to the peptide structure to modify.
 *  @set1: first input set to use for the assignment.
 *  @set2: second input set to use for the assignment.
 *  @d: central distance value of the restraint.
 *  @dmin: lower distance offset of the restraint.
 *  @dplus: upper distance offset of the restraint.
 *
 * returns:
 *  integer indicating whether (1) or not (0) the operation succeeded.
 */
int assign_set_distance (peptide_t *P,
                         assign_set_t *set1,
                         assign_set_t *set2,
                         double d, double dmin, double dplus) {
  /* declare required variables:
   *  @i1, @i2: atom indices to use for the restraint.
   *  @l, @u: restraint lower and upper bounds.
   */
  unsigned int i1, i2;
  double l, u;

  /* check if any assignment set pointer is null. */
  if (!set1 || !set2)
    throw("one or more null-pointer atom selectors");

  /* check if any assignment set is empty. */
  if (assign_is_empty(set1, 1) ||
      assign_is_empty(set2, 2))
    throw("one or more empty atom selectors");

  /* check if any assignment set is ambiguous. */
  if (assign_is_ambiguous(set1, 1) ||
      assign_is_ambiguous(set2, 2))
    throw("one or more ambiguous atom selectors");

  /* extract the atom indices. */
  i1 = set1->v[0];
  i2 = set2->v[0];

  /* compute the restraint bounds. */
  l = d - dmin;
  u = d + dplus;

  /* check the restraint bounds. */
  if (l < 0.0 || u < 0.0 || l > u)
    throw("invalid distance restraint [%.3lf,%.3lf]", l, u);

  /* add the distance restraint to the peptide, in the form of a bond.
   */
  if (!peptide_bond_add(P,
        P->atoms[i1].res_id, P->atoms[i1].name,
        P->atoms[i2].res_id, P->atoms[i2].name, 1))
    throw("unable to add distance restraint");

  /* store the value of the restraint. */
  P->bonds[P->n_bonds - 1].len = value_interval(l, u);

  /* return success. */
  return 1;
}

/* assign_set_dihedral(): add a dihedral angle restraint to a peptide
 * structure.
 *
 * arguments:
 *  @P: pointer to the peptide structure to modify.
 *  @set1: first input set to use for the assignment.
 *  @set2: second input set to use for the assignment.
 *  @set3: third input set to use for the assignment.
 *  @set4: fourth input set to use for the assignment.
 *  @phi: central angle value of the restraint.
 *  @dphi: symmetric angle offset of the restraint.
 *
 * returns:
 *  integer indicating whether (1) or not (0) the operation succeeded.
 */
int assign_set_dihedral (peptide_t *P,
                         assign_set_t *set1,
                         assign_set_t *set2,
                         assign_set_t *set3,
                         assign_set_t *set4,
                         double phi, double dphi) {
  /* declare required variables:
   *  @i1, @i2, @i3, @i4: atom indices to use for the restraint.
   *  @l, @u: restraint lower and upper bounds.
   */
  unsigned int i1, i2, i3, i4;
  double l, u;

  /* check if any assignment set pointer is null. */
  if (!set1 || !set2 || !set3 || !set4)
    throw("one or more null-pointer atom selectors");

  /* check if any assignment set is empty. */
  if (assign_is_empty(set1, 1) ||
      assign_is_empty(set2, 2) ||
      assign_is_empty(set3, 3) ||
      assign_is_empty(set4, 4))
    throw("one or more empty atom selectors");

  /* check if any assignment set is ambiguous. */
  if (assign_is_ambiguous(set1, 1) ||
      assign_is_ambiguous(set2, 2) ||
      assign_is_ambiguous(set3, 3) ||
      assign_is_ambiguous(set4, 4))
    throw("one or more ambiguous atom selectors");

  /* extract the atom indices. */
  i1 = set1->v[0];
  i2 = set2->v[0];
  i3 = set3->v[0];
  i4 = set4->v[0];

  /* compute the restraint bounds. */
  l = phi - dphi / 2.0;
  u = phi + dphi / 2.0;

  /* check the restraint bounds. */
  if (l < -180.0 || u > 180.0 || l > u)
    throw("invalid dihedral restraint [%.3lf,%.3lf]", l, u);

  /* add the dihedral restraint to the peptide, in the form of
   * an improper torsion angle.
   */
  if (!peptide_improper_add(P,
        P->atoms[i1].res_id, P->atoms[i1].name,
        P->atoms[i2].res_id, P->atoms[i2].name,
        P->atoms[i3].res_id, P->atoms[i3].name,
        P->atoms[i4].res_id, P->atoms[i4].name))
    throw("unable to add dihedral restraint");

  /* store the value of the restraint. */
  P->impropers[P->n_impropers - 1].ang = value_interval(l, u);

  /* return success. */
  return 1;
}

/* assign_set_from_file(): parse a restraint assignments file in order to
 * add distance and dihedral restraints to a peptide data structure.
 *
 * arguments:
 *  @P: pointer to the peptide structure to modify.
 *  @fname: input restraints filename string.
 *
 * returns:
 *  integer indicating whether (1) or not (0) the operation succeeded.
 */
int assign_set_from_file (peptide_t *P, const char *fname) {
  /* declare required variables:
   *  @fh: input file handle.
   */
  FILE *fh;

  /* open the input file. */
  fh = fopen(fname, "r");

  /* check that the file was opened. */
  if (!fh)
    throw("unable to open '%s' for reading", fname);

  /* parse the input file. */
  if (!assign_parse(fh, P)) {
    fclose(fh);
    throw("unable to parse restraints from '%s'", fname);
  }

  /* close the input file. */
  fclose(fh);

  /* return success. */
  return 1;
}

