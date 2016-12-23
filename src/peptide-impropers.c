
/* include the peptide headers. */
#include "peptide.h"
#include "peptide-atoms.h"

/* peptide_improper_find(): lookup an improper dihedral in a peptide
 * structure by its atom array indices.
 *
 * arguments:
 *  @P: pointer to the peptide structure to access.
 *  @id*: array indices of the query atoms.
 *
 * returns:
 *  index of the queried improper in the peptide, or -1 if no match.
 */
int peptide_improper_find (peptide_t *P,
                           unsigned int id1, unsigned int id2,
                           unsigned int id3, unsigned int id4) {
  /* declare required variables:
   *  @i: improper array index.
   *  @ids: atom indices.
   */
  unsigned int i, *ids;

  /* loop over the improper array. */
  for (i = 0; i < P->n_impropers; i++) {
    /* get the atom indices in the impropers. */
    ids = P->impropers[i].atom_id;

    /* return true if the current improper is a match. */
    if (ids[0] == id1 && ids[1] == id2 &&
        ids[2] == id3 && ids[3] == id4)
      return (int) i;

    /* return true if the current improper is a reverse match. */
    if (ids[0] == id4 && ids[1] == id3 &&
        ids[2] == id2 && ids[3] == id1)
      return (int) i;
  }

  /* return false. */
  return -1;
}

/* peptide_improper_add(): add an improper dihedral angle between four atoms
 * in a peptide structure.
 *
 * arguments:
 *  @P: pointer to the peptide structure to modify.
 *  @resid*: residue sequence indices of the atoms.
 *  @name*: names of the atoms.
 *
 * returns:
 *  integer indicating whether (1) or not (0,-1) the angle was successfully
 *  added to the peptide.
 */
int peptide_improper_add (peptide_t *P,
                          unsigned int resid1, const char *name1,
                          unsigned int resid2, const char *name2,
                          unsigned int resid3, const char *name3,
                          unsigned int resid4, const char *name4) {
  /* declare required variables:
   *  @ia, @ib, @ic, @id: atom indices.
   *  @i: angle index.
   */
  int ia, ib, ic, id;
  unsigned int i;

  /* locate the index of the first atom in the angle. */
  ia = peptide_atom_find(P, resid1, name1);
  if (ia < 0)
    throw("atom %u.%s does not exist", resid1, name1);

  /* locate the index of the second atom in the angle. */
  ib = peptide_atom_find(P, resid2, name2);
  if (ib < 0)
    throw("atom %u.%s does not exist", resid2, name2);

  /* locate the index of the third atom in the angle. */
  ic = peptide_atom_find(P, resid3, name3);
  if (ic < 0)
    throw("atom %u.%s does not exist", resid3, name3);

  /* locate the index of the fourth atom in the angle. */
  id = peptide_atom_find(P, resid4, name4);
  if (id < 0)
    throw("atom %u.%s does not exist", resid4, name4);

  /* fail if the peptide already contains the improper. */
  if (peptide_improper_find(P, ia, ib, ic, id) >= 0)
    throw("improper (%u.%s, %u.%s, %u.%s, %u.%s) already exists",
          resid1 + 1, name1, resid2 + 1, name2,
          resid3 + 1, name3, resid4 + 1, name4);

  /* increment the improper array length. */
  i = P->n_impropers;
  P->n_impropers++;

  /* reallocate the improper array. */
  P->impropers = (peptide_dihed_t*)
    realloc(P->impropers, P->n_impropers * sizeof(peptide_dihed_t));

  /* check if reallocation failed. */
  if (!P->impropers)
    throw("unable to reallocate improper array");

  /* store the atom indices. */
  P->impropers[i].atom_id[0] = ia;
  P->impropers[i].atom_id[1] = ib;
  P->impropers[i].atom_id[2] = ic;
  P->impropers[i].atom_id[3] = id;

  /* initialize the angle. */
  P->impropers[i].ang = value_undefined();

  /* initialize the probability parameters. */
  P->impropers[i].mu = 0.0;
  P->impropers[i].kappa = 0.0;

  /* return success. */
  return 1;
}

/* peptide_improper_delete(): delete an existing improper dihedral angle
 * from a peptide structure.
 *
 * arguments:
 *  @P: pointer to the peptide structure to modify.
 *  @resid*: residue sequence indices of the atoms.
 *  @name*: names of the atoms.
 *
 * returns:
 *  integer indicating whether (1) or not (0) the operation succeeded.
 */
int peptide_improper_delete (peptide_t *P,
                             unsigned int resid1, const char *name1,
                             unsigned int resid2, const char *name2,
                             unsigned int resid3, const char *name3,
                             unsigned int resid4, const char *name4) {
  /* declare required variables:
   *  @ia, @ib, @ic, @id: atom indices.
   *  @i: improper index.
   */
  int i, ia, ib, ic, id;

  /* locate the index of the first atom in the improper. */
  ia = peptide_atom_find(P, resid1, name1);
  if (ia < 0)
    return 1;

  /* locate the index of the second atom in the improper. */
  ib = peptide_atom_find(P, resid2, name2);
  if (ib < 0)
    return 1;

  /* locate the index of the third atom in the improper. */
  ic = peptide_atom_find(P, resid3, name3);
  if (ic < 0)
    return 1;

  /* locate the index of the fourth atom in the improper. */
  id = peptide_atom_find(P, resid4, name4);
  if (id < 0)
    return 1;

  /* return if the peptide does not contain the improper. */
  i = peptide_improper_find(P, ia, ib, ic, id);
  if (i < 0)
    return 1;

  /* swap the atom indices within the improper. */
  if ((unsigned int) i != P->n_impropers - 1)
    memcpy(P->impropers + i,
           P->impropers + (P->n_impropers - 1),
           sizeof(peptide_dihed_t));

  /* decrement the improper count and return success. */
  P->n_impropers--;
  return 1;
}

/* peptide_improper_delete_any(): delete any improper dihedral angles from
 * a peptide structure that contain a specified atom.
 *
 * arguments:
 *  @P: pointer to the peptide structure to modify.
 *  @resid: residue sequence index of the atom.
 *  @name: name of the atom.
 *
 * returns:
 *  integer indicating whether (1) or not (0) the operation succeeded.
 */
int peptide_improper_delete_any (peptide_t *P,
                                 unsigned int resid,
                                 const char *name) {
  /* declare required variables:
   *  @i: improper index.
   *  @ids: atom indices.
   */
  unsigned int i, *ids;

  /* loop over the impropers in the peptide. */
  for (i = 0; P->n_impropers && i < P->n_torsions; i++) {
    /* get the atom indices. */
    ids = P->impropers[i].atom_id;

    /* skip if no atom is a match. */
    if (!(strcmp(P->atoms[ids[0]].name, name) == 0 &&
          P->atoms[ids[0]].res_id == resid) &&
        !(strcmp(P->atoms[ids[1]].name, name) == 0 &&
          P->atoms[ids[1]].res_id == resid) &&
        !(strcmp(P->atoms[ids[2]].name, name) == 0 &&
          P->atoms[ids[2]].res_id == resid) &&
        !(strcmp(P->atoms[ids[3]].name, name) == 0 &&
          P->atoms[ids[3]].res_id == resid))
      continue;

    /* swap the atom indices within the last improper. */
    if (i != P->n_impropers - 1)
      memcpy(P->impropers + i,
             P->impropers + (P->n_impropers - 1),
             sizeof(peptide_dihed_t));

    /* decrement the improper count. */
    P->n_impropers--;
  }

  /* return success. */
  return 1;
}

/* peptide_field_impropers(): compute the probability parameters for
 * each improper in a peptide structure.
 *
 * arguments:
 *  @P: pointer to the peptide structure to modify.
 *  @tol: angular tolerance to add into intervals.
 *
 * returns:
 *  integer indicating whether (1) or not (0) the operation succeeded.
 */
int peptide_field_impropers (peptide_t *P, double tol) {
  /* declare required variables:
   *  @i: peptide dihedral array index.
   *  @L, @U: dihedral value interval.
   */
  double L, U, mu, kappa;
  unsigned int i;

  /* define the 0.025-0.975 quantile range for the standard normal. */
  const double zl = -1.95996938;
  const double zu =  1.95996938;
  const double Z = pow(zu - zl, 2.0);

  /* loop over the dihedrals in the peptide. */
  for (i = 0; i < P->n_impropers; i++) {
    /* get the current dihedral angle. */
    L = P->impropers[i].ang.l * (M_PI / 180.0) - tol;
    U = P->impropers[i].ang.u * (M_PI / 180.0) + tol;

    /* rectify the bounds. */
    if (L < -M_PI) L = -M_PI;
    if (U >  M_PI) U =  M_PI;

    /* compute the parameters. */
    kappa = Z * pow(U - L, -2.0);
    mu = L - zl / sqrt(kappa);

    /* store the computed parameters. */
    P->impropers[i].mu = mu;
    P->impropers[i].kappa = kappa;
  }

  /* return success. */
  return 1;
}

/* peptide_graph_impropers(): update the edge set of a graph structure
 * using peptide improper dihedral angle information.
 *
 * arguments:
 *  @P: pointer to the peptide structure to access.
 *  @G: pointer to the graph structure to modify.
 *
 * returns:
 *  integer indicating whether (1) or not (0) the operation succeeded.
 */
int peptide_graph_impropers (peptide_t *P, graph_t *G) {
  /* declare required variables:
   *  @omega: angle value.
   *  @d01, @d02, @d12, @d13, @d23, @d03: distance values.
   *  @i: peptide improper array index.
   *  @atoms: atom array indices.
   */
  value_t omega, d01, d02, d12, d13, d23, d03;
  unsigned int i, *atoms;

  /* loop over the impropers in the peptide. */
  for (i = 0; i < P->n_impropers; i++) {
    /* get the current improper information. */
    atoms = P->impropers[i].atom_id;

    /* scale and bound the improper. */
    omega = value_scal(P->impropers[i].ang, M_PI / 180.0);
    omega = value_bound(omega, value_interval(-M_PI, M_PI));

    /* get the known distances between each pair of atoms. */
    d01 = graph_get_edge(G, atoms[0], atoms[1]);
    d02 = graph_get_edge(G, atoms[0], atoms[2]);
    d12 = graph_get_edge(G, atoms[1], atoms[2]);
    d13 = graph_get_edge(G, atoms[1], atoms[3]);
    d23 = graph_get_edge(G, atoms[2], atoms[3]);

    /* check that the required distances are defined. */
    if (value_is_undefined(d01) ||
        value_is_undefined(d02) ||
        value_is_undefined(d12) ||
        value_is_undefined(d13) ||
        value_is_undefined(d23)) {
      /* output a warning message... */
      warn("undefined distance in improper (%s%u) %s-%s-%s-%s",
           resid_get_code3(P->res[P->atoms[atoms[0]].res_id]),
           P->atoms[atoms[0]].res_id + 1,
           P->atoms[atoms[0]].name,
           P->atoms[atoms[1]].name,
           P->atoms[atoms[2]].name,
           P->atoms[atoms[3]].name);

      /* ... and skip to the next improper. */
      continue;
    }

    /* compute the edge weight from the torsion parameters. */
    d03 = value_from_dihedral(d01, d02, d12, d13, d23, omega);

    /* attempt to refine the graph edge associated with the torsion. */
    if (!graph_refine_edge(G, atoms[0], atoms[3], d03,
                           &P->impropers[i].ang,
                           VALUE_IS_DIHEDRAL))
      throw("unable to refine graph edge from %u.%s to %u.%s",
            P->atoms[atoms[0]].res_id + 1, P->atoms[atoms[0]].name,
            P->atoms[atoms[3]].res_id + 1, P->atoms[atoms[3]].name);
  }

  /* return success. */
  return 1;
}

