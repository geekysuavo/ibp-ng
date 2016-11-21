
/* include the peptide headers. */
#include "peptide.h"
#include "peptide-atoms.h"

/* peptide_torsion_find(): lookup a torsional dihedral in a peptide
 * structure by its atom array indices.
 *
 * arguments:
 *  @P: pointer to the peptide structure to access.
 *  @id*: array indices of the query atoms.
 *
 * returns:
 *  index of the queried torsion in the peptide, or -1 if no match.
 */
int peptide_torsion_find (peptide_t *P,
                          unsigned int id1, unsigned int id2,
                          unsigned int id3, unsigned int id4) {
  /* declare required variables:
   *  @i: torsion array index.
   *  @ids: atom indices.
   */
  unsigned int i, *ids;

  /* loop over the torsion array. */
  for (i = 0; i < P->n_torsions; i++) {
    /* get the atom indices in the torsion. */
    ids = P->torsions[i].atom_id;

    /* return true if the current torsion is a match. */
    if (ids[0] == id1 && ids[1] == id2 &&
        ids[2] == id3 && ids[3] == id4)
      return (int) i;

    /* return true if the current torsion is a reverse match. */
    if (ids[0] == id4 && ids[1] == id3 &&
        ids[2] == id2 && ids[3] == id1)
      return (int) i;
  }

  /* return false. */
  return -1;
}

/* peptide_torsion_add(): add a torsional dihedral angle between four atoms
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
int peptide_torsion_add (peptide_t *P,
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

  /* fail if the peptide already contains the torsion. */
  if (peptide_torsion_find(P, ia, ib, ic, id) >= 0)
    throw("torsion (%u.%s, %u.%s, %u.%s, %u.%s) already exists",
          resid1 + 1, name1, resid2 + 1, name2,
          resid3 + 1, name3, resid4 + 1, name4);

  /* increment the torsion array length. */
  i = P->n_torsions;
  P->n_torsions++;

  /* reallocate the torsion array. */
  P->torsions = (peptide_dihed_t*)
    realloc(P->torsions, P->n_torsions * sizeof(peptide_dihed_t));

  /* check if reallocation failed. */
  if (!P->torsions)
    throw("unable to reallocate torsion array");

  /* store the atom indices. */
  P->torsions[i].atom_id[0] = ia;
  P->torsions[i].atom_id[1] = ib;
  P->torsions[i].atom_id[2] = ic;
  P->torsions[i].atom_id[3] = id;

  /* initialize the angle. */
  P->torsions[i].ang = value_undefined();

  /* initialize the probability parameters. */
  P->torsions[i].mu = 0.0;
  P->torsions[i].kappa = 0.0;

  /* return success. */
  return 1;
}

/* peptide_torsion_delete(): delete an existing torsional dihedral angle
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
int peptide_torsion_delete (peptide_t *P,
                            unsigned int resid1, const char *name1,
                            unsigned int resid2, const char *name2,
                            unsigned int resid3, const char *name3,
                            unsigned int resid4, const char *name4) {
  /* declare required variables:
   *  @ia, @ib, @ic, @id: atom indices.
   *  @i: torsion index.
   */
  int i, ia, ib, ic, id;

  /* locate the index of the first atom in the torsion. */
  ia = peptide_atom_find(P, resid1, name1);
  if (ia < 0)
    return 1;

  /* locate the index of the second atom in the torsion. */
  ib = peptide_atom_find(P, resid2, name2);
  if (ib < 0)
    return 1;

  /* locate the index of the third atom in the torsion. */
  ic = peptide_atom_find(P, resid3, name3);
  if (ic < 0)
    return 1;

  /* locate the index of the fourth atom in the torsion. */
  id = peptide_atom_find(P, resid4, name4);
  if (id < 0)
    return 1;

  /* return if the peptide does not contain the torsion. */
  i = peptide_torsion_find(P, ia, ib, ic, id);
  if (i < 0)
    return 1;

  /* swap the atom indices within the torsion. */
  if ((unsigned int) i != P->n_torsions - 1)
    memcpy(P->torsions + i,
           P->torsions + (P->n_torsions - 1),
           sizeof(peptide_dihed_t));

  /* decrement the torsion count and return success. */
  P->n_torsions--;
  return 1;
}

/* peptide_torsion_delete_any(): delete any torsional dihedral angles from
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
int peptide_torsion_delete_any (peptide_t *P,
                                unsigned int resid,
                                const char *name) {
  /* declare required variables:
   *  @i: torsion index.
   *  @ids: atom indices.
   */
  unsigned int i, *ids;

  /* loop over the torsions in the peptide. */
  for (i = 0; P->n_torsions && i < P->n_torsions; i++) {
    /* get the atom indices. */
    ids = P->torsions[i].atom_id;

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

    /* swap the atom indices within the last torsion. */
    if (i != P->n_torsions - 1)
      memcpy(P->torsions + i,
             P->torsions + (P->n_torsions - 1),
             sizeof(peptide_dihed_t));

    /* decrement the torsion count. */
    P->n_torsions--;
  }

  /* return success. */
  return 1;
}

/* peptide_field_torsions(): compute the probability parameters for
 * each torsion in a peptide structure.
 *
 * FIXME: currently, an approximation is used that assumes relatively
 * high precision. of course, this will break down when [-180,180] is
 * used to generate a completely uninformative interval.
 *
 * one hacky solution would be to place [-x,x] into the ibp-protein.par
 * file, where x >> 180, to force the resulting von mises distribution
 * be closer to a circular uniform distribution.
 *
 * arguments:
 *  @P: pointer to the peptide structure to modify.
 *  @tol: angular tolerance to add into intervals.
 *
 * returns:
 *  integer indicating whether (1) or not (0) the operation succeeded.
 */
int peptide_field_torsions (peptide_t *P, double tol) {
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
  for (i = 0; i < P->n_torsions; i++) {
    /* get the current dihedral angle. */
    L = P->torsions[i].ang.l * (M_PI / 180.0) - tol;
    U = P->torsions[i].ang.u * (M_PI / 180.0) + tol;

    /* rectify the bounds. */
    if (L < -M_PI) L = -M_PI;
    if (U >  M_PI) U =  M_PI;

    /* compute the parameters. */
    kappa = Z * pow(U - L, -2.0);
    mu = L - zl / sqrt(kappa);

    /* store the computed parameters. */
    P->torsions[i].mu = mu;
    P->torsions[i].kappa = kappa;
  }

  /* return success. */
  return 1;
}

/* peptide_graph_torsions(): update the edge set of a graph structure
 * using peptide torsional dihedral angle information.
 *
 * arguments:
 *  @P: pointer to the peptide structure to access.
 *  @G: pointer to the graph structure to modify.
 *
 * returns:
 *  integer indicating whether (1) or not (0) the operation succeeded.
 */
int peptide_graph_torsions (peptide_t *P, graph_t *G) {
  /* declare required variables:
   *  @omega: angle value.
   *  @d01, @d02, @d12, @d13, @d23, @d03: distance values.
   *  @i: peptide torsion array index.
   *  @atoms: atom array indices.
   */
  value_t omega, d01, d02, d12, d13, d23, d03;
  unsigned int i, *atoms;

  /* loop over the torsions in the peptide. */
  for (i = 0; i < P->n_torsions; i++) {
    /* get the current torsion information. */
    atoms = P->torsions[i].atom_id;

    /* scale and bound the torsion. */
    omega = value_scal(P->torsions[i].ang, M_PI / 180.0);
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
        value_is_undefined(d23))
      throw("undefined distance in dihedral (%s%u) %s-%s-%s-%s",
            resid_get_code3(P->res[P->atoms[atoms[0]].res_id]),
            P->atoms[atoms[0]].res_id + 1,
            P->atoms[atoms[0]].name,
            P->atoms[atoms[1]].name,
            P->atoms[atoms[2]].name,
            P->atoms[atoms[3]].name);

    /* compute the edge weight from the torsion parameters. */
    d03 = value_from_dihedral(d01, d02, d12, d13, d23, omega);

    /* attempt to refine the graph edge associated with the torsion. */
    if (!graph_refine_edge(G, atoms[0], atoms[3], d03))
      throw("unable to refine graph edge from %u.%s to %u.%s",
            P->atoms[atoms[0]].res_id + 1, P->atoms[atoms[0]].name,
            P->atoms[atoms[3]].res_id + 1, P->atoms[atoms[3]].name);
  }

  /* return success. */
  return 1;
}

