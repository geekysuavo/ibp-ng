
/* include the peptide headers. */
#include "peptide.h"
#include "peptide-atoms.h"

/* peptide_angle_find(): lookup an angle in a peptide structure by its
 * atom array indices.
 *
 * arguments:
 *  @P: pointer to the peptide structure to access.
 *  @id*: array indices of the query atoms.
 *
 * returns:
 *  index of the queried angle in the peptide, or -1 if no match.
 */
int peptide_angle_find (peptide_t *P,
                        unsigned int id1,
                        unsigned int id2,
                        unsigned int id3) {
  /* declare required variables:
   *  @i: angle array index.
   *  @ids: atom indices.
   */
  unsigned int i, *ids;

  /* loop over the angle array. */
  for (i = 0; i < P->n_angles; i++) {
    /* get the atom indices in the bond. */
    ids = P->angles[i].atom_id;

    /* return true if the current angle is a match. */
    if (ids[1] == id2 &&
        ((ids[0] == id1 && ids[2] == id3) ||
         (ids[0] == id3 && ids[2] == id1)))
      return (int) i;
  }

  /* return false. */
  return -1;
}

/* peptide_angle_add(): add an angle between three atoms in
 * a peptide structure.
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
int peptide_angle_add (peptide_t *P,
                       unsigned int resid1, const char *name1,
                       unsigned int resid2, const char *name2,
                       unsigned int resid3, const char *name3) {
  /* declare required variables:
   *  @ia, @ib, @ic: atom indices.
   *  @i: angle index.
   */
  unsigned int i;
  int ia, ib, ic;

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

  /* fail if the peptide already contains the angle. */
  if (peptide_angle_find(P, ia, ib, ic) >= 0)
    throw("angle (%u.%s, %u.%s, %u.%s) already exists",
          resid1 + 1, name1,
          resid2 + 1, name2,
          resid3 + 1, name3);

  /* increment the angle array length. */
  i = P->n_angles;
  P->n_angles++;

  /* reallocate the angle array. */
  P->angles = (peptide_angle_t*)
    realloc(P->angles, P->n_angles * sizeof(peptide_angle_t));

  /* check if reallocation failed. */
  if (!P->angles)
    throw("unable to reallocate angle array");

  /* store the atom indices. */
  P->angles[i].atom_id[0] = ia;
  P->angles[i].atom_id[1] = ib;
  P->angles[i].atom_id[2] = ic;

  /* initialize the angle. */
  P->angles[i].ang = value_undefined();

  /* initialize the probability parameters. */
  P->angles[i].mu = 0.0;
  P->angles[i].kappa = 0.0;

  /* return success. */
  return 1;
}

/* peptide_angle_delete(): delete an existing angle from a peptide structure.
 *
 * arguments:
 *  @P: pointer to the peptide structure to modify.
 *  @resid*: residue sequence indices of the atoms.
 *  @name*: names of the atoms.
 *
 * returns:
 *  integer indicating whether (1) or not (0) the operation succeeded.
 */
int peptide_angle_delete (peptide_t *P,
                          unsigned int resid1, const char *name1,
                          unsigned int resid2, const char *name2,
                          unsigned int resid3, const char *name3) {
  /* declare required variables:
   *  @ia, @ib, @ic: atom indices.
   *  @i: angle index.
   */
  int i, ia, ib, ic;

  /* locate the index of the first atom in the angle. */
  ia = peptide_atom_find(P, resid1, name1);
  if (ia < 0)
    return 1;

  /* locate the index of the second atom in the angle. */
  ib = peptide_atom_find(P, resid2, name2);
  if (ib < 0)
    return 1;

  /* locate the index of the third atom in the angle. */
  ic = peptide_atom_find(P, resid3, name3);
  if (ic < 0)
    return 1;

  /* return if the peptide does not contain the angle. */
  i = peptide_angle_find(P, ia, ib, ic);
  if (i < 0)
    return 1;

  /* swap the atom indices with the last angle. */
  if ((unsigned int) i != P->n_angles - 1)
    memcpy(P->angles + i,
           P->angles + (P->n_angles - 1),
           sizeof(peptide_angle_t));

  /* decrement the angle count return success. */
  P->n_angles--;
  return 1;
}

/* peptide_angle_delete_any(): delete any angles from a peptide structure
 * that contain a specified atom.
 *
 * arguments:
 *  @P: pointer to the peptide structure to modify.
 *  @resid: residue sequence index of the atom.
 *  @name: name of the atom.
 *
 * returns:
 *  integer indicating whether (1) or not (0) the operation succeeded.
 */
int peptide_angle_delete_any (peptide_t *P,
                              unsigned int resid,
                              const char *name) {
  /* declare required variables:
   *  @i: angle index.
   *  @ids: atom indices.
   */
  unsigned int i, *ids;

  /* loop over the angles in the peptide. */
  for (i = 0; P->n_angles && i < P->n_angles; i++) {
    /* get the atom indices. */
    ids = P->angles[i].atom_id;

    /* skip if no atom is a match. */
    if (!(strcmp(P->atoms[ids[0]].name, name) == 0 &&
          P->atoms[ids[0]].res_id == resid) &&
        !(strcmp(P->atoms[ids[1]].name, name) == 0 &&
          P->atoms[ids[1]].res_id == resid) &&
        !(strcmp(P->atoms[ids[2]].name, name) == 0 &&
          P->atoms[ids[2]].res_id == resid))
      continue;

    /* swap the atom indices with the last angle. */
    if (i != P->n_angles - 1)
      memcpy(P->angles + i,
             P->angles + (P->n_angles - 1),
             sizeof(peptide_angle_t));

    /* decrement the angle count. */
    P->n_angles--;
  }

  /* return success. */
  return 1;
}

/* peptide_field_angles(): compute the probability parameters for
 * each angle in a peptide structure.
 *
 * arguments:
 *  @P: pointer to the peptide structure to modify.
 *  @tol: angular tolerance to add into intervals.
 *
 * returns:
 *  integer indicating whether (1) or not (0) the operation succeeded.
 */
int peptide_field_angles (peptide_t *P, double tol) {
  /* declare required variables:
   *  @i: peptide angle array index.
   *  @L, @U: angle value interval.
   */
  double L, U, mu, kappa;
  unsigned int i;

  /* define the 0.025-0.975 quantile range for the standard normal.
   * this is an approximation, since angles are von mises, but
   * since our angles will be precise (exact), it works.
   */
  const double zl = -1.95996938;
  const double zu =  1.95996938;
  const double Z = pow(zu - zl, 2.0);

  /* loop over the angles in the peptide. */
  for (i = 0; i < P->n_angles; i++) {
    /* get the current bond angle. */
    L = P->angles[i].ang.l * (M_PI / 180.0) - tol;
    U = P->angles[i].ang.u * (M_PI / 180.0) + tol;

    /* rectify the bounds. */
    if (L <  0.0) L =  0.0;
    if (U > M_PI) U = M_PI;

    /* compute the parameters. */
    kappa = Z * pow(U - L, -2.0);
    mu = L - zl / sqrt(kappa);

    /* store the computed parameters. */
    P->angles[i].mu = mu;
    P->angles[i].kappa = kappa;
  }

  /* return success. */
  return 1;
}

/* peptide_graph_angles(): update the edge set of a graph structure
 * using peptide angle information.
 *
 * arguments:
 *  @P: pointer to the peptide structure to access.
 *  @G: pointer to the graph structure to modify.
 *
 * returns:
 *  integer indicating whether (1) or not (0) the operation succeeded.
 */
int peptide_graph_angles (peptide_t *P, graph_t *G) {
  /* declare required variables:
   *  @theta: angle value.
   *  @d01, @d12, @d02: distance values.
   *  @i: peptide angle array index.
   *  @atoms: atom array indices.
   */
  value_t theta, d01, d12, d02;
  unsigned int i, *atoms;

  /* loop over the angles in the peptide. */
  for (i = 0; i < P->n_angles; i++) {
    /* get the current angle information. */
    atoms = P->angles[i].atom_id;

    /* scale and bound the angle. */
    theta = value_scal(P->angles[i].ang, M_PI / 180.0);
    theta = value_bound(theta, value_interval(0.0, M_PI));

    /* get the known distances between each pair of atoms. */
    d01 = graph_get_edge(G, atoms[0], atoms[1]);
    d12 = graph_get_edge(G, atoms[1], atoms[2]);

    /* check that the required distances are defined. */
    if (value_is_undefined(d01) ||
        value_is_undefined(d12))
      throw("undefined distance in angle (%s%u) %s-%s-%s",
            peptide_get_resname(P, P->atoms[atoms[0]].res_id),
            P->atoms[atoms[0]].res_id + 1,
            P->atoms[atoms[0]].name,
            P->atoms[atoms[1]].name,
            P->atoms[atoms[2]].name);

    /* compute the edge weight from the angle parameters. */
    d02 = value_from_angle(d01, d12, theta);

    /* attempt to refine the graph edge associated with the angle. */
    if (!graph_refine_edge(G, atoms[0], atoms[2], d02,
                           &P->angles[i].ang,
                           VALUE_IS_ANGLE))
      throw("unable to refine graph edge from %u.%s to %u.%s",
            P->atoms[atoms[0]].res_id + 1, P->atoms[atoms[0]].name,
            P->atoms[atoms[2]].res_id + 1, P->atoms[atoms[2]].name);
  }

  /* return success. */
  return 1;
}

