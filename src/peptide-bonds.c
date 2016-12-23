
/* include the peptide headers. */
#include "peptide.h"
#include "peptide-atoms.h"

/* peptide_bond_find(): lookup a bond in a peptide structure by its atom
 * array indices.
 *
 * arguments:
 *  @P: pointer to the peptide structure to access.
 *  @id*: array indices of the query atoms.
 *
 * returns:
 *  index of the queried bond in the peptide, or -1 if no match.
 */
int peptide_bond_find (peptide_t *P, unsigned int id1, unsigned int id2) {
  /* declare required variables:
   *  @i: bond array index.
   *  @ids: atom indices.
   */
  unsigned int i, *ids;

  /* loop over the bond array. */
  for (i = 0; i < P->n_bonds; i++) {
    /* get the atom indices in the bond. */
    ids = P->bonds[i].atom_id;

    /* return true if the current bond is a match. */
    if ((ids[0] == id1 && ids[1] == id2) ||
        (ids[0] == id2 && ids[1] == id1))
      return (int) i;
  }

  /* return false. */
  return -1;
}

/* peptide_bond_add(): add a bond between two atoms in a peptide structure.
 *
 * arguments:
 *  @P: pointer to the peptide structure to modify.
 *  @resid1: residue sequence index of the first atom.
 *  @resid2: residue sequence index of the second atom.
 *  @name1: name of the first atom.
 *  @name2: name of the second atom.
 *  @is_virtual: virtual state of the atom. (cf. peptide.h)
 *
 * returns:
 *  integer indicating whether (1) or not (0,-1) the bond was successfully
 *  added to the peptide.
 */
int peptide_bond_add (peptide_t *P,
                      unsigned int resid1, const char *name1,
                      unsigned int resid2, const char *name2,
                      unsigned int is_virtual) {
  /* declare required variables:
   *  @ia, @ib: atom indices.
   *  @i: bond index.
   */
  unsigned int i;
  int ia, ib;

  /* locate the index of the first atom in the bond. */
  ia = peptide_atom_find(P, resid1, name1);
  if (ia < 0)
    return -1;

  /* locate the index of the second atom in the bond. */
  ib = peptide_atom_find(P, resid2, name2);
  if (ib < 0)
    return -1;

  /* fail if the peptide already contains the bond. */
  if (peptide_bond_find(P, ia, ib) >= 0)
    throw("bond (%u.%s, %u.%s) already exists",
          resid1 + 1, name1, resid2 + 1, name2);

  /* increment the bond array length. */
  i = P->n_bonds;
  P->n_bonds++;

  /* reallocate the bond array. */
  P->bonds = (peptide_bond_t*)
    realloc(P->bonds, P->n_bonds * sizeof(peptide_bond_t));

  /* check if reallocation failed. */
  if (!P->bonds)
    throw("unable to reallocate bond array");

  /* store the atom indices. */
  P->bonds[i].atom_id[0] = ia;
  P->bonds[i].atom_id[1] = ib;

  /* initialize the bond length. */
  P->bonds[i].len = value_undefined();

  /* store the virtual state. */
  P->bonds[i].is_virtual = is_virtual;

  /* initialize the probability parameters. */
  P->bonds[i].mu = 0.0;
  P->bonds[i].kappa = 0.0;

  /* return success. */
  return 1;
}

/* peptide_bond_delete(): delete an existing bond from a peptide structure.
 *
 * arguments:
 *  @P: pointer to the peptide structure to modify.
 *  @resid1: residue sequence index of the first atom.
 *  @resid2: residue sequence index of the second atom.
 *  @name1: name of the first atom.
 *  @name2: name of the second atom.
 *
 * returns:
 *  integer indicating whether (1) or not (0) the operation succeeded.
 */
int peptide_bond_delete (peptide_t *P,
                         unsigned int resid1, const char *name1,
                         unsigned int resid2, const char *name2) {
  /* declare required variables:
   *  @ia, @ib: atom indices.
   *  @i: bond index.
   */
  int i, ia, ib;

  /* locate the index of the first atom in the bond. */
  ia = peptide_atom_find(P, resid1, name1);
  if (ia < 0)
    return 1;

  /* locate the index of the second atom in the bond. */
  ib = peptide_atom_find(P, resid2, name2);
  if (ib < 0)
    return 1;

  /* return if the peptide does not contain the bond. */
  i = peptide_bond_find(P, ia, ib);
  if (i < 0)
    return 1;

  /* swap the atom indices with the last bond. */
  if ((unsigned int) i != P->n_bonds - 1)
    memcpy(P->bonds + i,
           P->bonds + (P->n_bonds - 1),
           sizeof(peptide_bond_t));

  /* decrement the bond count and return success. */
  P->n_bonds--;
  return 1;
}

/* peptide_bond_delete_any(): delete any bonds from a peptide structure
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
int peptide_bond_delete_any (peptide_t *P,
                             unsigned int resid,
                             const char *name) {
  /* declare required variables:
   *  @i: bond index.
   *  @ids: atom indices.
   */
  unsigned int i, *ids;

  /* loop over the bonds in the peptide. */
  for (i = 0; P->n_bonds && i < P->n_bonds; i++) {
    /* get the atom indices. */
    ids = P->bonds[i].atom_id;

    /* skip if neither atom is a match. */
    if (!(strcmp(P->atoms[ids[0]].name, name) == 0 &&
          P->atoms[ids[0]].res_id == resid) &&
        !(strcmp(P->atoms[ids[1]].name, name) == 0 &&
          P->atoms[ids[1]].res_id == resid))
      continue;

    /* swap the atom indices with the last bond. */
    if (i != P->n_bonds - 1)
      memcpy(P->bonds + i,
             P->bonds + (P->n_bonds - 1),
             sizeof(peptide_bond_t));

    /* decrement the bond count. */
    P->n_bonds--;
  }

  /* return success. */
  return 1;
}

/* peptide_field_bonds(): compute the probability parameters for
 * each bond in a peptide structure.
 *
 * arguments:
 *  @P: pointer to the peptide structure to modify.
 *  @tol: distance tolerance to add into intervals.
 *
 * returns:
 *  integer indicating whether (1) or not (0) the operation succeeded.
 */
int peptide_field_bonds (peptide_t *P, double tol) {
  /* declare required variables:
   *  @i: peptide bond array index.
   *  @L, @U: bond length interval.
   */
  double L, U, mu, kappa;
  unsigned int i;

  /* define the 0.025-0.975 quantile range for the standard normal. */
  const double zl = -1.95996938;
  const double zu =  1.95996938;
  const double Z = pow(zu - zl, 2.0);

  /* loop over the bonds in the peptide. */
  for (i = 0; i < P->n_bonds; i++) {
    /* get the current bond length. */
    L = P->bonds[i].len.l - tol;
    U = P->bonds[i].len.u + tol;

    /* rectify the bounds. */
    if (L <= 0.0) L = tol;
    if (U <= L) U = L + 2.0 * tol;

    /* compute the parameters. */
    if (P->bonds[i].is_virtual) {
      /* log-normal. */
      kappa = Z * pow(log(U / L), -2.0);
      mu = L * exp(-zl / sqrt(kappa));
    }
    else {
      /* normal. */
      kappa = Z * pow(U - L, -2.0);
      mu = L - zl / sqrt(kappa);
    }

    /* store the computed parameters. */
    P->bonds[i].mu = mu;
    P->bonds[i].kappa = kappa;
  }

  /* return success. */
  return 1;
}

/* peptide_graph_bonds(): update the edge set of a graph structure
 * using peptide bond information.
 *
 * arguments:
 *  @P: pointer to the peptide structure to access.
 *  @G: pointer to the graph structure to modify.
 *
 * returns:
 *  integer indicating whether (1) or not (0) the operation succeeded.
 */
int peptide_graph_bonds (peptide_t *P, graph_t *G) {
  /* declare required variables:
   *  @i: peptide bond array index.
   *  @atoms: atom array indices.
   *  @len: bond length value.
   */
  unsigned int i, *atoms;
  value_t len;

  /* loop over the bonds in the peptide. */
  for (i = 0; i < P->n_bonds; i++) {
    /* get the current bond information. */
    atoms = P->bonds[i].atom_id;
    len = P->bonds[i].len;

    /* attempt to refine the graph edge associated with the bond. */
    if (!graph_refine_edge(G, atoms[0], atoms[1], len,
                           &P->bonds[i].len,
                           VALUE_IS_DISTANCE))
      throw("unable to refine graph edge from %u.%s to %u.%s",
            P->atoms[atoms[0]].res_id + 1, P->atoms[atoms[0]].name,
            P->atoms[atoms[1]].res_id + 1, P->atoms[atoms[1]].name);
  }

  /* return success. */
  return 1;
}

