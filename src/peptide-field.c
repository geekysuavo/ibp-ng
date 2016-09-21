
/* include all peptide headers. */
#include "peptide.h"
#include "peptide-atoms.h"
#include "peptide-bonds.h"
#include "peptide-angles.h"
#include "peptide-torsions.h"
#include "peptide-impropers.h"

/* peptide_field(): back-calculate the probability parameters
 * (or equivalently, the force field parameters) from a peptide
 * data structure.
 *
 * arguments:
 *  @P: pointer to the peptide structure to modify.
 *  @tol: distance tolerance to use for exact quantities.
 *
 * returns:
 *  integer indicating whether (1) or not (0) the operation succeeded.
 */
int peptide_field (peptide_t *P, double tol) {
  /* check that the structure pointer is valid. */
  if (!P)
    throw("input structure pointer is null");

  /* back-calculate bond parameters. */
  if (!peptide_field_bonds(P, tol))
    throw("unable to recompute bond parameters");

  /* back-calculate angle parameters. */
  if (!peptide_field_angles(P, tol))
    throw("unable to recompute angle parameters");

  /* back-calculate torsion parameters. */
  if (!peptide_field_torsions(P, tol))
    throw("unable to recompute torsional dihedral parameters");

  /* back-calculate improper parameters. */
  if (!peptide_field_impropers(P, tol))
    throw("unable to recompute improper dihedral parameters");

  /* return success. */
  return 1;
}

