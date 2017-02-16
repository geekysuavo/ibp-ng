
/* include the peptide header. */
#include "peptide.h"

/* peptide_add_residue1(): add a residue to a peptide structure based
 * on the single-letter residue code.
 *
 * arguments:
 *  @P: pointer to the peptide structure to modify.
 *  @res: residue single-letter code to add.
 *
 * returns:
 *  integer indicating whether (1) or not (0) the operation succeeded.
 */
int peptide_add_residue1 (peptide_t *P, char res) {
  /* declare required variables:
   *  @ires: residue index 
   */
  unsigned int ires;

  /* lookup the residue index from the code. */
  ires = resid_lookup1(res);

  /* check that the residue is valid. */
  if (!ires)
    throw("invalid residue code '%c'", res);

  /* increment the array length. */
  P->n_res++;

  /* reallocate the array. */
  P->res = (unsigned int*) realloc(P->res, P->n_res * sizeof(unsigned int));
  if (!P->res)
    throw("unable to reallocate sequence array");

  /* store the new residue index. */
  P->res[P->n_res - 1] = ires;

  /* return success. */
  return 1;
}

/* peptide_add_residue3(): add a residue to a peptide structure based
 * on the three-letter residue code.
 *
 * arguments:
 *  @P: pointer to the peptide structure to modify.
 *  @res: residue three-letter code to add.
 *
 * returns:
 *  integer indicating whether (1) or not (0) the operation succeeded.
 */
int peptide_add_residue3 (peptide_t *P, char *res) {
  /* declare required variables:
   *  @ires: residue index 
   *  @i: sequence index.
   */
  unsigned int i, ires;

  /* lookup the residue index from the code. */
  ires = resid_lookup3(res);

  /* check that the residue is valid. */
  if (!ires)
    throw("invalid residue code '%s'", res);

  /* increment the array length. */
  i = P->n_res;
  P->n_res++;

  /* reallocate the array. */
  P->res = (unsigned int*) realloc(P->res, P->n_res * sizeof(unsigned int));
  if (!P->res)
    throw("unable to reallocate sequence array");

  /* store the new residue index. */
  P->res[i] = ires;

  /* return success. */
  return 1;
}

/* peptide_add_sidechain(): request that a given sidechain in a peptide
 * be made atomistically explicit during modeling.
 *
 * arguments:
 *  @P: pointer to the peptide structure to modify.
 *  @res: index of the residue to make explicit.
 *
 * returns:
 *  integer indicating whether (1) or not (0) the function succeeded.
 */
int peptide_add_sidechain (peptide_t *P, unsigned int res) {
  /* declare required variables:
   *  @i: new residue insertion array index.
   *  @j: residue array shift loop counter.
   */
  unsigned int i;
  int j;

  /* check that the residue index is in bounds. */
  if (res >= P->n_res)
    throw("residue index %u out of bounds [1,%u]", res + 1, P->n_res);

  /* return if the peptide already contains the sidechain. */
  if (peptide_has_sidechain(P, res))
    return 1;

  /* increase the sidechain array count. */
  P->n_sc++;

  /* reallocate the sidechain array. */
  P->sc = (unsigned int*) realloc(P->sc, P->n_sc * sizeof(unsigned int));

  /* check if reallocation failed. */
  if (!P->sc)
    throw("unable to reallocate sidechain array");

  /* loop until the insertion point is located. */
  for (i = 0; i < P->n_sc - 1; i++) {
    /* return if the residue has already been made explicit. */
    if (P->sc[i] == res)
      return 1;

    /* break if the insertion point has been identified. */
    if (P->sc[i] > res)
      break;
  }

  /* shift the larger residue indices down in the array to make room. */
  for (j = P->n_sc - 1; j > (int) i; j--)
    P->sc[j] = P->sc[j - 1];

  /* store the new residue index. */
  P->sc[i] = res;

  /* return success. */
  return 1;
}

/* peptide_has_sidechain(): determine whether a peptide structure has
 * a given sidechain explicitly defined.
 *
 * arguments:
 *  @P: pointer to the peptide structure to query.
 *  @res: index of the query residue.
 *
 * returns:
 *  integer indicating whether (1) or not (0) the specified residue
 *  has explicit sidechain atoms.
 */
int peptide_has_sidechain (peptide_t *P, unsigned int res) {
  /* declare required variables:
   *  @imin, @imax: lower and upper binary search bounds.
   *  @i: midpoint sidechain array index of the search.
   */
  int imin, i, imax;

  /* handle the trivial cases. */
  if (!P || P->n_sc == 0)
    return 0;

  /* initialize the binary search bounds. */
  imin = i = 0;
  imax = P->n_sc - 1;

  /* loop until the bounds have converged. */
  while (imin <= imax) {
    /* locate the midpoint of the search bounds. */
    i = (imin + imax) / 2;

    /* return true on equality, or refine the search bounds. */
    if (P->sc[i] == res)
      return 1;
    else if (P->sc[i] < res)
      imin = i + 1;
    else
      imax = i - 1;
  }

  /* residue index not found. return false. */
  return 0;
}

/* peptide_get_reschar(): get the residue name character of a given indexed
 * residue within a peptide structure.
 *
 * arguments:
 *  @P: pointer to the peptide structure to access.
 *  @res: residue index to retrieve name information for.
 *
 * returns:
 *  constant residue name character.
 */
char peptide_get_reschar (peptide_t *P, unsigned int res) {
  /* always return the single-letter residue code. */
  return (P && res < P->n_res ? resid_get_code1(P->res[res]) : 'X');
}

/* peptide_get_resname(): get the residue name string of a given indexed
 * residue within a peptide structure.
 *
 * arguments:
 *  @P: pointer to the peptide structure to access.
 *  @res: residue index to retrieve name information for.
 *
 * returns:
 *  constant residue name string, or NULL on failure.
 */
const char *peptide_get_resname (peptide_t *P, unsigned int res) {
  /* always return the three-letter residue code. */
  return (P && res < P->n_res ? resid_get_code3(P->res[res]) : NULL);
}

/* peptide_get_restype(): get the residue type string of a given indexed
 * residue within a peptide structure.
 *
 * arguments:
 *  @P: pointer to the peptide structure to access.
 *  @res: residue index to retrieve type information for.
 *
 * returns:
 *  constant residue type string, or NULL on failure.
 */
const char *peptide_get_restype (peptide_t *P, unsigned int res) {
  /* check that all arguments are valid. */
  if (!P || res >= P->n_res)
    return NULL;

  /* determine whether the residue has an explicit sidechain. */
  if (peptide_has_sidechain(P, res)) {
    /* yes, return the three-letter residue code. */
    return resid_get_code3(P->res[res]);
  }
  else {
    /* no, return a stock residue type. */
    if (res == 0)
      return "BB1";
    else if (res == 1)
      return "BB2";
    else if (res == P->n_res - 1)
      return "BBN";
    else
      return "BBI";
  }
}

