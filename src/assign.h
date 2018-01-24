
/* ensure once-only inclusion. */
#pragma once

/* include the peptide and string headers. */
#include "peptide.h"
#include "str.h"

/* assign_set_t: structure for holding a peptide assignment atom set.
 */
typedef struct {
  /* @v: array of atom indices contained within the set.
   * @n: number of atom indices contained by the set.
   * @n_max: maximum size of the atom index set.
   */
  unsigned int *v, n, n_max;

  /* @P: pointer to a peptide structure to access for all operations.
   */
  peptide_t *P;
}
assign_set_t;

/* function declarations: */

void assign_free (assign_set_t *set);

assign_set_t *assign_none (peptide_t *P);

assign_set_t *assign_all (peptide_t *P);

assign_set_t *assign_atomid (peptide_t *P, unsigned int id);

assign_set_t *assign_resid (peptide_t *P, unsigned int id);

assign_set_t *assign_resname (peptide_t *P, const char *name);

assign_set_t *assign_name (peptide_t *P, const char *name);

assign_set_t *assign_type (peptide_t *P, char type);

assign_set_t *assign_not (assign_set_t *set1);

assign_set_t *assign_and (assign_set_t *set1, assign_set_t *set2);

assign_set_t *assign_or (assign_set_t *set1, assign_set_t *set2);

int assign_set_distance (peptide_t *P,
                         assign_set_t *set1,
                         assign_set_t *set2,
                         double d, double dmin, double dplus);

int assign_set_dihedral (peptide_t *P,
                         assign_set_t *set1,
                         assign_set_t *set2,
                         assign_set_t *set3,
                         assign_set_t *set4,
                         double phi, double dphi);

int assign_set_from_file (peptide_t *P, const char *fname);

