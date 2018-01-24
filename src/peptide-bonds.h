
/* ensure once-only inclusion. */
#pragma once

/* function declarations (peptide-bonds.c): */

int peptide_bond_add (peptide_t *P,
                      unsigned int resid1, const char *name1,
                      unsigned int resid2, const char *name2,
                      unsigned int is_virtual);

int peptide_bond_delete (peptide_t *P,
                         unsigned int resid1, const char *name1,
                         unsigned int resid2, const char *name2);

int peptide_bond_delete_any (peptide_t *P,
                             unsigned int resid,
                             const char *name);

int peptide_field_bonds (peptide_t *P, double tol);

int peptide_graph_bonds (peptide_t *P, graph_t *G);

