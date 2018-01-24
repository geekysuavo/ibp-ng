
/* ensure once-only inclusion. */
#pragma once

/* function declarations (peptide-impropers.c): */

int peptide_improper_add (peptide_t *P,
                          unsigned int resid1, const char *name1,
                          unsigned int resid2, const char *name2,
                          unsigned int resid3, const char *name3,
                          unsigned int resid4, const char *name4);

int peptide_improper_delete (peptide_t *P,
                             unsigned int resid1, const char *name1,
                             unsigned int resid2, const char *name2,
                             unsigned int resid3, const char *name3,
                             unsigned int resid4, const char *name4);

int peptide_improper_delete_any (peptide_t *P,
                                 unsigned int resid,
                                 const char *name);

int peptide_field_impropers (peptide_t *P, double tol);

int peptide_graph_impropers (peptide_t *P, graph_t *G);

