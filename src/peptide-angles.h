
/* ensure once-only inclusion. */
#pragma once

/* function declarations (peptide-angles.c): */

int peptide_angle_add (peptide_t *P,
                       unsigned int resid1, const char *name1,
                       unsigned int resid2, const char *name2,
                       unsigned int resid3, const char *name3);

int peptide_angle_delete (peptide_t *P,
                          unsigned int resid1, const char *name1,
                          unsigned int resid2, const char *name2,
                          unsigned int resid3, const char *name3);

int peptide_angle_delete_any (peptide_t *P,
                              unsigned int resid,
                              const char *name);

int peptide_field_angles (peptide_t *P, double tol);

int peptide_graph_angles (peptide_t *P, graph_t *G);

