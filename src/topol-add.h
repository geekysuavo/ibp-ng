
/* ensure once-only inclusion. */
#ifndef __IBPNG_TOPOL_ADD_H__
#define __IBPNG_TOPOL_ADD_H__

/* function declarations (topol-add.c): */

int topol_add_mass (topol_t *top, const char *type, double mass);

int topol_add_residue (topol_t *top, const char *name, unsigned int patch);

int topol_add_atom (topol_t *top,
                    const char *resname,
                    const char *name,
                    const char *type,
                    double charge,
                    topol_mode_t mode, int off);

int topol_add_bond (topol_t *top,
                    const char *resname,
                    const char *a, int aoff,
                    const char *b, int boff,
                    topol_mode_t mode);

int topol_add_angle (topol_t *top,
                     const char *resname,
                     const char *a, int aoff,
                     const char *b, int boff,
                     const char *c, int coff,
                     topol_mode_t mode);

int topol_add_torsion (topol_t *top,
                       const char *resname,
                       const char *a, int aoff,
                       const char *b, int boff,
                       const char *c, int coff,
                       const char *d, int doff,
                       topol_mode_t mode);

int topol_add_improper (topol_t *top,
                        const char *resname,
                        const char *a, int aoff,
                        const char *b, int boff,
                        const char *c, int coff,
                        const char *d, int doff,
                        topol_mode_t mode);

#endif /* !__IBPNG_TOPOL_ADD_H__ */

