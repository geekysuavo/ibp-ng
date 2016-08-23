
/* ensure once-only inclusion. */
#ifndef __IBPNG_ENUM_WRITE_H__
#define __IBPNG_ENUM_WRITE_H__

/* function declarations (DCD): */

int enum_write_dcd_open (enum_t *E);

void enum_write_dcd_close (enum_t *E);

int enum_write_dcd (enum_t *E, enum_thread_t *th);

/* function declarations (PDB): */

int enum_write_pdb_open (enum_t *E);

int enum_write_pdb (enum_t *E, enum_thread_t *th);

#endif /* !__IBPNG_ENUM_WRITE_H__ */

