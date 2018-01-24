
/* ensure once-only inclusion. */
#pragma once

/* function declarations (DCD): */

int enum_write_dcd_open (enum_t *E);

void enum_write_dcd_close (enum_t *E);

int enum_write_dcd (enum_t *E, enum_thread_t *th);

/* function declarations (PDB): */

int enum_write_pdb_open (enum_t *E);

int enum_write_pdb (enum_t *E, enum_thread_t *th);

