
/* ensure once-only inclusion. */
#pragma once

/* include the peptide and hash headers. */
#include "peptide.h"
#include "dmdgp-hash.h"

/* function declarations: */

int dmdgp_write (const char *fname, peptide_t *P, graph_t *G);

