
/* ensure once-only inclusion. */
#ifndef __IBPNG_DMDGP_H__
#define __IBPNG_DMDGP_H__

/* include the peptide and hash headers. */
#include "peptide.h"
#include "dmdgp-hash.h"

/* function declarations: */

int dmdgp_write (const char *fname, peptide_t *P, graph_t *G);

#endif  /* !__IBPNG_DMDGP_H__ */

