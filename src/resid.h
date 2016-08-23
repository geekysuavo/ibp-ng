
/* ensure once-only inclusion. */
#ifndef __IBPNG_RESID_H__
#define __IBPNG_RESID_H__

/* include the required standard c library headers. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* include the string header. */
#include "str.h"

/* function declarations: */

unsigned int resid_lookup1 (char code);

unsigned int resid_lookup3 (const char *code);

char resid_get_code1 (unsigned int i);

const char *resid_get_code3 (unsigned int i);

unsigned int resid_get_nr (unsigned int i);

#endif /* !__IBPNG_RESID_H__ */

