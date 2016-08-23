
/* ensure once-only inclusion. */
#ifndef __IBPNG_STR_H__
#define __IBPNG_STR_H__

/* include required standard c library headers. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

/* function declarations: */

char *strdup (const char *s);

char *strndup (const char *s, unsigned int n);

char *strtolower (const char *s);

char *strtoupper (const char *s);

#endif  /* !__IBPNG_STR_H__ */

