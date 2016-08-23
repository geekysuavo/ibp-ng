
/* ensure once-only inclusion. */
#ifndef __IBPNG_TRACE_H__
#define __IBPNG_TRACE_H__

/* include required standard c library headers. */
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <errno.h>

/* verbosity level definitions. */
#define IBPNG_SILENT        0
#define IBPNG_VERBOSE_WARN  1
#define IBPNG_VERBOSE_INFO  2

/* raise(): macro function to add traceback information onto the
 * application's error stack.
 */
#define raise(...) \
  traceback_throw(__FILE__, __LINE__, __VA_ARGS__)

/* throw(): macro function to add traceback information onto the
 * application's error stack and return zero.
 */
#define throw(...) \
  return traceback_throw(__FILE__, __LINE__, __VA_ARGS__)

/* die(): macro function to add traceback information onto the
 * application's error stack, print the formatted error stack to
 * standard error, and return an error for main().
 */
#define die(...) \
  { traceback_throw(__FILE__, __LINE__, __VA_ARGS__); \
    traceback_print(); \
    goto death; }

/* info(): macro function to print a verbose informational message to
 * standard error.
 */
#define info(...) \
  verbose_print(IBPNG_VERBOSE_INFO, __FILE__, __LINE__, __VA_ARGS__)

/* warn(): macro function to print a verbose warning message to
 * standard error.
 */
#define warn(...) \
  verbose_print(IBPNG_VERBOSE_WARN, __FILE__, __LINE__, __VA_ARGS__)

/* function declarations: */

void verbosity_set (unsigned int level);

unsigned int traceback_length (void);

void traceback_print (void);

void traceback_clear (void);

int traceback_throw (const char *f,
                     const unsigned int l,
                     const char *format, ...);

void verbose_print (const int lev,
                    const char *f,
                    const unsigned int l,
                    const char *format, ...);

#endif  /* !__IBPNG_TRACE_H__ */

