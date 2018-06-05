
/* ensure once-only inclusion. */
#pragma once

/* include required standard c library headers. */
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <errno.h>

/* verbosity level definitions. */
#define IBPNG_SILENT         0
#define IBPNG_VERBOSE_WARN   1
#define IBPNG_VERBOSE_INFO   2
#define IBPNG_VERBOSE_DEBUG  3

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

/* warn(): macro function to print a verbose warning message to
 * standard error.
 */
#define warn(...) \
  verbose_print(IBPNG_VERBOSE_WARN, __FILE__, __LINE__, __VA_ARGS__)

/* info(): macro function to print a verbose informational message to
 * standard error.
 */
#define info(...) \
  verbose_print(IBPNG_VERBOSE_INFO, __FILE__, __LINE__, __VA_ARGS__)

/* debug(): macro function to print a verbose debugging message to
 * standard error.
 */
#define debug(...) \
  verbose_print(IBPNG_VERBOSE_DEBUG, __FILE__, __LINE__, __VA_ARGS__)

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

