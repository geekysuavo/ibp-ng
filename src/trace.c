
/* include the traceback header. */
#include "trace.h"

/* traceback_t: global data structure for holding stack trace information.
 *
 * as errors are emitted from functions using the throw() macro function,
 * they will be stored in the error stack until the error propagates back
 * to a die() macro function, which will print the stack trace and end
 * program execution.
 */
typedef struct {
  /* @line: the line number where the error was emitted.
   * @file: the source file that emitted the error.
   * @num: the emitted error code.
   */
  unsigned int line;
  char *file;
  int num;

  /* @msg: the custom error message string, or NULL if none. */
  char *msg;
}
traceback_t;

/* tb: a variable known only to the traceback_* functions herein, which
 * holds the current stack trace of the currently executing application.
 */
traceback_t *tb;
int n_tb = -1;

/* tb_verb: a local variable used by traceback_* functions to determine
 * whether or not to output verbose information messages and warnings.
 */
int tb_verb = 0;

/* traceback_init(): initialize the traceback array for use.
 *
 * this function is called automatically by traceback_throw(), which is
 * itself called by raise(), throw() and die().
 */
void traceback_init (void) {
  /* check if the traceback array has been initialized. */
  if (n_tb >= 0)
    return;

  /* initialize the traceback array and set its size to zero. */
  tb = NULL;
  n_tb = 0;
}

/* verbosity_set(): set the verbosity level of the application.
 *
 * arguments:
 *  @level: verbosity level to set.
 */
void verbosity_set (unsigned int level) {
  /* set the verbosity level. */
  tb_verb = (int) level;

  /* print some info. */
  warn("verbose warnings enabled");
  info("verbose messages enabled");
}

/* traceback_length(): return the number of entries in the traceback array.
 *
 * returns:
 *  number of entries (if any) in the traceback array.
 */
unsigned int traceback_length (void) {
  /* ensure the traceback array is initialized. */
  traceback_init();

  /* return the traceback length. */
  return (unsigned int) n_tb;
}

/* traceback_print(): print the contents of the traceback array to standard
 * error.
 */
void traceback_print (void) {
  /* declare required variables:
   *  @i: traceback array index.
   */
  int i;

  /* ensure the traceback array is initialized. */
  traceback_init();

  /* loop over the traceback array entries. */
  for (i = 0; i < n_tb; i++) {
    /* print the traceback line header. */
    fprintf(stderr, "[%d] %s:%u:\n", i, tb[i].file, tb[i].line);

    /* check if a custom message was set at the time of the error. */
    if (tb[i].msg)
      fprintf(stderr, "  %s\n", tb[i].msg);

    /* check if errno was set at the time of the error. */
    if (tb[i].num)
      fprintf(stderr, "  (%s)\n", strerror(tb[i].num));
  }
}

/* traceback_clear(): clear the contents of the traceback array.
 *
 * this function is useful for catching errors that are expected and
 * non-fatal to application execution.
 */
void traceback_clear (void) {
  /* declare required variables:
   *  @i: traceback array index.
   */
  int i;

  /* ensure the traceback array is initialized. */
  traceback_init();

  /* loop over the traceback array entries. */
  for (i = 0; i < n_tb; i++) {
    /* free the allocated strings. */
    if (tb[i].file) free(tb[i].file);
    if (tb[i].msg) free(tb[i].msg);
  }

  /* free the traceback array. */
  free(tb);

  /* prepare for the next error. */
  errno = 0;
  tb = NULL;
  n_tb = 0;
}

/* traceback_throw(): push a new frame of error information onto the
 * traceback error stack.
 *
 * arguments:
 *  @f: filename from where the error was emitted.
 *  @l: line number in the file where the error was emitted.
 *  @format: printf-style format string for custom messages.
 *  @...: arguments that correspond to the @format string.
 *
 * returns:
 *  zero. this is used by the throw() macro function to indicate failure
 *  from any function returning an integer success code.
 */
int traceback_throw (const char *f,
                     const unsigned int l,
                     const char *format, ...) {
  /* declare required variables:
   *  @n_msg: length of the filename and message strings.
   *  @i: traceback array index of the appended error.
   *  @n_print: amount of characters written by vsnprintf.
   *  @vl: variable arguments list for custom error messages.
   */
  int i, n_msg, n_print;
  va_list vl;

  /* ensure the traceback array is initialized. */
  traceback_init();

  /* increment the traceback array size. */
  i = n_tb;
  n_tb++;

  /* reallocate the traceback array. */
  tb = (traceback_t*) realloc(tb, n_tb * sizeof(traceback_t));

  /* check if the traceback array was successfully reallocated. */
  if (!tb) {
    /* it was not. output an extra error message... */
    fprintf(stderr, "[x] %s:%u:\n", __FILE__, __LINE__);
    fprintf(stderr, "  unable to reallocate traceback array\n");

    /* ...and return failure. */
    return 0;
  }

  /* add the filename string to the traceback entry. */
  n_msg = 2 * strlen(f);
  tb[i].file = (char*) malloc(n_msg * sizeof(char));
  if (tb[i].file)
    strcpy(tb[i].file, f);

  /* add the numeric information to the traceback entry. */
  tb[i].num = errno;
  tb[i].line = l;

  /* build the custom error message string. */
  if (format) {
    /* begin with a guess of the output string length. */
    n_msg = strlen(format);
    tb[i].msg = NULL;

    /* loop until the whole string was printed. */
    do {
      /* allocate memory for the message. */
      n_msg *= 2;
      tb[i].msg = (char*) realloc(tb[i].msg, n_msg * sizeof(char));

      /* check if the message string was successfully reallocated. */
      if (!tb[i].msg) {
        /* it was not. output an extra error message. */
        fprintf(stderr, "[x] %s:%u:\n", __FILE__, __LINE__);
        fprintf(stderr, "  unable to reallocate error message string\n");

        /* ...and return failure. */
        return 0;
      }

      /* write the formatted message string. */
      va_start(vl, format);
      n_print = vsnprintf(tb[i].msg, n_msg, format, vl);
      va_end(vl);
    }
    while (n_print >= n_msg);
  }

  /* clear the numeric error code, as it's now safely stored. */
  errno = 0;

  /* return failure. */
  return 0;
}

/* verbose_print(): print a message string to standard output without
 * changing the state of the traceback error stack. used by the info()
 * and warn() macro functions.
 *
 * arguments:
 *  @lev: verbosity level of the message.
 *  @f: filename from where the error was emitted.
 *  @l: line number in the file where the error was emitted.
 *  @format: printf-style format string for custom messages.
 *  @...: arguments that correspond to the @format string.
 */
void verbose_print (const int lev,
                    const char *f,
                    const unsigned int l,
                    const char *format, ...) {
  /* declare required variables:
   *  @vl: variable arguments list for custom error messages.
   */
  va_list vl;

  /* return if the verbosity level is too low. */
  if (tb_verb < lev)
    return;

  /* print the first part of the message. */
  fprintf(stderr, "[%s] %s:%u:\n  ",
          lev == IBPNG_VERBOSE_WARN  ? "WARN" :
          lev == IBPNG_VERBOSE_INFO  ? "INFO" :
          lev == IBPNG_VERBOSE_DEBUG ? "DBUG" :
          "????", f, l);

  /* print the second part of the message. */
  va_start(vl, format);
  vfprintf(stderr, format, vl);
  va_end(vl);

  /* flush the output buffer. */
  fprintf(stderr, "\n");
  fflush(stderr);
}

