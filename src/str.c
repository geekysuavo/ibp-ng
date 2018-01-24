
/* include the string functions header. */
#include "str.h"

/* strtolower(): duplicate a string while converting to lowercase.
 *
 * arguments:
 *  @s: string to duplicate.
 *
 * returns:
 *  pointer to a newly duplicated string, or NULL on failure.
 */
char *strtolower (const char *s) {
  /* declare required variables:
   *  @snew: newly allocated string.
   *  @i: character index.
   */
  char *snew;
  unsigned int i;

  /* duplicate the string. */
  snew = strdup(s);
  if (!snew)
    return NULL;

  /* change the string to lowercase. */
  for (i = 0; i < strlen(snew); i++)
    snew[i] = (char) tolower ((int) snew[i]);

  /* return the duplicate string. */
  return snew;
}

/* strtoupper(): duplicate a string while converting to uppercase.
 *
 * arguments:
 *  @s: string to duplicate.
 *
 * returns:
 *  pointer to a newly duplicated string, or NULL on failure.
 */
char *strtoupper (const char *s) {
  /* declare required variables:
   *  @snew: newly allocated string.
   *  @i: character index.
   */
  char *snew;
  unsigned int i;

  /* duplicate the string. */
  snew = strdup(s);
  if (!snew)
    return NULL;

  /* change the string to lowercase. */
  for (i = 0; i < strlen(snew); i++)
    snew[i] = (char) toupper ((int) snew[i]);

  /* return the duplicate string. */
  return snew;
}

