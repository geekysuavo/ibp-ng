
/* include the string functions header. */
#include "str.h"

/* strdup(): duplicate a string.
 *
 * arguments:
 *  @s: string to duplicate.
 *
 * returns:
 *  pointer to a newly duplicated string, or NULL on failure.
 */
char *strdup (const char *s) {
  /* declare required variables:
   *  @snew: newly allocated string.
   */
  char *snew;

  /* return null if the input string is null. */
  if (!s)
    return NULL;

  /* allocate memory for the duplicate string. */
  snew = (char*) malloc((strlen(s) + 1) * sizeof(char));

  /* return null if allocation failed. */
  if (!snew)
    return NULL;

  /* copy the string contents. */
  strcpy(snew, s);

  /* return the duplicate string. */
  return snew;
}

/* strndup(): duplicate a string of a specified length. this function
 * is dangerous. use it wisely.
 *
 * arguments:
 *  @s: string to duplicate.
 *  @n: string length.
 *
 * returns:
 *  pointer to a newly duplicated string, or NULL on failure.
 */
char *strndup (const char *s, unsigned int n) {
  /* declare required variabies:
   *  @snew: newly allocated string.
   */
  char *snew;

  /* return null if the input string is null. */
  if (!s)
    return NULL;

  /* allocate memory for the duplicate string. */
  snew = (char*) malloc((n + 1) * sizeof(char));

  /* return null if allocation failed. */
  if (!snew)
    return NULL;

  /* copy the string contents. */
  strncpy(snew, s, n);
  snew[n] = '\0';

  /* return the duplicate string. */
  return snew;
}

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

