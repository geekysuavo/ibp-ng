
/* ensure once-only inclusion. */
#ifndef __IBPNG_DMDGP_HASH_H__
#define __IBPNG_DMDGP_HASH_H__

/* include the traceback and string headers. */
#include "trace.h"
#include "str.h"

/* dmdgp_hash_t: data structure for storing a list of integers that are
 * each associated with a much smaller (i.e. highly overlapping) set of
 * string values.
 */
typedef struct {
  /* @n: number of key-value pairs in the hash.
   */
  unsigned int n;

  /* @keys: array of key strings.
   * @nums: number of values for each key.
   * @vals: array of value integer arrays.
   */
  char **keys;
  unsigned int *nums;
  unsigned int **vals;
}
dmdgp_hash_t;

/* function declarations: */

dmdgp_hash_t *dmdgp_hash_new (void);

void dmdgp_hash_free (dmdgp_hash_t *hash);

int dmdgp_hash_add (dmdgp_hash_t *hash, const char *key, unsigned int val);

int dmdgp_hash_write (dmdgp_hash_t *hash, const char *fmt, FILE *fh);

#endif  /* !__IBPNG_DMDGP_HASH_H__ */

