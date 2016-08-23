
/* include the dmdgp hash header. */
#include "dmdgp-hash.h"

/* dmdgp_hash_new(): allocate a new, empty hash structure.
 *
 * returns:
 *  pointer to a newly allocated and initialized hash data structure, or
 *  NULL if allocation failed.
 */
dmdgp_hash_t *dmdgp_hash_new (void) {
  /* declare required variables:
   *  @hash: output structure pointer.
   */
  dmdgp_hash_t *hash;

  /* allocate a new structure pointer. */
  hash = (dmdgp_hash_t*) malloc(sizeof(dmdgp_hash_t));

  /* check if allocation failed. */
  if (!hash) {
    /* raise an exception and return NULL. */
    raise("unable to allocate hash structure pointer");
    return NULL;
  }

  /* initialize the structure contents. */
  hash->n = 0;
  hash->keys = NULL;
  hash->nums = NULL;
  hash->vals = NULL;

  /* return the structure pointer. */
  return hash;
}

/* dmdgp_hash_free(): free all memory associated with a hash structure.
 *
 * arguments:
 *  @hash: pointer to the hash structure to free.
 */
void dmdgp_hash_free (dmdgp_hash_t *hash) {
  /* declare required variables:
   *  @i: array loop counter.
   */
  unsigned int i;

  /* return if the structure pointer is null. */
  if (!hash) return;

  /* check if the structure contains elements. */
  if (hash->n) {
    /* loop over the element arrays. */
    for (i = 0; i < hash->n; i++) {
      /* free the key string and value array. */
      free(hash->keys[i]);
      free(hash->vals[i]);
    }

    /* free the arrays. */
    free(hash->keys);
    free(hash->nums);
    free(hash->vals);
  }

  /* finally, free the structure pointer. */
  free(hash);
}

/* dmdgp_hash_add(): add a key-value pair to a dmdgp hash structure.
 *
 * arguments:
 *  @hash: pointer to the hash structure to modify.
 *  @key: key string of the new pair.
 *  @val: value integer of the new pair.
 *
 * returns:
 *  integer indicating whether (1) or not (0) the operation succeeded.
 */
int dmdgp_hash_add (dmdgp_hash_t *hash, const char *key, unsigned int val) {
  /* declare required variables:
   *  @i: key index.
   *  @j: value index.
   */
  unsigned int i, j;

  /* loop over the keys to determine whether the new key exists. */
  for (i = 0; i < hash->n; i++) {
    /* break if the key is found. */
    if (strcmp(hash->keys[i], key) == 0)
      break;
  }

  /* check whether the key was found. */
  if (i >= hash->n) {
    /* no. resize the hash to include the new key. */
    i = hash->n;
    hash->n++;

    /* reallocate the arrays. */
    hash->keys = (char**) realloc(hash->keys, hash->n * sizeof(char*));
    hash->nums = (unsigned int*)
      realloc(hash->nums, hash->n * sizeof(unsigned int));
    hash->vals = (unsigned int**)
      realloc(hash->vals, hash->n * sizeof(unsigned int*));

    /* check if reallocation failed. */
    if (!hash->keys || !hash->nums || !hash->vals)
      throw("unable to reallocate hash arrays");

    /* store the new key. */
    hash->keys[i] = strdup(key);
    if (!hash->keys[i])
      throw("unable to store new hash key '%s'", key);

    /* initialize the values portion. */
    hash->nums[i] = 0;
    hash->vals[i] = NULL;
  }

  /* resize the values array of the current key. */
  j = hash->nums[i];
  hash->nums[i]++;

  /* reallocate the values array. */
  hash->vals[i] = (unsigned int*)
    realloc(hash->vals[i], hash->nums[i] * sizeof(unsigned int));

  /* check if reallocation failed. */
  if (!hash->vals[i])
    throw("unable to reallocate hash value array");

  /* store the new value. */
  hash->vals[i][j] = val;

  /* return success. */
  return 1;
}

/* dmdgp_hash_write(): write the contents of a hash structure to a file.
 *
 * arguments:
 *  @hash: pointer to the hash structure to access.
 *  @fmt: output format string for each value.
 *  @fh: output file handle.
 *
 * returns:
 *  integer indicating whether (1) or not (0) the write succeeded.
 */
int dmdgp_hash_write (dmdgp_hash_t *hash, const char *fmt, FILE *fh) {
  /* declare required variables:
   *  @i: key index.
   *  @j: value index.
   */
  unsigned int i, j;

  /* loop over the hash keys. */
  for (i = 0; i < hash->n; i++) {
    /* print the key. */
    fprintf(fh, "%-4s ", hash->keys[i]);

    /* loop over the values of the key. */
    for (j = 0; j < hash->nums[i]; j++) {
      /* print the value. */
      fprintf(fh, fmt, hash->vals[i][j]);
    }

    /* print a newline. */
    fprintf(fh, "\n");
  }

  /* return success. */
  return 1;
}

