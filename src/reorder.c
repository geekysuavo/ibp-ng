
/* include the repetition ordering header. */
#include "reorder.h"

/* function declarations from parsing code: */
int reorder_parse (FILE *fh, reorder_t *ord);

/* reorder_new(): allocate a new repetition ordering data structure.
 *
 * returns:
 *  pointer to a newly allocated empty data structure.
 */
reorder_t *reorder_new (void) {
  /* declare required variables:
   *  @ord: output structure pointer.
   */
  reorder_t *ord;

  /* allocate a new structure pointer. */
  ord = (reorder_t*) malloc(sizeof(reorder_t));

  /* check if allocation failed. */
  if (!ord) {
    /* raise an exception and return null. */
    raise("unable to allocate reorder structure pointer");
    return NULL;
  }

  /* initialize the linkages. */
  ord->prev = ord->next = NULL;

  /* initialize the link name. */
  ord->name = NULL;

  /* initialize the atoms. */
  ord->atoms = NULL;
  ord->n_atoms = 0;

  /* return the structure pointer. */
  return ord;
}

/* reorder_free(): free all allocated memory associated with a repetition
 * ordering data structure.
 *
 * arguments:
 *  @ord: pointer to the data structure to free.
 */
void reorder_free (reorder_t *ord) {
  /* declare required variables:
   *  @lnk, @nxt: reorder pointers for traversal.
   *  @i: general purpose loop counter.
   */
  reorder_t *lnk, *nxt;
  unsigned int i;

  /* return if the structure pointer is null. */
  if (!ord) return;

  /* move to the first link in the chain. */
  lnk = ord;
  while (lnk->prev)
    lnk = lnk->prev;

  /* traverse the chain to free all links. */
  while (lnk) {
    /* free the link name. */
    if (lnk->name)
      free(lnk->name);

    /* check if the link contains atoms. */
    if (lnk->n_atoms) {
      /* free the atom name strings. */
      for (i = 0; i < lnk->n_atoms; i++) {
        if (lnk->atoms[i].name)
          free(lnk->atoms[i].name);
      }

      /* free the atom array. */
      free(lnk->atoms);
      lnk->n_atoms = 0;
    }

    /* save the next link. */
    nxt = lnk->next;

    /* free the structure pointer and move to the next link. */
    free(lnk);
    lnk = nxt;
  }
}

/* reorder_new_from_file(): allocate a new repetition ordering data structure
 * by reading information from a file.
 *
 * arguments:
 *  @fname: input ordering filename string.
 *
 * returns:
 *  pointer to a newly allocated and filled data structure.
 */
reorder_t *reorder_new_from_file (const char *fname) {
  /* declare required variables:
   *  @ord: output structure pointer.
   *  @fh: input file handle.
   */
  reorder_t *ord;
  FILE *fh;

  /* allocate the structure pointer. */
  ord = reorder_new();
  if (!ord)
    return NULL;

  /* open the input file. */
  fh = fopen(fname, "r");

  /* check that the file was opened successfully. */
  if (!fh) {
    /* no. raise an exception and return null. */
    raise("unable to open '%s' for reading", fname);
    reorder_free(ord);
    return NULL;
  }

  /* attempt to parse the input file. */
  if (!reorder_parse(fh, ord)) {
    /* failure. raise an exception. */
    raise("unable to parse reorder information");

    /* clean up. */
    reorder_free(ord);
    fclose(fh);

    /* return null. */
    return NULL;
  }

  /* close the input file. */
  fclose(fh);

  /* return the structure pointer. */
  return ord;
}

/* reorder_add_residue(): add an empty residue into a repetition order.
 *
 * arguments:
 *  @ord: pointer to the reorder structure to modify.
 *  @name: name of the new residue.
 *
 * returns:
 *  integer indicating whether (1) or not (0) the operation succeeded.
 */
int reorder_add_residue (reorder_t *ord, const char *name) {
  /* declare required variables:
   *  @lnk: reorder pointer for traversal.
   */
  reorder_t *lnk;

  /* fail if the structure pointer is null. */
  if (!ord)
    throw("reorder structure pointer is null");

  /* locate the last element of the structure. */
  lnk = ord;
  while (lnk->next)
    lnk = lnk->next;

  /* check if the element has a name. */
  if (lnk->name) {
    /* yes. allocate a new link. */
    lnk->next = reorder_new();
    if (!lnk->next)
      return 0;

    /* move to the newly allocated link. */
    lnk->next->prev = lnk;
    lnk = lnk->next;
  }

  /* store the link name. */
  lnk->name = (char*) name;

  /* return success. */
  return 1;
}

/* reorder_add_atom(): add an atom into the last residue of a repetition
 * order.
 *
 * arguments:
 *  @ord: pointer to the reorder structure to modify.
 *  @name: name of the atom to add into the order.
 *  @offset: residue offset of the atom to add.
 *  @optional: whether the atom is optional.
 *
 * returns:
 *  integer indicating whether (1) or not (0) the operation succeeded.
 */
int reorder_add_atom (reorder_t *ord, const char *name,
                      int offset, int optional) {
  /* declare required variables:
   *  @lnk: reorder pointer for traversal.
   *  @i: atom array index.
   */
  reorder_t *lnk;
  unsigned int i;

  /* fail if the structure pointer is null. */
  if (!ord)
    throw("reorder structure pointer is null");

  /* locate the last element of the structure. */
  lnk = ord;
  while (lnk->next)
    lnk = lnk->next;

  /* increase the size of the atom array. */
  i = lnk->n_atoms;
  lnk->n_atoms++;

  /* reallocate the atom array. */
  lnk->atoms = (reorder_atom_t*)
    realloc(lnk->atoms, lnk->n_atoms * sizeof(reorder_atom_t));

  /* check that reallocation succeeded. */
  if (!lnk->atoms)
    throw("unable to reallocate atom array");

  /* store the new atom in the array. */
  lnk->atoms[i].name = (char*) name;
  lnk->atoms[i].off = offset;
  lnk->atoms[i].opt = optional;

  /* return success. */
  return 1;
}

/* reorder_get_residue(): locate a particular named repetition order group
 * from a reorder structure.
 *
 * arguments:
 *  @ord: pointer to the reorder structure to access.
 *  @name: residue/group name string.
 *
 * returns:
 *  pointer to the identified reorder structure, or NULL if no match.
 */
reorder_t *reorder_get_residue (reorder_t *ord, const char *name) {
  /* declare required variables:
   *  @lnk: reorder pointer for traversal.
   */
  reorder_t *lnk;

  /* fail if the structure pointer is null. */
  if (!ord)
    return NULL;

  /* locate the first element of the structure. */
  lnk = ord;
  while (lnk->prev)
    lnk = lnk->prev;

  /* loop over the elements of the structure. */
  while (lnk) {
    /* check if the current link matches our name string. */
    if (strcmp(lnk->name, name) == 0)
      return lnk;

    /* move to the next link. */
    lnk = lnk->next;
  }

  /* no match, return null. */
  return NULL;
}

