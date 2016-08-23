
/* ensure once-only inclusion. */
#ifndef __IBPNG_REORDER_H__
#define __IBPNG_REORDER_H__

/* include the traceback and string function headers. */
#include "trace.h"
#include "str.h"

/* reorder_atom_t: structure for holding a single atom within a repetition
 * ordering sequence.
 */
typedef struct {
  /* @name: name string of the atom.
   * @off: residue index offset of the atom.
   * @opt: whether the atom is optional or not.
   */
  char *name;
  int off, opt;
}
reorder_atom_t;

/* reorder_t: structure for holding repetition ordering information.
 */
typedef struct _reorder_t reorder_t;
struct _reorder_t {
  /* @prev, @next: pointers to the previous and next entries in the
   * reorder [linked list] data structure.
   */
  reorder_t *prev, *next;

  /* @name: reorder name string.
   */
  char *name;

  /* @atoms: array of atoms in the current reorder entry.
   * @n_atoms: number of atoms in the reorder entry.
   */
  reorder_atom_t *atoms;
  unsigned int n_atoms;
};

/* function declarations: */

reorder_t *reorder_new (void);

reorder_t *reorder_new_from_file (const char *fname);

void reorder_free (reorder_t *ord);

int reorder_add_residue (reorder_t *ord, const char *name);

int reorder_add_atom (reorder_t *ord, const char *name,
                      int offset, int optional);

reorder_t *reorder_get_residue (reorder_t *ord, const char *name);

#endif  /* !__IBPNG_REORDER_H__ */

