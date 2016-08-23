
/* include the molecular parameters header. */
#include "param.h"

/* function declarations from parsing code: */
int param_parse (FILE *fh, param_t *par);

/* param_new(): allocate a new molecular parameters data structure.
 *
 * returns:
 *  pointer to a newly allocated empty data structure.
 */
param_t *param_new (void) {
  /* declare required variables:
   *  @par: output structure pointer.
   */
  param_t *par;

  /* allocate a new structure pointer. */
  par = (param_t*) malloc(sizeof(param_t));

  /* check if allocation failed. */
  if (!par) {
    /* raise an exception and return null. */
    raise("unable to allocate parameter structure pointer");
    return NULL;
  }

  /* initialize the bond array. */
  par->bonds = NULL;
  par->n_bonds = 0;

  /* initialize the angle array. */
  par->angles = NULL;
  par->n_angles = 0;

  /* initialize the torsion array. */
  par->torsions = NULL;
  par->n_torsions = 0;

  /* initialize the improper array. */
  par->impropers = NULL;
  par->n_impropers = 0;

  /* initialize the radius array. */
  par->radii = NULL;
  par->n_radii = 0;
  par->vdw_scale = 1.0;

  /* return the structure pointer. */
  return par;
}

/* param_free(): free all allocated memory associated with a molecular
 * parameters data structure.
 *
 * arguments:
 *  @par: pointer to the data structure to free.
 */
void param_free (param_t *par) {
  /* declare required variables:
   *  @i: general-purpose loop counter.
   */
  unsigned int i;

  /* return if the structure pointer is null. */
  if (!par) return;

  /* check if the bond array is allocated. */
  if (par->n_bonds) {
    /* loop over the array entries. */
    for (i = 0; i < par->n_bonds; i++) {
      /* free the strings of the entry. */
      free(par->bonds[i].a);
      free(par->bonds[i].b);
    }

    /* free the array. */
    free(par->bonds);
    par->n_bonds = 0;
  }

  /* check if the angle array is allocated. */
  if (par->n_angles) {
    /* loop over the array entries. */
    for (i = 0; i < par->n_angles; i++) {
      /* free the strings of the entry. */
      free(par->angles[i].a);
      free(par->angles[i].b);
      free(par->angles[i].c);
    }

    /* free the array. */
    free(par->angles);
    par->n_angles = 0;
  }

  /* check if the torsion array is allocated. */
  if (par->n_torsions) {
    /* loop over the array entries. */
    for (i = 0; i < par->n_torsions; i++) {
      /* free the strings of the entry. */
      free(par->torsions[i].a);
      free(par->torsions[i].b);
      free(par->torsions[i].c);
      free(par->torsions[i].d);
    }

    /* free the array. */
    free(par->torsions);
    par->n_torsions = 0;
  }

  /* check if the improper array is allocated. */
  if (par->n_impropers) {
    /* loop over the array entries. */
    for (i = 0; i < par->n_impropers; i++) {
      /* free the strings of the entry. */
      free(par->impropers[i].a);
      free(par->impropers[i].b);
      free(par->impropers[i].c);
      free(par->impropers[i].d);
    }

    /* free the array. */
    free(par->impropers);
    par->n_impropers = 0;
  }

  /* check if the radius array is allocated. */
  if (par->n_radii) {
    /* loop over the array entries. */
    for (i = 0; i < par->n_radii; i++)
      free(par->radii[i].a);

    /* free the array. */
    free(par->radii);
    par->n_radii = 0;
  }

  /* finally, free the structure pointer. */
  free(par);
}

/* param_new_from_file(): allocate a new molecular parameters data structure
 * by reading information from a "toppar"-style file.
 *
 * arguments:
 *  @fname: input filename string.
 *
 * returns:
 *  pointer to a newly allocated and filled data structure.
 */
param_t *param_new_from_file (const char *fname, double vdw_scale) {
  /* declare required variables:
   *  @par: output structure pointer.
   *  @fh: input file handle.
   */
  param_t *par;
  FILE *fh;

  /* allocate the structure pointer. */
  par = param_new();
  if (!par)
    return NULL;

  /* store the scaling factor. */
  par->vdw_scale = vdw_scale;

  /* open the input file. */
  fh = fopen(fname, "r");

  /* check that the file was opened successfully. */
  if (!fh) {
    /* no. raise an exception and return null. */
    raise("unable to open '%s' for reading", fname);
    param_free(par);
    return NULL;
  }

  /* attempt to parse the input file. */
  if (!param_parse(fh, par)) {
    /* failure. raise an exception. */
    raise("unable to parse parameter information");

    /* clean up. */
    param_free(par);
    fclose(fh);

    /* return null. */
    return NULL;
  }

  /* close the input file. */
  fclose(fh);

  /* return the structure pointer. */
  return par;
}

