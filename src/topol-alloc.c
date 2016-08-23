
/* include the molecular topology header. */
#include "topol.h"

/* function declarations from parsing code: */
int topol_parse (FILE *fh, topol_t *top);

/* topol_new(): allocate a new molecular topology data structure.
 *
 * returns:
 *  pointer to a newly allocated empty data structure.
 */
topol_t *topol_new (void) {
  /* declare required variables:
   *  @top: output structure pointer.
   */
  topol_t *top;

  /* allocate a new structure pointer. */
  top = (topol_t*) malloc(sizeof(topol_t));

  /* check if allocation failed. */
  if (!top) {
    /* raise an exception and return null. */
    raise("unable to allocate topology structure pointer");
    return NULL;
  }

  /* initialize the residue array. */
  top->res = NULL;
  top->n_res = 0;

  /* initialize the mass array. */
  top->mass = NULL;
  top->n_mass = 0;

  /* initialize the autogeneration options. */
  top->auto_angles = 0;
  top->auto_dihedrals = 0;

  /* return the structure pointer. */
  return top;
}

/* topol_free(): free all allocated memory associated with a molecular
 * topology data structure.
 *
 * arguments:
 *  @top: pointer to the data structure to free.
 */
void topol_free (topol_t *top) {
  /* declare required variables:
   *  @i, @j: general-purpose loop counters.
   */
  unsigned int i, j;

  /* return if the structure pointer is null. */
  if (!top) return;

  /* check if the residue array is allocated. */
  if (top->n_res) {
    /* loop over the residue array entries. */
    for (i = 0; i < top->n_res; i++) {
      /* free the residue name. */
      if (top->res[i].name)
        free(top->res[i].name);

      /* check if the atom array is allocated. */
      if (top->res[i].n_atoms) {
        /* loop over the atom array entries. */
        for (j = 0; j < top->res[i].n_atoms; j++) {
          /* free the atom strings. */
          free(top->res[i].atoms[j].name);
          free(top->res[i].atoms[j].type);
        }

        /* free the array. */
        free(top->res[i].atoms);
      }

      /* check if the bond array is allocated. */
      if (top->res[i].n_bonds) {
        /* loop over the bond array entries. */
        for (j = 0; j < top->res[i].n_bonds; j++) {
          /* free the bond strings. */
          free(top->res[i].bonds[j].atoms[0]);
          free(top->res[i].bonds[j].atoms[1]);
        }

        /* free the array. */
        free(top->res[i].bonds);
      }

      /* check if the angle array is allocated. */
      if (top->res[i].n_angles) {
        /* loop over the angle array entries. */
        for (j = 0; j < top->res[i].n_angles; j++) {
          /* free the angle strings. */
          free(top->res[i].angles[j].atoms[0]);
          free(top->res[i].angles[j].atoms[1]);
          free(top->res[i].angles[j].atoms[2]);
        }

        /* free the array. */
        free(top->res[i].angles);
      }

      /* check if the torsion array is allocated. */
      if (top->res[i].n_torsions) {
        /* loop over the torsion array entries. */
        for (j = 0; j < top->res[i].n_torsions; j++) {
          /* free the torsion strings. */
          free(top->res[i].torsions[j].atoms[0]);
          free(top->res[i].torsions[j].atoms[1]);
          free(top->res[i].torsions[j].atoms[2]);
          free(top->res[i].torsions[j].atoms[3]);
        }

        /* free the array. */
        free(top->res[i].torsions);
      }

      /* check if the improper array is allocated. */
      if (top->res[i].n_impropers) {
        /* loop over the improper array entries. */
        for (j = 0; j < top->res[i].n_impropers; j++) {
          /* free the improper strings. */
          free(top->res[i].impropers[j].atoms[0]);
          free(top->res[i].impropers[j].atoms[1]);
          free(top->res[i].impropers[j].atoms[2]);
          free(top->res[i].impropers[j].atoms[3]);
        }

        /* free the array. */
        free(top->res[i].impropers);
      }

      /* reset the counters. */
      top->res[i].n_atoms = 0;
      top->res[i].n_bonds = 0;
      top->res[i].n_angles = 0;
      top->res[i].n_torsions = 0;
      top->res[i].n_impropers = 0;
    }

    /* free the residue array. */
    free(top->res);
    top->n_res = 0;
  }

  /* check if the mass array is allocated. */
  if (top->n_mass) {
    /* loop over the mass array entries. */
    for (i = 0; i < top->n_mass; i++) {
      /* free the mass entry string. */
      if (top->mass[i].type)
        free(top->mass[i].type);
    }

    /* free the mass array. */
    free(top->mass);
    top->n_mass = 0;
  }

  /* finally, free the structure pointer. */
  free(top);
}

/* topol_new_from_file(): allocate a new molecular topology data structure
 * by reading information from a "toppar"-style file.
 *
 * arguments:
 *  @fname: input parameter filename string.
 *
 * returns:
 *  pointer to a newly allocated and filled data structure.
 */
topol_t *topol_new_from_file (const char *fname) {
  /* declare required variables:
   *  @top: output structure pointer.
   *  @fh: input file handle.
   */
  topol_t *top;
  FILE *fh;

  /* allocate the structure pointer. */
  top = topol_new();
  if (!top)
    return NULL;

  /* open the input file. */
  fh = fopen(fname, "r");

  /* check that the file was opened successfully. */
  if (!fh) {
    /* no. raise an exception and return null. */
    raise("unable to open '%s' for reading", fname);
    topol_free(top);
    return NULL;
  }

  /* attempt to parse the input file. */
  if (!topol_parse(fh, top)) {
    /* failure. raise an exception. */
    raise("unable to parse topology information");

    /* clean up. */
    topol_free(top);
    fclose(fh);

    /* return null. */
    return NULL;
  }

  /* close the input file. */
  fclose(fh);

  /* return the structure pointer. */
  return top;
}

