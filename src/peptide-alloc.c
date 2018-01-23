
/* include the peptide headers. */
#include "peptide.h"
#include "peptide-alloc.h"

/* peptide_new(): allocate a new empty peptide data structure.
 *
 * returns:
 *  pointer to an allocated and initialized peptide_t structure. the pointer
 *  must be freed after the use of seq_free().
 */
peptide_t *peptide_new (void) {
  /* declare required variables:
   *  @P: pointer to a newly allocated structure.
   */
  peptide_t *P;

  /* allocate the peptide structure pointer. */
  P = (peptide_t*) malloc(sizeof(peptide_t));

  /* check if the pointer was allocated successfully. */
  if (!P) {
    /* raise an exception and return null. */
    raise("unable to allocate peptide structure pointer");
    return NULL;
  }

  /* initialize the residue array. */
  P->res = NULL;
  P->n_res = 0;

  /* initialize the sidechain array. */
  P->sc = NULL;
  P->n_sc = 0;

  /* initialize the atom array. */
  P->atoms = NULL;
  P->n_atoms = 0;

  /* initialize the bond array. */
  P->bonds = NULL;
  P->n_bonds = 0;

  /* initialize the angle array. */
  P->angles = NULL;
  P->n_angles = 0;

  /* initialize the torsion array. */
  P->torsions = NULL;
  P->n_torsions = 0;

  /* initialize the improper array. */
  P->impropers = NULL;
  P->n_impropers = 0;

  /* return the newly allocated peptide. */
  return P;
}

/* peptide_free(): free all allocated memory associated with a peptide.
 *
 * arguments:
 *  @P: pointer to the peptide data structure to free.
 */
void peptide_free (peptide_t *P) {
  /* declare required variables:
   *  @i: atom loop counter.
   */
  unsigned int i;

  /* return if the structure pointer is null. */
  if (!P) return;

  /* free the elements of the residue names array. */
  for (i = 0; i < P->n_res; i++)
    free(P->res[i]);

  /* free the residue names array. */
  if (P->n_res)
    free(P->res);

  /* free the sidechain array. */
  if (P->n_sc)
    free(P->sc);

  /* loop over the atoms in the peptide. */
  for (i = 0; i < P->n_atoms; i++) {
    /* free the atom strings. */
    free(P->atoms[i].name);
    free(P->atoms[i].type);
  }

  /* free the atom array. */
  if (P->n_atoms)
    free(P->atoms);

  /* free the bond array. */
  if (P->n_bonds)
    free(P->bonds);

  /* free the angle array. */
  if (P->n_angles)
    free(P->angles);

  /* free the torsion array. */
  if (P->n_torsions)
    free(P->torsions);

  /* free the improper array. */
  if (P->n_impropers)
    free(P->impropers);

  /* reset the array counts. */
  P->n_res = 0;
  P->n_sc = 0;
  P->n_atoms = 0;
  P->n_bonds = 0;
  P->n_angles = 0;
  P->n_torsions = 0;
  P->n_impropers = 0;

  /* finally, free the structure pointer. */
  free(P);
}

/* peptide_new_from_file(): allocate and fill a new peptide data structure
 * from a file, which is automatically determined to be either PSF, FASTA,
 * or PDB.
 *
 * arguments:
 *  @fname: the filename of the input file.
 *  @opts: extra string containing special options.
 *
 * returns:
 *  pointer to an allocated, initialized and filled structure, or NULL on
 *  failure.
 */
peptide_t *peptide_new_from_file (const char *fname, const char *opts) {
  /* declare required variables:
   *  @buf: buffer for file type determination.
   *  @P: output peptide structure pointer.
   *  @fh: file handle for buffer reading.
   */
  peptide_t *P;
  char buf[32];
  FILE *fh;
  int ret;

  /* initialize the buffer string. */
  strcpy(buf, "");

  /* attempt to open the input file. */
  fh = fopen(fname, "r");

  /* check if file opening failed. */
  if (!fh) {
    /* raise an exception and return null. */
    raise("unable to open '%s' for reading", fname);
    return NULL;
  }

  /* read a small chunk of character data from the file. */
  if (!fgets(buf, 32, fh)) {
    /* read failure. raise an exception and return null. */
    raise("unable to read header of '%s'", fname);
    fclose(fh);
    return NULL;
  }

  /* rewind the input file. */
  if (fseek(fh, 0, SEEK_SET)) {
    /* seek failure. raise an exception and return null. */
    raise("unable to rewind '%s'", fname);
    fclose(fh);
    return NULL;
  }

  /* allocate the structure pointer. */
  P = peptide_new();

  /* check if the structure pointer was allocated. */
  if (!P) {
    /* nope. close the input file and return null. */
    fclose(fh);
    return NULL;
  }

  /* guess the file type from the header chunk. */
  if (strncmp(buf, "PSF", 3) == 0) {
    /* interpret the file using the XPLOR PSF format. */
    ret = psf_parse(fh, P);
  }
  else if (strncmp(buf, "data_cns_mtf", 12) == 0) {
    /* interpret the file using the CNS format. */
    ret = cns_parse(fh, P);
  }
  else if (strncmp(buf, ">", 1) == 0 ||
           strncmp(buf, ";", 1) == 0) {
    /* interpret the file using the FASTA format. */
    if (opts && strlen(opts) > 0)
      ret = fasta_parse(fh, atoi(opts), P);
    else
      ret = fasta_parse(fh, 1, P);
  }
  else if (strncmp(buf, "HEADER", 6) == 0 ||
           strncmp(buf, "REMARK", 6) == 0) {
    /* interpret the file using the RCSB PDB format. */
    if (opts && strlen(opts) > 0)
      ret = pdb_parse(fh, opts[0], P);
    else
      ret = pdb_parse(fh, 'A', P);
  }
  else {
    /* unknown file type. return null. */
    raise("file '%s' has an unsupported type", fname);
    ret = 0;
  }

  /* close the input file. */
  fclose(fh);

  /* check for parsing errors. */
  if (!ret) {
    /* raise an exception and return null. */
    raise("failed to parse '%s'", fname);
    peptide_free(P);
    return NULL;
  }

  /* check that residues were parsed. */
  if (P->n_res == 0) {
    /* raise an exception and return null. */
    raise("unable to read sequence from '%s'", fname);
    peptide_free(P);
    return NULL;
  }

  /* return the completed peptide structure pointer. */
  return P;
}

