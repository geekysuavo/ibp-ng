
/* include the residue header. */
#include "resid.h"

/* residue_t: structure for mapping between one-letter codes,
 * three-letter codes, and common names of a peptide residue.
 */
typedef struct {
  /* @code1: single-letter code.
   * @code3: three-letter code.
   * @name: common name.
   */
  char code1;
  char code3[4];
  char *name;

  /* @nr: number of atoms in the sidechain "R-group".
   */
  unsigned int nr;
}
residue_t;

/* residues: mapping between identifiers of all twenty common amino acids.
 */
const residue_t residues[] = {
/*  0 */  { '\0', "   ", NULL,             0 },
/*  1 */  { 'A',  "ALA", "Alanine",        4 },
/*  2 */  { 'C',  "CYS", "Cysteine",       5 },
/*  3 */  { 'D',  "ASP", "Aspartate",      6 },
/*  4 */  { 'E',  "GLU", "Glutamate",      9 },
/*  5 */  { 'F',  "PHE", "Phenylalanine", 14 },
/*  6 */  { 'G',  "GLY", "Glycine",        1 },
/*  7 */  { 'H',  "HIS", "Histidine",     11 },
/*  8 */  { 'I',  "ILE", "Isoleucine",    13 },
/*  9 */  { 'K',  "LYS", "Lysine",        16 },
/* 10 */  { 'L',  "LEU", "Leucine",       13 },
/* 11 */  { 'M',  "MET", "Methionine",    11 },
/* 12 */  { 'N',  "ASN", "Asparagine",     8 },
/* 13 */  { 'P',  "PRO", "Proline",        9 },
/* 14 */  { 'Q',  "GLN", "Glutamine",     11 },
/* 15 */  { 'R',  "ARG", "Arginine",      18 },
/* 16 */  { 'S',  "SER", "Serine",         5 },
/* 17 */  { 'T',  "THR", "Threonine",      8 },
/* 18 */  { 'V',  "VAL", "Valine",        10 },
/* 19 */  { 'W',  "TRP", "Tryptophan",    16 },
/* 20 */  { 'Y',  "TYR", "Tyrosine",      15 },
/* 21 */  { '\0', "   ", NULL,             0 }
};

/* resid_lookup1(): lookup the index of an amino acid
 * based on its single-letter code.
 *
 * arguments:
 *  @code: single-letter code query.
 *
 * returns:
 *  index of the queried amino acid in residues[].
 */
unsigned int resid_lookup1 (char code) {
  /* declare required variables:
   *  @i: residue mapping index.
   */
  unsigned int i;

  /* search the residue array. */
  for (i = 1; residues[i].name; i++) {
    /* break the search if a match is found. */
    if (residues[i].code1 == (char) toupper((int) code))
      return i;
  }

  /* return a mismatch residue index. */
  return 0;
}

/* resid_lookup3(): lookup the index of an amino acid
 * based on its three-letter code.
 *
 * arguments:
 *  @code: three-letter code query.
 *
 * returns:
 *  index of the queried amino acid in residues[].
 */
unsigned int resid_lookup3 (const char *code) {
  /* declare required variables:
   *  @i: residue mapping index.
   *  @ucode: uppercase code string.
   */
  unsigned int i;
  char *ucode;

  /* build the uppercase code string. */
  ucode = strtoupper(code);
  if (!ucode)
    return 0;

  /* search the residue array. */
  for (i = 1; residues[i].name; i++) {
    /* break the search if a match is found. */
    if (strncmp(residues[i].code3, ucode, 3) == 0) {
      /* free the uppercase string and return. */
      free(ucode);
      return i;
    }
  }

  /* return a mismatch residue index. */
  free(ucode);
  return 0;
}

/* resid_get_code1(): get the single-letter code of a residue from its index.
 *
 * arguments:
 *  @i: residue index.
 *
 * returns:
 *  residue single-letter code, or '\0' if no match.
 */
char resid_get_code1 (unsigned int i) {
  /* return the value. */
  return (i >= 1 && i <= 20 ? residues[i].code1 : residues[0].code1);
}

/* resid_get_code3(): get the three-letter code of a residue from its index.
 *
 * arguments:
 *  @i: residue index.
 *
 * returns:
 *  residue three-letter code, or "   " if no match.
 */
const char *resid_get_code3 (unsigned int i) {
  /* return the value. */
  return (i >= 1 && i <= 20 ? residues[i].code3 : residues[0].code3);
}

/* resid_get_nr(): get the R-group atom count of a residue from its index.
 *
 * arguments:
 *  @i: residue index.
 *
 * returns:
 *  residue sidechain group atom count, or 0 if no match.
 */
unsigned int resid_get_nr (unsigned int i) {
  /* return the value. */
  return (i >= 1 && i <= 20 ? residues[i].nr : residues[0].nr);
}

