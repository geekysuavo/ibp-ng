
/* enable verbose error reporting. */
%error-verbose
%debug

%{
/* include the parsing and peptide headers. */
#include "fasta-parse.h"
#include "peptide.h"

/* provide a link to the input file handle. */
extern FILE *fasta_io_in;

/* provide variables for parsing sequence data into peptides:
 *  @fasta_idx: index of the currently parsed sequence.
 *  @P: output peptide structure pointer.
 */
unsigned int fasta_idx, idx;
peptide_t *P;

/* flex function declarations: */
void fasta_io_error (const char *msg);
void fasta_io_clean (void);
int fasta_io_lex (void);

/* fasta function declarations: */
const char *fasta_lookup (char code);

/* struct fasta_lut_entry: structure for mapping between single-letter
 * and three-letter codes used by fasta-format sequence files.
 */
struct fasta_lut_entry {
  /* @code1: single-letter code.
   * @code3: three-letter code.
   */
  char code1;
  char code3[4];
};

/* fasta_lut: mapping between all single-letter and three-letter amino
 * acid codes accepted by canonical fasta sequences.
 */
const struct fasta_lut_entry fasta_lut[] = {
/*  0 */  { '\0', "   " },
/*  1 */  { 'A',  "ALA" },
/*  2 */  { 'C',  "CYS" },
/*  3 */  { 'D',  "ASP" },
/*  4 */  { 'E',  "GLU" },
/*  5 */  { 'F',  "PHE" },
/*  6 */  { 'G',  "GLY" },
/*  7 */  { 'H',  "HIS" },
/*  8 */  { 'I',  "ILE" },
/*  9 */  { 'K',  "LYS" },
/* 10 */  { 'L',  "LEU" },
/* 11 */  { 'M',  "MET" },
/* 12 */  { 'N',  "ASN" },
/* 13 */  { 'P',  "PRO" },
/* 14 */  { 'Q',  "GLN" },
/* 15 */  { 'R',  "ARG" },
/* 16 */  { 'S',  "SER" },
/* 17 */  { 'T',  "THR" },
/* 18 */  { 'V',  "VAL" },
/* 19 */  { 'W',  "TRP" },
/* 20 */  { 'Y',  "TYR" },
/* 21 */  { '\0', "   " }};
%}

/* define the yylval type union. */
%union {
  char cval;
}

/* define all recognized tokens. */
%token T_EOL T_BEGIN T_RESIDUE
%token T_UNKNOWN

/* define the types of parsed values. */
%type<cval> T_RESIDUE

%%

/* global listing of all entries. */
entries: entry
 | entries entry
 ;

/* entry instance tokens. */
entry: begin sequence optional_newlines {
  /* increment the sequence index. */
  fasta_idx++;
};

/* sequence beginning instance tokens. */
begin: T_BEGIN ;

sequence: sequence residues T_EOL
 | residues T_EOL
 ;

/* sequence residue list tokens. */
residues: residues residue
 | residue
 ;

/* single residue instance tokens. */
residue: T_RESIDUE {
  /* check if we are to parse the current sequence. */
  if (fasta_idx == idx) {
    /* search for the residue three-letter code. */
    const char *resname = fasta_lookup($1);
    if (!resname)
      YYERROR;

    /* add the residue into the peptide sequence. */
    if (!peptide_add_residue(P, resname))
      YYERROR;
  }
};

/* list of newline characters. */
newlines: newlines T_EOL
 | T_EOL
 ;

/* optional rule for matching spacing between sequences. */
optional_newlines:
 | newlines
 ;

%%

/* fasta_parse(): extract peptide sequence information from a FASTA file.
 *
 * arguments:
 *  @fh: input file handle to read from.
 *  @sidx: sequence index in the input file.
 *  @pep: pointer to the parameter data structure to modify.
 *
 * returns:
 *  integer indicating whether (1) or not (0) parsing succeeded.
 *
 * warning:
 *  this function is absolutely NOT re-entrant.
 */
int fasta_parse (FILE *fh, unsigned int sidx, peptide_t *pep) {
  /* set the input file handle and output structure pointer. */
  fasta_io_in = fh;
  P = pep;

  /* initialize the sequence index. */
  fasta_idx = 1;
  idx = sidx;

  /* check the sequence index. */
  if (idx < 1)
    throw("sequence index %u out of bounds [1,inf)", sidx);

  /* attempt to parse the input file. */
  if (fasta_io_parse()) {
    /* free the scanner buffer and return failure. */
    fasta_io_clean();
    return 0;
  }

  /* free the scanner buffer. */
  fasta_io_clean();

  /* return success. */
  errno = 0;
  return 1;
}

/* fasta_io_error(): error reporting interface for parsing FASTA files.
 *
 * arguments:
 *  @msg: error message string.
 */
void fasta_io_error (const char *msg) {
  /* raise an exception, to be handled by fasta_parse(). */
  raise(msg);
}

/* fasta_lookup(): lookup the three-letter code of an amino acid
 * from its single-letter code character.
 *
 * arguments:
 *  @code: single-letter code query.
 *
 * returns:
 *  string containing the three-letter code, or null on mismatch.
 */
const char *fasta_lookup (char code) {
  /* search the lookup table. */
  for (unsigned int i = 1; fasta_lut[i].code1; i++) {
    /* break the search if a match is found. */
    if (fasta_lut[i].code1 == (char) toupper((int) code))
      return fasta_lut[i].code3;
  }

  /* no match found. */
  return NULL;
}

