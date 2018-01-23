
/* enable verbose error reporting. */
%error-verbose
%debug

%{
/* include the parsing and peptide headers. */
#include "cns-parse.h"
#include "peptide.h"

/* provide a link to the input file handle. */
extern FILE *cns_io_in;

/* provide variables for parsing sequence data into peptides:
 *  @P: output peptide structure pointer.
 */
peptide_t *P;

/* flex function declarations: */
void cns_io_error (const char *msg);
void cns_io_clean (void);
int cns_io_lex (void);
%}

/* define the yylval type union. */
%union {
  int ival;
  float fval;
  char *sval;
}

/* define all recognized tokens. */
%token T_BEGIN T_LOOP T_INT T_FLOAT T_WORD T_ATOM T_OTHER T_QUOT
%token T_UNKNOWN

/* define the types of parsed values. */
%type<ival> T_INT
%type<fval> T_FLOAT
%type<sval> T_WORD

%%

/* general cns format: list of blocks, either parsed as atom blocks
 * or ignored as "other" blocks.
 */
cns: T_BEGIN T_OTHER blocks ;
blocks: blk | blocks blk ;
blk: T_LOOP blk_atom | T_LOOP blk_other ;

/* atom block: atom field header and list of atom entries. */
blk_atom: fields atoms ;

/* "other" block: any field header followed by integers. */
blk_other: others ints ;

/* atom field header: one or more atom field lines. */
fields: T_ATOM | fields T_ATOM ;

/* atoms: list of atom entries. */
atoms: atom | atoms atom ;

/* atom: token group identifying a single atom entry. */
atom: T_INT T_WORD T_WORD T_WORD T_WORD T_WORD T_FLOAT T_FLOAT {
  /* parse the third field as an integer. */
  int atom_ires = atoi($3 + 1);

  /* add the residue into the peptide sequence. */
  if ($1 > 0 && atom_ires >= (int) P->n_res &&
      !peptide_add_residue(P, $4)) {
    /* free the allocated strings. */
    free($2);
    free($3);
    free($4);
    free($5);
    free($6);

    /* return failure. */
    YYERROR;
  }

  /* free the allocated strings. */
  free($2);
  free($3);
  free($4);
  free($5);
  free($6);
};

/* others: list of "other" field headers. */
others: T_OTHER | others T_OTHER ;

/* ints: list of one or more integers. */
ints: T_INT | ints T_INT ;

%%

/* cns_parse(): extract peptide sequence information from a CNS PSF file.
 *
 * arguments:
 *  @fh: input file handle to read from.
 *  @pep: pointer to the parameter data structure to modify.
 *
 * returns:
 *  integer indicating whether (1) or not (0) parsing succeeded.
 *
 * warning:
 *  this function is absolutely NOT re-entrant.
 */
int cns_parse (FILE *fh, peptide_t *pep) {
  /* set the input file handle and output structure pointer. */
  cns_io_in = fh;
  P = pep;

  /* attempt to parse the input file. */
  if (cns_io_parse()) {
    /* free the scanner buffer and return failure. */
    cns_io_clean();
    return 0;
  }

  /* free the scanner buffer. */
  cns_io_clean();

  /* return success. */
  errno = 0;
  return 1;
}

/* cns_io_error(): error reporting interface for parsing CNS PSF files.
 *
 * arguments:
 *  @msg: error message string.
 */
void cns_io_error (const char *msg) {
  /* raise an exception, to be handled by cns_parse(). */
  raise(msg);
}

