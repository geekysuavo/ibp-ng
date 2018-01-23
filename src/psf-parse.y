
/* enable verbose error reporting. */
%error-verbose
%debug

%{
/* include the parsing and peptide headers. */
#include "psf-parse.h"
#include "peptide.h"

/* provide a link to the input file handle. */
extern FILE *psf_io_in;

/* provide variables for parsing sequence data into peptides:
 *  @P: output peptide structure pointer.
 */
peptide_t *P;

/* flex function declarations: */
void psf_io_error (const char *msg);
void psf_io_clean (void);
int psf_io_lex (void);
%}

/* define the yylval type union. */
%union {
  int ival;
  float fval;
  char *sval;
}

/* define all recognized tokens. */
%token T_BEGIN T_REMARK T_INT T_FLOAT T_WORD T_ATOM T_OTHER
%token T_UNKNOWN

/* define the types of parsed values. */
%type<ival> T_INT
%type<fval> T_FLOAT
%type<sval> T_WORD

%%

/* general psf format: list of blocks, either parsed as atom blocks
 * or ignored as "other" blocks.
 */
psf: T_BEGIN blocks ;
blocks: blk | blocks blk ;
blk: blk_atom | blk_other ;

/* atom block: atom header and list of atom entries. */
blk_atom: T_ATOM atoms ;

/* "other" block: any header followed by remarks, integers, or nothing. */
blk_other: T_OTHER | T_OTHER remarks | T_OTHER ints ;

/* atoms: list of atom entries. */
atoms: atom | atoms atom ;

/* atom: token group identifying a single atom entry. */
atom: T_INT T_INT T_WORD T_WORD T_WORD T_FLOAT T_FLOAT T_INT {
  /* add the residue into the peptide sequence. */
  if ($1 > 0 && $2 >= (int) P->n_res && !peptide_add_residue(P, $3)) {
    /* free the allocated strings. */
    free($3);
    free($4);
    free($5);

    /* return failure. */
    YYERROR;
  }

  /* free the allocated strings. */
  free($3);
  free($4);
  free($5);
};

/* ints: list of one or more integers. */
ints: T_INT | ints T_INT ;

/* remarks: list of remark entries. */
remarks: T_REMARK | remarks T_REMARK ;

%%

/* psf_parse(): extract peptide sequence information from an XPLOR PSF file.
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
int psf_parse (FILE *fh, peptide_t *pep) {
  /* set the input file handle and output structure pointer. */
  psf_io_in = fh;
  P = pep;

  /* attempt to parse the input file. */
  if (psf_io_parse()) {
    /* free the scanner buffer and return failure. */
    psf_io_clean();
    return 0;
  }

  /* free the scanner buffer. */
  psf_io_clean();

  /* return success. */
  errno = 0;
  return 1;
}

/* psf_io_error(): error reporting interface for parsing XPLOR PSF files.
 *
 * arguments:
 *  @msg: error message string.
 */
void psf_io_error (const char *msg) {
  /* raise an exception, to be handled by psf_parse(). */
  raise(msg);
}

