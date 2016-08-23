
/* enable verbose error reporting. */
%error-verbose
%debug

%{
/* include the repetition ordering headers. */
#include "reorder-parse.h"
#include "reorder.h"

/* provide a link to the input file handle and output structure pointer.
 */
extern FILE *reorder_io_in;
reorder_t *RO;

/* flex function declarations: */
void reorder_io_error (const char *msg);
void reorder_io_clean (void);
int reorder_io_lex (void);
%}

/* define the yylval type union. */
%union {
  char *sval;
  int ival;
}

/* define all recognized tokens. */
%token T_REORDER T_END T_OPT T_PREV T_CURR T_NEXT T_WORD
%token T_UNKNOWN

/* define the types of parsed values. */
%type<ival> reorder_modifier reorder_offset
%type<sval> T_WORD

%%

/* global listing of all reorder information. */
reorders: reorder
 | reorders reorder
 ;

/* single reorder block. */
reorder: reorder_begin reorder_entries T_END ;

/* reorder block initial tokens. */
reorder_begin: T_REORDER T_WORD {
  /* add a new residue into the reorder. */
  if (!reorder_add_residue(RO, $2)) {
    /* free the allocated word and return failure. */
    free($2);
    YYERROR;
  }
};

/* list of acceptable reorder entries. */
reorder_entries: reorder_entry
 | reorder_entries reorder_entry
 ;

/* reorder entry tokens. */
reorder_entry: reorder_modifier reorder_offset T_WORD {
  /* add a new atom to the last residue of the reorder. */
  if (!reorder_add_atom(RO, $3, $2, $1)) {
    /* free the allocated word and return failure. */
    free($3);
    YYERROR; 
  }
};

/* reorder modifier tokens. */
reorder_modifier: { $$ = 0; }
 | T_OPT          { $$ = 1; }
 ;

/* reorder offset tokens. */
reorder_offset: { $$ =  0; }
 | T_CURR       { $$ =  0; }
 | T_PREV       { $$ = -1; }
 | T_NEXT       { $$ =  1; }
 ;

%%

/* reorder_parse(): extract repetition ordering information from a file.
 *
 * arguments:
 *  @fh: input file handle to read from.
 *  @ord: pointer to the reorder data structure to modify.
 *
 * returns:
 *  integer indicating whether (1) or not (0) parsing succeeded.
 *
 * warning:
 *  this function is absolutely NOT re-entrant.
 */
int reorder_parse (FILE *fh, reorder_t *ord) {
  /* set the input file handle and output structure pointer. */
  reorder_io_in = fh;
  RO = ord;

  /* attempt to parse the input file. */
  if (reorder_io_parse()) {
    /* free the scanner buffer and return failure. */
    reorder_io_clean();
    return 0;
  }

  /* free the scanner buffer. */
  reorder_io_clean();

  /* return success. */
  errno = 0;
  return 1;
}

/* reorder_io_error(): error reporting interface for parsing repetition
 * ordering files.
 *
 * arguments:
 *  @msg: error message string.
 */
void reorder_io_error (const char *msg) {
  /* raise an exception, to be handled by reorder_parse(). */
  raise(msg);
}

