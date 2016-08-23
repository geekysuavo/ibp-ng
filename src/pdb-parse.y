
/* enable verbose error reporting. */
%error-verbose
%debug

%{
/* include the parsing and peptide headers. */
#include "pdb-parse.h"
#include "peptide.h"

/* provide a link to the input file handle. */
extern FILE *pdb_io_in;

/* provide variables for parsing sequence data into peptides:
 *  @P: output peptide structure pointer.
 */
char pdb_chain;
peptide_t *P;

/* flex function declarations: */
void pdb_io_error (const char *msg);
void pdb_io_clean (void);
int pdb_io_lex (void);
%}

/* define the yylval type union. */
%union {
  int ival;
  char *sval;
}

/* define all recognized tokens. */
%token T_SEQRES T_INT T_WORD
%token T_EOL

/* define the types of parsed values. */
%type<ival> T_INT
%type<sval> T_WORD words

%%

/* pdb sub-format: list of SEQRES lines, flanked by a bunch of other junk. */
pdb: lines ;
lines: line | lines line ;
line: others T_EOL | seqres T_EOL ;

/* others: list of anything not matching a SEQRES rule. */
others: other | others other ;

/* other: either an integer or a word. simple enough. */
other: T_INT
 | T_WORD { free($1); }
 ;

/* seqres: a complete line of SEQRES data. */
seqres: T_SEQRES T_INT T_WORD T_INT words {
  /* declare required variables:
   *  @i: string offset.
   */
  unsigned int i;

  /* check if the current chain is a match. */
  if (($3)[0] == pdb_chain) {
    /* loop for each residue in the sequence string. */
    for (i = 0; i < strlen($5); i += 4) {
      /* attempt to add the residue to the peptide sequence. */
      if (!peptide_add_residue3(P, $5 + i)) {
        /* clean up. */
        free($3);
        free($5);

        /* return failure. */
        YYERROR;
      }
    }
  }

  /* free the allocated strings. */
  free($3);
  free($5);
};

/* words: list of words to be joined together. */
words: T_WORD { $$ = $1; } | words T_WORD {
  /* compute the length of the output string. */
  int n = strlen($1) + strlen($2) + 2;

  /* allocate a new string to hold the concatenation. */
  $$ = malloc(n * sizeof(char));
  if (!($$))
    YYERROR;

  /* concatenate the two strings. */
  strcpy($$, "");
  strcat($$, $1);
  strcat($$, " ");
  strcat($$, $2);

  /* free the input strings. */
  free($1);
  free($2);
};

%%

/* pdb_parse(): extract peptide sequence information from an RCSB PDB file.
 *
 * arguments:
 *  @fh: input file handle to read from.
 *  @chain: chain character of the desired sequence.
 *  @pep: pointer to the parameter data structure to modify.
 *
 * returns:
 *  integer indicating whether (1) or not (0) parsing succeeded.
 *
 * warning:
 *  this function is absolutely NOT re-entrant.
 */
int pdb_parse (FILE *fh, char chain, peptide_t *pep) {
  /* set the input file handle and output structure pointer. */
  pdb_io_in = fh;
  P = pep;

  /* store the chain character. */
  pdb_chain = chain;

  /* attempt to parse the input file. */
  if (pdb_io_parse()) {
    /* free the scanner buffer and return failure. */
    pdb_io_clean();
    return 0;
  }

  /* free the scanner buffer. */
  pdb_io_clean();

  /* return success. */
  errno = 0;
  return 1;
}

/* pdb_io_error(): error reporting interface for parsing RCSB PDB files.
 *
 * arguments:
 *  @msg: error message string.
 */
void pdb_io_error (const char *msg) {
  /* raise an exception, to be handled by pdb_parse(). */
  raise(msg);
}

