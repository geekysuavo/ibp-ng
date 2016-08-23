
/* enable verbose error reporting. */
%error-verbose
%debug

%{
/* include the assignment header. */
#include "assign.h"

/* include the parsing header. */
#include "assign-parse.h"

/* provide a link to the input file handle. */
extern FILE *assign_io_in;
peptide_t *pep;

/* flex function declarations: */
void assign_io_error (const char *msg);
void assign_io_clean (void);
int assign_io_lex (void);
%}

/* define the yylval type union. */
%union {
  assign_set_t *setval;
  int ival;
  float fval;
  char *sval;
  char cval;
}

/* define all recognized tokens. */
%token T_ASSIGN T_ID T_RESIDUE T_RESNAME T_NAME T_TYPE
%token T_OR T_AND T_NOT T_ALL T_NONE
%token T_INT T_FLOAT T_WORD
%token T_PAREN_OPEN T_PAREN_CLOSE T_EOL
%token T_UNKNOWN

/* define the types of parsed values. */
%type<setval> sel expr rule all none atomid residue resname name type not
%type<ival> T_INT
%type<fval> T_FLOAT
%type<sval> T_WORD
%type<cval> T_TYPE

%%

/* assignment file structure: a simple list of assignment statements. */
root: opteols assigns;

/* eols: list of newline tokens. */
eols: T_EOL | eols T_EOL;

/* opteols: optional list of newline tokens. */
opteols:
 | eols
 ;

/* assigns: list of assignment statements. */
assigns: assign | assigns assign ;

/* assign: either a distance or dihedral assignment statement. */
assign: distance | dihedral ;

/* distance: assigns a distance restraint between two atoms. */
distance: T_ASSIGN sel sel T_FLOAT T_FLOAT T_FLOAT eols {
  /* add the parsed distance restraint to the peptide structure. */
  int ret = assign_set_distance(pep, $2, $3, $4, $5, $6);

  /* free the atom selectors. */
  assign_free($2);
  assign_free($3);

  /* check for errors. */
  if (!ret) YYERROR;
};

/* dihedral: assigns a dihedral restraint among four atoms. */
dihedral: T_ASSIGN sel sel sel sel T_FLOAT T_FLOAT T_FLOAT T_INT eols {
  /* add the parsed dihedral restraint to the peptide structure. */
  int ret = assign_set_dihedral(pep, $2, $3, $4, $5, $7, $8);

  /* free the atom selectors. */
  assign_free($2);
  assign_free($3);
  assign_free($4);
  assign_free($5);

  /* check for errors. */
  if (!ret) YYERROR;
};

/* sel: parenthetical atom selector rule used in assignments. */
sel: T_PAREN_OPEN expr T_PAREN_CLOSE opteols { $$ = $2; };

/* expr: atom selector expression rule. */
expr: rule /* use the default ($$ = $1) here. */
 | expr T_OR rule {
  /* perform a set "or" operation between two assignments. */
  $$ = assign_or($1, $3);

  /* free the input assignments. */
  assign_free($1);
  assign_free($3);

  /* check for errors. */
  if (!($$)) YYERROR;
}
 | expr T_AND rule {
  /* perform a set "and" operation between two assignments. */
  $$ = assign_and($1, $3);

  /* free the input assignments. */
  assign_free($1);
  assign_free($3);

  /* check for errors. */
  if (!($$)) YYERROR;
};

/* rule: simple assignment rule for building atom selections. */
rule: all
 | none
 | atomid
 | residue
 | resname
 | name
 | type
 | sel
 | not
 ;

/* all: initialize a set of all possible atoms. */
all: T_ALL {
  /* create the assignment set. */
  $$ = assign_all(pep);
  if (!($$))
    YYERROR;
};

/* none: initialize an empty set of atoms. */
none: T_NONE {
  /* create the assignment set. */
  $$ = assign_none(pep);
  if (!($$))
    YYERROR;
};

/* not: negate a set with the set of all possible atoms. */
not: T_NOT sel {
  /* perform the negation and free the input set. */
  $$ = assign_not($2);
  assign_free($2);

  /* check for errors. */
  if (!($$)) YYERROR;
};

/* atomid: initialize a set based on internal atom index. */
atomid: T_ID T_INT {
  /* initialize the assignment set. */
  $$ = assign_atomid(pep, $2);
  if (!($$))
    YYERROR;
};

/* residue: initialize a set based on residue index. */
residue: T_RESIDUE T_INT {
  /* initialize the assignment set. */
  $$ = assign_resid(pep, $2);
  if (!($$))
    YYERROR;
};

/* resname: initialize a set based on residue name. */
resname: T_RESNAME T_WORD {
  /* initialize the assignment set and free the input word. */
  $$ = assign_resname(pep, $2);
  free($2);

  /* check for errors. */
  if (!($$))
    YYERROR;
};

/* name: initialize a set based on atom name. */
name: T_NAME T_WORD {
  /* initialize the assignment set and free the input word. */
  $$ = assign_name(pep, $2);
  free($2);

  /* check for errors. */
  if (!($$)) YYERROR;
};

/* type: initialize a set of all atoms having a specific type. */
type: T_TYPE {
  /* create the assignment set. */
  $$ = assign_type(pep, $1);
  if (!($$))
    YYERROR;
};

%%

/* assign_parse(): extract distance and dihedral restraint information
 * from XPLOR-format restraint files.
 *
 * arguments:
 *  @fh: input file handle to read from.
 *  @P: pointer to the peptide structure to parse assignments for.
 *
 * returns:
 *  integer indicating whether (1) or not (0) parsing succeeded.
 *
 * warning:
 *  this function is absolutely NOT re-entrant.
 */
int assign_parse (FILE *fh, peptide_t *P) {
  /* set the input file handle and peptide structure pointer. */
  assign_io_in = fh;
  pep = P;

  /* attempt to parse the input file. */
  if (assign_io_parse()) {
    /* free the scanner buffer and return failure. */
    assign_io_clean();
    return 0;
  }

  /* free the scanner buffer. */
  assign_io_clean();

  /* return success. */
  errno = 0;
  return 1;
}

/* assign_io_error(): error reporting interface for parsing restraint files.
 *
 * arguments:
 *  @msg: error message string.
 */
void assign_io_error (const char *msg) {
  /* raise an exception, to be handled by assign_parse(). */
  raise(msg);
}

