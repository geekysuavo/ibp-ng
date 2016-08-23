
/* enable verbose error reporting. */
%error-verbose
%debug

%{
/* include the molecular topology header. */
#include "topol.h"
#include "topol-add.h"
#include "topol-auto.h"

/* include the parsing header. */
#include "topol-parse.h"

/* provide a link to the input file handle and output structure pointer.
 */
extern FILE *topol_io_in;
topol_t *T;

/* flex function declarations: */
void topol_io_error (const char *msg);
void topol_io_clean (void);
int topol_io_lex (void);
%}

/* define the yylval type union. */
%union {
  topol_mode_t tval;
  unsigned int bval;
  int ival;
  float fval;
  char *sval;
}

/* define all recognized tokens. */
%token T_RESIDUE T_PATCH T_GROUP T_ATOM T_BOND T_ANGLE T_DIHEDRAL T_IMPROPER
%token T_AUTO T_TYPE T_CHARGE T_MASS T_ADD T_DELETE T_MODIFY T_END
%token T_BOOL T_WORD T_NUM T_PM T_EQ
%token T_UNKNOWN

/* define the types of parsed values. */
%type<tval> mode
%type<bval> T_BOOL
%type<ival> residue_token off T_PM
%type<sval> type id T_WORD
%type<fval> charge T_NUM

%%

/* global listing of all topology information. */
root: entries ;

/* list of top-level topology statements. */
entries: entry
 | entries entry
 ;

/* acceptable top-level topology statements. */
entry: residue | mass | autogen;

/* single residue topology. */
residue: residue_begin residue_entries T_END {
  /* finish the last residue topology information. */
  if (!topol_autogen(T))
    YYERROR;
};

/* residue topology initial tokens. */
residue_begin: residue_token T_WORD {
  /* add a new residue to the topology structure. */
  int ret = topol_add_residue(T, $2, $1);

  /* free the allocated strings. */
  free($2);

  /* check for errors. */
  if (!ret) YYERROR;
};

/* tokens for initializing new residues. */
residue_token
 : T_RESIDUE { $$ = 0; }
 | T_PATCH   { $$ = 1; }
 ;

/* list of acceptable residue entries. */
residue_entries: residue_entry
 | residue_entries residue_entry
 ;

/* residue entry types. */
residue_entry: group
 | atom
 | bond
 | angle
 | dihedral
 | improper
 ;

/* group token. */
group: T_GROUP;

/* atom instance tokens:
$$    $1   $2     $3  $4 $5   $6     $7 */
atom: mode T_ATOM off id type charge T_END {
  /* store the atom information in the last residue topology. */
  int ret = topol_add_atom(T, NULL, $4, $5, $6, $1, $3);

  /* free the allocated strings. */
  free($4);
  free($5);

  /* check for errors. */
  if (!ret) YYERROR;
};

/* bond instance tokens:
$$    $1   $2     $3  $4 $5  $6 */
bond: mode T_BOND off id off id {
  /* store the bond information in the last residue topology. */
  int ret = topol_add_bond(T, NULL, $4, $3, $6, $5, $1);

  /* free the allocated strings. */
  free($4);
  free($6);

  /* check for errors. */
  if (!ret) YYERROR;
};

/* angle instance tokens:
$$     $1   $2      $3  $4 $5  $6 $7  $8 */
angle: mode T_ANGLE off id off id off id {
  /* store the angle information in the last residue topology. */
  int ret = topol_add_angle(T, NULL, $4, $3, $6, $5, $8, $7, $1);

  /* free the allocated strings. */
  free($4);
  free($6);
  free($8);

  /* check for errors. */
  if (!ret) YYERROR;
};

/* dihedral instance tokens:
$$        $1   $2         $3  $4 $5  $6 $7  $8 $9  $10 */
dihedral: mode T_DIHEDRAL off id off id off id off id {
  /* store the dihedral information in the last residue topology. */
  int ret = topol_add_torsion(T, NULL, $4, $3, $6, $5, $8, $7, $10, $9, $1);

  /* free the allocated strings. */
  free($4);
  free($6);
  free($8);
  free($10);

  /* check for errors. */
  if (!ret) YYERROR;
};

/* improper instance tokens:
$$        $1   $2         $3  $4 $5  $6 $7  $8 $9  $10 */
improper: mode T_IMPROPER off id off id off id off id {
  /* store the improper information in the last residue topology. */
  int ret = topol_add_improper(T, NULL, $4, $3, $6, $5, $8, $7, $10, $9, $1);

  /* free the allocated strings. */
  free($4);
  free($6);
  free($8);
  free($10);

  /* check for errors. */
  if (!ret) YYERROR;
};

/* single mass instance. */
mass: T_MASS T_WORD T_NUM {
  /* store the mass information in the topology structure. */
  int ret = topol_add_mass(T, $2, $3);

  /* free the allocated strings. */
  free($2);

  /* check for errors. */
  if (!ret) YYERROR;
};

/* atom identifier tokens. */
id: T_WORD;

/* patch mode tokens. */
mode:       { $$ = TOPOL_MODE_ADD; }
 | T_ADD    { $$ = TOPOL_MODE_ADD; }
 | T_MODIFY { $$ = TOPOL_MODE_MODIFY; }
 | T_DELETE { $$ = TOPOL_MODE_DELETE; }
 ;

/* patch offset tokens. */
off:    { $$ = 0; }
 | T_PM { $$ = $1; }
 ;

/* single type parameter. */
type: {
  /* return an empty string. */
  $$ = strdup("");
}
 | T_TYPE T_EQ T_WORD {
  /* return the type string. */
  $$ = $3;
};

/* single charge parameter. */
charge: {
  /* return an invalid charge. */
  $$ = NAN;
}
 | T_CHARGE T_EQ T_NUM {
  /* return the charge value. */
  $$ = $3;
};

/* autogeneration directive. */
autogen: T_AUTO T_ANGLE T_EQ T_BOOL T_DIHEDRAL T_EQ T_BOOL T_END {
  /* set the autogeneration options. */
  T->auto_angles = $4;
  T->auto_dihedrals = $7;
};

%%

/* topol_parse(): extract molecular topology information from a file.
 *
 * arguments:
 *  @fh: input file handle to read from.
 *  @top: pointer to the topology data structure to modify.
 *
 * returns:
 *  integer indicating whether (1) or not (0) parsing succeeded.
 *
 * warning:
 *  this function is absolutely NOT re-entrant.
 */
int topol_parse (FILE *fh, topol_t *top) {
  /* set the input file handle and output structure pointer. */
  topol_io_in = fh;
  T = top;

  /* attempt to parse the input file. */
  if (topol_io_parse()) {
    /* free the scanner buffer and return failure. */
    topol_io_clean();
    return 0;
  }

  /* free the scanner buffer. */
  topol_io_clean();

  /* return success. */
  errno = 0;
  return 1;
}

/* topol_io_error(): error reporting interface for parsing topology files.
 *
 * arguments:
 *  @msg: error message string.
 */
void topol_io_error (const char *msg) {
  /* raise an exception, to be handled by topol_parse(). */
  raise(msg);
}

