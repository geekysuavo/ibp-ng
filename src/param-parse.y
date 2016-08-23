
/* enable verbose error reporting. */
%error-verbose
%debug

%{
/* include the molecular parameter header. */
#include "param.h"

/* include the parsing header. */
#include "param-parse.h"

/* provide a link to the input file handle and output structure pointer.
 */
extern FILE *param_io_in;
param_t *P;

/* flex function declarations: */
void param_io_error (const char *msg);
void param_io_clean (void);
int param_io_lex (void);
%}

/* define the yylval type union. */
%union {
  value_t gval;
  float fval;
  char *sval;
}

/* define all recognized tokens. */
%token T_BOND T_ANGLE T_TORSION T_IMPROPER T_NONBOND
%token T_WORD T_NUM T_END T_OPEN T_CLOSE T_COMMA
%token T_UNKNOWN

/* define the types of parsed values. */
%type<gval> scalar interval value
%type<sval> T_WORD
%type<fval> T_NUM

%%

/* global listing of all parameter information. */
parameters: parameter
 | parameters parameter
 ;

/* single parameter instance. */
parameter: bond
 | angle
 | torsion
 | improper
 | nonbond
 ;

/* bond instance tokens. */
bond: T_BOND T_WORD T_WORD T_NUM value {
  /* store the information in the parameter structure. */
  if (!param_add_bond(P, $2, $3, $5))
    YYERROR;
};

/* angle instance tokens. */
angle: T_ANGLE T_WORD T_WORD T_WORD T_NUM value {
  /* store the information in the parameter structure. */
  if (!param_add_angle(P, $2, $3, $4, $6))
    YYERROR;
};

/* torsion instance tokens. */
torsion: T_TORSION T_WORD T_WORD T_WORD T_WORD T_NUM T_NUM value {
  /* store the information in the parameter structure. */
  if (!param_add_torsion(P, $2, $3, $4, $5, $8))
    YYERROR;
};

/* improper instance tokens. */
improper: T_IMPROPER T_WORD T_WORD T_WORD T_WORD T_NUM T_NUM value {
  /* store the information in the parameter structure. */
  if (!param_add_improper(P, $2, $3, $4, $5, $8))
    YYERROR;
};

/* nonbonded instance tokens. */
nonbond: T_NONBOND T_WORD T_NUM T_NUM T_NUM T_NUM {
  /* store the information in the parameter structure. */
  if (!param_add_radius(P, $2, $4))
    YYERROR;
};

/* generalized parameter values. */
value: interval | scalar ;

/* interval parameter values. */
interval: T_OPEN T_NUM T_COMMA T_NUM T_CLOSE {
  /* store the interval parameter value. */
  $$ = value_interval($2, $4);
};

/* scalar parameter values. */
scalar: T_NUM {
  /* store the scalar parameter value. */
  $$ = value_scalar($1);
};

%%

/* param_parse(): extract molecular parameter information from a file.
 *
 * arguments:
 *  @fh: input file handle to read from.
 *  @par: pointer to the parameter data structure to modify.
 *
 * returns:
 *  integer indicating whether (1) or not (0) parsing succeeded.
 *
 * warning:
 *  this function is absolutely NOT re-entrant.
 */
int param_parse (FILE *fh, param_t *par) {
  /* set the input file handle and output structure pointer. */
  param_io_in = fh;
  P = par;

  /* attempt to parse the input file. */
  if (param_io_parse()) {
    /* free the scanner buffer and return failure. */
    param_io_clean();
    return 0;
  }

  /* free the scanner buffer. */
  param_io_clean();

  /* return success. */
  errno = 0;
  return 1;
}

/* param_io_error(): error reporting interface for parsing parameter files.
 *
 * arguments:
 *  @msg: error message string.
 */
void param_io_error (const char *msg) {
  /* raise an exception, to be handled by param_parse(). */
  raise(msg);
}

