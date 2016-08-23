
%{
/* include the assignment header. */
#include "assign.h"

/* include the parsing header. */
#include "assign-parse.h"
%}

%option noyywrap
%option noinput
%option nounput
%x COMMENT
%%

(assi|ASSI)(gn|GN)? { return T_ASSIGN; }

(resi|RESI)(d|due|D|DUE)? { return T_RESIDUE; }
(resn|RESN)(ame|AME)? { return T_RESNAME; }
(name|NAME) { return T_NAME; }
(id|ID) { return T_ID; }

(or|OR)     { return T_OR; }
(and|AND)   { return T_AND; }
(not|NOT)   { return T_NOT; }
(all|ALL)   { return T_ALL; }
(none|NONE) { return T_NONE; }

(hydr|HYDR)(ogen|OGEN)? { assign_io_lval.cval = 'H'; return T_TYPE; }
(carb|CARB)(on|ON)?     { assign_io_lval.cval = 'C'; return T_TYPE; }
(nitr|NITR)(ogen|OGEN)? { assign_io_lval.cval = 'N'; return T_TYPE; }
(oxyg|oxgy)(gen|GEN)?   { assign_io_lval.cval = 'O'; return T_TYPE; }

[-+]?[0-9]+ {
  assign_io_lval.ival = atoi(yytext);
  return T_INT;
}

[-+]?([0-9]*\.?[0-9]+|[0-9]+\.)([eE][-+]?[0-9]+)? {
  assign_io_lval.fval = atof(yytext);
  return T_FLOAT;
}

[A-z][A-z0-9#?]* {
  assign_io_lval.sval = strdup(yytext);
  return T_WORD;
}

"(" { return T_PAREN_OPEN; }
")" { return T_PAREN_CLOSE; }
"{" { BEGIN(COMMENT); }

[ \t] {}
\n { return T_EOL; }
. { return T_UNKNOWN; }

<COMMENT>"}" { BEGIN(INITIAL); }
<COMMENT>\n {}
<COMMENT>[^}]+ {}

%%

/* assign_io_clean(): free any memory allocated for the assignment scanner.
 */
void assign_io_clean (void) {
  /* call the buffer deletion function. */
  assign_io__delete_buffer(YY_CURRENT_BUFFER);
  yy_init = 1;
}

