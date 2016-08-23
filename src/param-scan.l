
%{
/* include the molecular parameter header. */
#include "param.h"

/* include the parsing and string headers. */
#include "param-parse.h"
#include "str.h"
%}

%option noyywrap
%option noinput
%option nounput
%x COMMENT NBOND
%%

^(checkversion|CHECKVERSION).*$ {}
^(remarks|REMARKS).*$ {}
(set|SET).*$ {}
[ \t\n] {}
!.*$ {}

(bond|BOND) { return T_BOND; }
(angl|ANGL)[eE]? { return T_ANGLE; }
(dihe|DIHE)[A-z]{0,4} { return T_TORSION; }
(impr|IMPR)[A-z]{0,4} { return T_IMPROPER; }
(nonb|NONB)[A-z]{0,5} { return T_NONBOND; }
(nbon|NBON)[A-z]{0,2} { BEGIN(NBOND); }
(end|END) { return T_END; }

[a-zA-Z][a-zA-Z0-9]* {
  param_io_lval.sval = strdup(yytext);
  return T_WORD;
}

[-+]?([0-9]*\.?[0-9]+|[0-9]+\.) {
  param_io_lval.fval = atof(yytext);
  return T_NUM;
}

"[" { return T_OPEN; }
"]" { return T_CLOSE; }
"," { return T_COMMA; }

"{" { BEGIN(COMMENT); }

. { return T_UNKNOWN; }

<COMMENT>"}" { BEGIN(INITIAL); }
<COMMENT>\n {}
<COMMENT>[^}]+ {}

<NBOND>^[ ]*(end|END) { BEGIN(INITIAL); }
<NBOND>\n {}
<NBOND>. {}

%%

/* param_io_clean(): free any memory allocated for the parameter scanner.
 */
void param_io_clean (void) {
  /* call the buffer deletion function. */
  param_io__delete_buffer(YY_CURRENT_BUFFER);
  yy_init = 1;
}

