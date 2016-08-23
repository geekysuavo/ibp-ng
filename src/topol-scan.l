
%{
/* include the molecular topology header. */
#include "topol.h"

/* include the parsing and string headers. */
#include "topol-parse.h"
#include "str.h"
%}

%option noyywrap
%option noinput
%option nounput
%x COMMENT
%%

^(checkversion|CHECKVERSION).*$ {}
^(remark|REMARK).*$ {}
(set|SET).*$ {}
[ \t\n] {}
!.*$ {}

(auto|AUTO)[A-z]{0,8} { return T_AUTO; }
(resi|RESI)[A-z]{0,3} { return T_RESIDUE; }
(pres|PRES)[A-z]{0,4} { return T_PATCH; }
(grou|GROU)(p|P)? { return T_GROUP; }
(atom|ATOM) { return T_ATOM; }
(bond|BOND) { return T_BOND; }
(angl|ANGL)[A-z]{0,2} { return T_ANGLE; }
(dihe|DIHE)[A-z]{0,5} { return T_DIHEDRAL; }
(impr|IMPR)[A-z]{0,4} { return T_IMPROPER; }
(char|CHAR)[A-z]{0,2} { return T_CHARGE; }
(type|TYPE) { return T_TYPE; }
(mass|MASS) { return T_MASS; }
(end|END) { return T_END; }

(add|ADD) { return T_ADD; }
(mod|MOD)[A-z]{0,3} { return T_MODIFY; }
(del|DEL)[A-z]{0,3} { return T_DELETE; }

(true|True|TRUE) {
  topol_io_lval.bval = 1;
  return T_BOOL;
}

(false|False|FALSE) {
  topol_io_lval.bval = 0;
  return T_BOOL;
}

[a-zA-Z][a-zA-Z0-9]* {
  topol_io_lval.sval = strdup(yytext);
  return T_WORD;
}

[-+]?([0-9]*\.?[0-9]+|[0-9]+\.) {
  topol_io_lval.fval = atof(yytext);
  return T_NUM;
}

"-" {
  topol_io_lval.ival = 0;
  return T_PM;
}

"+" {
  topol_io_lval.ival = 1;
  return T_PM;
}

"=" { return T_EQ; }
"{" { BEGIN(COMMENT); }

. { return T_UNKNOWN; }

<COMMENT>"}" { BEGIN(INITIAL); }
<COMMENT>\n {}
<COMMENT>[^}]+ {}

%%

/* topol_io_clean(): free any memory allocated for the topology scanner.
 */
void topol_io_clean (void) {
  /* call the buffer deletion function. */
  topol_io__delete_buffer(YY_CURRENT_BUFFER);
  yy_init = 1;
}

