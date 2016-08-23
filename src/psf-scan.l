
%{
/* include the parsing and string headers. */
#include "psf-parse.h"
#include "str.h"
%}

%option noyywrap
%option noinput
%option nounput
%%

^PSF { return T_BEGIN; }
^[ ]*(rema|REMA)[A-z]{0,3}.*$ { return T_REMARK; }

[0-9]+ {
  psf_io_lval.ival = atoi(yytext);
  return T_INT;
}

[-+]?([0-9]*\.?[0-9]+|[0-9]+\.)([eE][-+]?[0-9]+)? {
  psf_io_lval.fval = atof(yytext);
  return T_FLOAT;
}

[0-9 ]+\!NATOM { return T_ATOM; }
[0-9 ]+\![A-Z]+ { return T_OTHER; }

[A-z][A-z0-9]* {
  psf_io_lval.sval = strdup(yytext);
  return T_WORD;
}

:.*$ {}
[ \t\n] {}
. { return T_UNKNOWN; }

%%

/* psf_io_clean(): free any memory allocated for the psf scanner.
 */
void psf_io_clean (void) {
  /* call the buffer deletion function. */
  psf_io__delete_buffer(YY_CURRENT_BUFFER);
  yy_init = 1;
}

