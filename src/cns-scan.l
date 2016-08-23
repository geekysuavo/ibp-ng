
%{
/* include the parsing and string headers. */
#include "cns-parse.h"
#include "str.h"
%}

%option noyywrap
%option noinput
%option nounput
%x COMMENT
%%

^data_cns_mtf { return T_BEGIN; }
^loop_ { return T_LOOP; }

[-+]?[0-9]+ {
  cns_io_lval.ival = atoi(yytext);
  return T_INT;
}

[-+]?([0-9]*\.?[0-9]+|[0-9]+\.)([eE][-+]?[0-9]+)? {
  cns_io_lval.fval = atof(yytext);
  return T_FLOAT;
}

^_cns_mtf_atom\..*$ { return T_ATOM; }
^_.*$ { return T_OTHER; }

'[A-z0-9 ][A-z0-9]*' {
  cns_io_lval.sval = strdup(yytext + 1);
  return T_WORD;
}

";" { BEGIN(COMMENT); }

:.*$ {}
[ \t\n] {}
. { return T_UNKNOWN; }

<COMMENT>";" { BEGIN(INITIAL); }
<COMMENT>\n {}
<COMMENT>[^;]+ {}

%%

/* cns_io_clean(): free any memory allocated for the cns scanner.
 */
void cns_io_clean (void) {
  /* call the buffer deletion function. */
  cns_io__delete_buffer(YY_CURRENT_BUFFER);
  yy_init = 1;
}

