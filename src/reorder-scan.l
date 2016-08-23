
%{
/* include the parsing and string headers. */
#include "reorder-parse.h"
#include "str.h"
%}

%option noyywrap
%option noinput
%option nounput
%x COMMENT
%%

[ \t\n,] {}
!.*$ {}

(re|RE)?(ord|ORD)(e|er|E|ER)? { return T_REORDER; }
(end|END) { return T_END; }

(opt|OPT)(ional|IONAL)? { return T_OPT; }
(prev|PREV)(ious|IOUS)? { return T_PREV; }
(cur|CUR)(rent|RENT)? { return T_CURR; }
(next|NEXT) { return T_NEXT; }

[a-zA-Z][a-zA-Z0-9]* {
  reorder_io_lval.sval = strdup(yytext);
  return T_WORD;
}

"*" { return T_OPT; }
"-" { return T_PREV; }
"+" { return T_NEXT; }
"{" { BEGIN(COMMENT); }

. { return T_UNKNOWN; }

<COMMENT>"}" { BEGIN(INITIAL); }
<COMMENT>\n {}
<COMMENT>[^}]+ {}

%%

/* reorder_io_clean(): free any memory allocated for the re-order scanner.
 */
void reorder_io_clean (void) {
  /* call the buffer deletion function. */
  reorder_io__delete_buffer(YY_CURRENT_BUFFER);
  yy_init = 1;
}

