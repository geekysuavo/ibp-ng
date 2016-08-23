
%{
/* include the parsing and string headers. */
#include "pdb-parse.h"
#include "str.h"
%}

%option noyywrap
%option noinput
%option nounput
%%

^SEQRES { return T_SEQRES; }

[0-9]+ {
  pdb_io_lval.ival = atoi(yytext);
  return T_INT;
}

[A-z][A-z0-9]* {
  pdb_io_lval.sval = strdup(yytext);
  return T_WORD;
}

\n { return T_EOL; }
[ \t] {}
. {}

%%

/* pdb_io_clean(): free any memory allocated for the pdb scanner.
 */
void pdb_io_clean (void) {
  /* call the buffer deletion function. */
  pdb_io__delete_buffer(YY_CURRENT_BUFFER);
  yy_init = 1;
}

