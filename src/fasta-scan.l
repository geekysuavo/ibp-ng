
%{
/* include the parsing and string headers. */
#include "fasta-parse.h"
#include "str.h"
%}

%option noyywrap
%option noinput
%option nounput
%x HEADER
%%

^;.*$ {}

">" {
  BEGIN(HEADER);
  return T_BEGIN;
}

[A-z] {
  fasta_io_lval.cval = yytext[0];
  return T_RESIDUE;
}

\n { return T_EOL; }
. { return T_UNKNOWN; }

<HEADER>\n { BEGIN(INITIAL); }
<HEADER>[^\n]+ {}

%%

/* fasta_io_clean(): free any memory allocated for the fasta scanner.
 */
void fasta_io_clean (void) {
  /* call the buffer deletion function. */
  fasta_io__delete_buffer(YY_CURRENT_BUFFER);
  yy_init = 1;
}

