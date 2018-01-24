
/* ensure once-only inclusion. */
#pragma once

/* function declarations (peptide-alloc.c): */

int fasta_parse (FILE *fh, unsigned int sidx, peptide_t *pep);

int pdb_parse (FILE *fh, char chain, peptide_t *pep);

int psf_parse (FILE *fh, peptide_t *pep);

int cns_parse (FILE *fh, peptide_t *pep);

