
/* ensure once-only inclusion. */
#ifndef __IBPNG_PEPTIDE_ALLOC_H__
#define __IBPNG_PEPTIDE_ALLOC_H__

/* function declarations (peptide-alloc.c): */

int fasta_parse (FILE *fh, unsigned int sidx, peptide_t *pep);

int pdb_parse (FILE *fh, char chain, peptide_t *pep);

int psf_parse (FILE *fh, peptide_t *pep);

int cns_parse (FILE *fh, peptide_t *pep);

#endif /* !__IBPNG_PEPTIDE_ALLOC_H__ */

