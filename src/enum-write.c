
/* include the enumerator headers. */
#include "enum.h"
#include "enum-node.h"

/* all functions defined in this source file must follow the function
 * pointer specifications outlined for enumerator output writing systems.
 *
 * for more details, consult:
 *  - enum_write_open_function   \
 *  - enum_write_data_function    }-- all in enum.h
 *  - enum_write_close_function  /
 */

/* * * * * * * * * * * * * * DCD: * * * * * * * * * * * * * */

/* enum_write_dcd_open(): called to open a DCD output system.
 */
int enum_write_dcd_open (enum_t *E) {
  /* declare required variables:
   *  @fd: temporary file descriptor.
   *  @sz: buffer size for data records.
   *  @ibuf, @dbuf, @cbuf: data buffers.
   */
  int fd, sz;
  int ibuf[9];
  double dbuf;
  char cbuf[160];

  /* attempt to open the output file. */
  fd = open(E->fname, O_CREAT | O_TRUNC | O_WRONLY,
            S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP | S_IROTH);

  /* check if the file failed to open. */
  if (fd < 0)
    throw("unable to open file '%s'", E->fname);

  /* success! store the file descriptor. */
  E->fd = fd;

  /* compute the size of the first data record. */
  sz = 4 * sizeof(char) + 18 * sizeof(int) + sizeof(double);

  /* initialize the buffers. */
  dbuf = 1.0;
  memset(ibuf, 0, 9 * sizeof(int));
  memset(cbuf, 0, 32 * sizeof(char));
  strcpy(cbuf, "CORD");

  /* write the first data record. */
  write(E->fd, &sz, sizeof(int));
  write(E->fd, cbuf, 4 * sizeof(char));
  write(E->fd, ibuf, 9 * sizeof(int));
  write(E->fd, &dbuf, sizeof(double));
  write(E->fd, ibuf, 9 * sizeof(int));
  write(E->fd, &sz, sizeof(int));

  /* compute the size of the second data record. */
  sz = 160 * sizeof(char) + sizeof(int);

  /* initialize the buffers. */
  memset(cbuf, 0, 160 * sizeof(char));
  strcpy(cbuf, "Created by ibp-ng");
  ibuf[0] = 2;

  /* write the second data record. */
  write(E->fd, &sz, sizeof(int));
  write(E->fd, ibuf, sizeof(int));
  write(E->fd, cbuf, 160 * sizeof(char));
  write(E->fd, &sz, sizeof(int));

  /* compute the size of the third data record and
   * initialize the buffer.
   */
  sz = sizeof(int);
  ibuf[0] = E->G->n_orig;

  /* write the third data record. */
  write(E->fd, &sz, sizeof(int));
  write(E->fd, ibuf, sizeof(int));
  write(E->fd, &sz, sizeof(int));

  /* return success. */
  return 1;
}

/* enum_write_dcd_close(): called to close a DCD output system.
 */
void enum_write_dcd_close (enum_t *E) {
  /* close the file descriptor. */
  if (E->fd >= 0)
    close(E->fd);

  /* re-init the file descriptor. */
  E->fd = -1;
}

/* enum_write_dcd(): called to write a structure to (32-bit) DCD output.
 */
int enum_write_dcd (enum_t *E) {
  /* declare required variables:
   *  @xyz: current coordinate value, casted to a float.
   *  @sz: size of each coordinate (x, y, z) record.
   */
  float xyz;
  int sz;

  /* compute and write the first size integer. */
  sz = E->G->n_orig * sizeof(float);
  write(E->fd, &sz, sizeof(int));

  /* locally store the thread length and originality array. */
  const unsigned int n = E->G->nv;
  const unsigned int max = E->G->n_order;
  const unsigned int *rev = E->G->ordrev;

  /* x: write the current thread coordinates. */
  for (unsigned int i = 0; i < n; i++) {
    /* skip duplicate atoms. */
    if (rev[i] >= max) continue;

    /* write the value. */
    xyz = (float) E->soln[rev[i]].x;
    write(E->fd, &xyz, sizeof(float));
  }

  /* write a pair of size integers. */
  write(E->fd, &sz, sizeof(int));
  write(E->fd, &sz, sizeof(int));

  /* y: write the current thread coordinates. */
  for (unsigned int i = 0; i < n; i++) {
    /* skip duplicate atoms. */
    if (rev[i] >= max) continue;

    /* write the value. */
    xyz = (float) E->soln[rev[i]].y;
    write(E->fd, &xyz, sizeof(float));
  }

  /* write another pair of size integers. */
  write(E->fd, &sz, sizeof(int));
  write(E->fd, &sz, sizeof(int));

  /* z: write the current thread coordinates. */
  for (unsigned int i = 0; i < n; i++) {
    /* skip duplicate atoms. */
    if (rev[i] >= max) continue;

    /* write the value. */
    xyz = (float) E->soln[rev[i]].z;
    write(E->fd, &xyz, sizeof(float));
  }

  /* write one final size integer. */
  write(E->fd, &sz, sizeof(int));

  /* return success. */
  return 1;
}

/* * * * * * * * * * * * * * PDB: * * * * * * * * * * * * * */

/* enum_write_pdb_open(): called to open a PDB output system.
 */
int enum_write_pdb_open (enum_t *E) {
  /* attempt to create the output directory. */
  if (mkdir(E->fname, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH))
    throw("unable to create directory '%s'", E->fname);

  /* return success. */
  return 1;
}

/* enum_write_pdb(): called to write a structure to PDB output.
 */
int enum_write_pdb (enum_t *E) {
  /* declare required variables:
   * @n: number of unique tree node pointers in the path.
   */
  unsigned int isol, i, l, n;
  peptide_atom_t *atom;
  char *fname;
  FILE *fh;

  /* allocate the filename string. */
  fname = (char*) malloc((strlen(E->fname) + 64) * sizeof(char));
  if (!fname)
    throw("unable to allocate filename string");

  /* get the solution index. */
  isol = E->nsol;

  /* construct the filename string. */
  sprintf(fname, "%s/%08u.pdb", E->fname, isol);

  /* open the output file. */
  fh = fopen(fname, "w");
  if (!fh) {
    /* free allocated memory and return failure. */
    raise("unable to open '%s' for writing", fname);
    free(fname);
    return 0;
  }

  /* print header information. */
  fprintf(fh, "HEADER    ibp-ng\n");
  fprintf(fh, "TITLE     ibp-ng solution %-24u\n", isol);
  fprintf(fh, "%-78s\n", "COMPND");

  /* output the residue sequence information. */
  fprintf(fh, "SEQRES %-3u %c %4u  ", (l = 1), 'A', E->P->n_res);
  for (i = 0; i < E->P->n_res; i++) {
    fprintf(fh, "%s ", peptide_get_resname(E->P, i));

    if (((i + 1) % 13) == 0 && i < E->P->n_res - 1)
      fprintf(fh, "\nSEQRES %-3u %c %4u  ", l++, 'A', E->P->n_res);
    else if (i == E->P->n_res - 1)
      fprintf(fh, "\n");
  }

  /* loop over the current thread state. */
  fprintf(fh, "%-6s    %-4u\n", "MODEL", 1);
  for (i = n = 0; i < E->G->nv; i++) {
    /* skip duplicate atoms. */
    if (E->G->ordrev[i] >= E->G->n_order)
      continue;

    /* write the current atom information. */
    atom = E->P->atoms + i;
    fprintf(fh, "%-6s%5u %-4s %3s %c%4u    %8.3lf%8.3lf%8.3lf"
                "%6.2lf%6.2lf           %c\n",
            "ATOM", n++, atom->name,
            peptide_get_resname(E->P, atom->res_id),
            'A', atom->res_id + 1,
            E->soln[E->G->ordrev[i]].x,
            E->soln[E->G->ordrev[i]].y,
            E->soln[E->G->ordrev[i]].z,
            1.0, 0.0,
            atom->type[0]);
  }

  /* print footer information. */
  fprintf(fh, "ENDMDL\nEND\n");

  /* clean up and return success. */
  fclose(fh);
  free(fname);
  return 1;
}

