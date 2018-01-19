
/* ensure once-only inclusion. */
#ifndef __IBPNG_OPTS_H__
#define __IBPNG_OPTS_H__

/* include the traceback and string headers. */
#include "trace.h"
#include "str.h"

/* opts_t: option task structure for storing all parsed command line
 * options for a given application instance.
 */
typedef struct {
  /* declare status variables for verbosity and help output:
   *  @help: whether or not to display the help message string.
   *  @verb: whether or not to output verbose messages.
   */
  unsigned int help, verb;

  /* declare variables for input filename storage:
   *  @fname_in: input filename string.
   *  @fname_top: topology filename string.
   *  @fname_par: parameter filename string.
   *  @fname_ord: reorder filename string.
   */
  char *fname_in;
  char *fname_top;
  char *fname_par;
  char *fname_ord;

  /* declare variables for output filename storage:
   *  @fname_out: output filename string.
   *  @fname_psf: psf filename string.
   *  @fname_dmdgp: dmdgp filename string.
   */
  char *fname_out;
  char *fname_psf;
  char *fname_dmdgp;

  /* declare variables for input and output clarifications:
   *  @idx_in: chain or index string for input file parsing.
   *  @fmt_out: format string for output file writing.
   */
  char *idx_in;
  char *fmt_out;

  /* declare variables for restraint filename storage:
   *  @fname_restr: array of restraint filename strings.
   *  @n_restr: number of restraint filename strings.
   */
  char **fname_restr;
  unsigned int n_restr;

  /* declare variables for sidechain index storage:
   *  @sidech: array of sidechain indices.
   *  @n_sidech: number of sidechain indices.
   */
  unsigned int *sidech, n_sidech;

  /* declare variables for pruning device name storage:
   *  @prune: array of pruning device names.
   *  @n_prune: number of pruning device names.
   */
  char **prune;
  unsigned int n_prune;

  /* declare variables for branch control:
   *  @thread_gpu: whether or not we should use the gpu.
   *  @thread_num: number of parallel threads to utilize.
   *  @branch_max: maximum number of branches per node.
   *  @branch_eps: smallest division for interval discretization.
   */
  unsigned int thread_gpu, thread_num;
  unsigned int branch_max;
  double branch_eps;

  /* declare variables for pruning control:
   *  @nsol_limit: maximum number of solutions to enumerate.
   *  @complete: whether or not to complete the graph edge set.
   *  @vdw_scale: atomic radius scaling factor for ddf lower-bounds.
   *  @ddf_tol: tolerance for acceptable out-of-bound errors.
   *  @rmsd_tol: rmsd for skipping structures.
   */
  unsigned int nsol_limit, complete;
  double vdw_scale;
  double ddf_tol;
  double rmsd_tol;
}
opts_t;

/* function declarations: */

opts_t *opts_new (void);

opts_t *opts_new_from_strings (int argc, char **argv);

void opts_free (opts_t *opts);

int opts_validate (opts_t *opts);

#endif  /* !__IBPNG_OPTS_H__ */

