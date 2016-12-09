
/* ensure once-only inclusion. */
#ifndef __IBPNG_ENUM_H__
#define __IBPNG_ENUM_H__

/* include the system headers. */
#include <sys/stat.h>
#include <unistd.h>
#include <fcntl.h>

/* include the peptide, graph and options headers. */
#include "peptide.h"
#include "graph.h"
#include "opts.h"

/* include the vector header. */
#include "vector.h"

/* predeclare enum_t and enum_node_t before defining them, in order
 * to allow the pruning function pointer specification below.
 */
typedef struct _enum_t enum_t;
typedef struct _enum_node_t enum_node_t;

/* enum_prune_init_fn: function pointer specification for 
 * initializing an enumerator pruning device.
 *
 * arguments:
 *  @E: pointer to the enumerator data structure to modify.
 *  @lev: level in the graph repetition order to initialize prunes for.
 *
 * returns:
 *  integer indicating whether (1) or not (0) the operation succeeded.
 */
typedef int (*enum_prune_init_fn) (struct _enum_t *E, unsigned int lev);

/* enum_prune_test_fn: function pointer specification for testing
 * prune-ability (i.e. infeasibility) by an enumerator pruning device.
 *
 * arguments:
 *  @E: pointer to the enumerator data structure to utilize.
 *  @end: pointer to the tree node to check.
 *  @data: optional payload data pointer.
 *
 * returns:
 *  integer indicating whether (1) or not (0) the tree node should be
 *  pruned by the pruning function.
 */
typedef int (*enum_prune_test_fn) (struct _enum_t *E,
                                   struct _enum_node_t *end,
                                   void *data);

/* enum_prune_report_fn: function pointer specification for
 * printing reporting information on a single pruning closure.
 *
 * arguments:
 *  @E: pointer to the enumerator data structure to access.
 *  @lev: level in the graph repetition order to report prunes for.
 *  @data: optional payload data pointer.
 */
typedef void (*enum_prune_report_fn) (struct _enum_t *E,
                                      unsigned int lev,
                                      void *data);

/* enum_write_open_fn: function pointer specification for initializing
 * an enumerator data output system.
 *
 * arguments:
 *  @E: pointer to the enumerator data structure to utilize.
 *
 * returns:
 *  integer indicating whether (1) or not (0) the operation succeeded.
 */
typedef int (*enum_write_open_fn) (struct _enum_t *E);

/* enum_write_data_fn: function pointer specification for writing
 * solution data through an enumerator data output system.
 *
 * arguments:
 *  @E: pointer to the enumerator data structure to utilize.
 *
 * returns:
 *  integer indicating whether (1) or not (0) the operation succeeded.
 */
typedef int (*enum_write_data_fn) (struct _enum_t *E);

/* enum_write_close_fn: function pointer specification for closing
 * an enumerator data output system.
 *
 * arguments:
 *  @E: pointer to the enumerator data structure to utilize.
 */
typedef void (*enum_write_close_fn) (struct _enum_t *E);

/* enum_node_t: tree data structure for holding the traversal state of
 * an iDMDGP solution enumerator.
 */
struct _enum_node_t {
  /* @pos: embedded coordinates of the current atom.
   * @energy: current energy at the node.
   */
  vector_t pos;
  double energy;

  /* @feas: whether the node's branch is currently considered feasible.
   */
  unsigned int feas;

  /* @lev: current level of the tree node.
   * @idx: index [0..nb-1] of the tree node.
   * @nb: number of branches from the tree node.
   */
  unsigned int lev, idx, nb;

  /* @prev: parent tree node pointer, or null for the root.
   * @next: array of child tree node pointers.
   */
  enum_node_t *prev;
  enum_node_t **next;
};

/* enum_t: structure for holding all state information required for the
 * complete enumeration (hence 'enum' for 'enumerator') of feasible
 * solutions from an iDMDGP instance graph.
 *
 * an enumerator stores pointers to the peptide and graph data structures
 * that solutions are being computed from.
 *
 * the enumerator state is stored in an explicit tree data structure, and
 * tree traversal is accomplished using a recursive post-order depth-first
 * search.
 */
struct _enum_t {
  /* @P: pointer to the related peptide data structure.
   * @G: pointer to the related graph data structure.
   */
  peptide_t *P;
  graph_t *G;

  /* @writeord: atom ordering to use when writing data.
   * @write_open: function pointer for opening the output system.
   * @write_data: function pointer for writing output data.
   * @write_close: function pointer for closing the output system.
   */
  unsigned int *writeord;
  enum_write_open_fn write_open;
  enum_write_data_fn write_data;
  enum_write_close_fn write_close;

  /* @term: flag to terminate the enumeration.
   */
  unsigned int term;

  /* @logW: logarithm of the number of leaves in the tree.
   * @nsol: number of solutions computed from the graph.
   * @nmax: maximum number of solutions to compute.
   * @fname: file/directory name string for storing outputs.
   * @fd: file descriptor for DCD-formatted output.
   */
  unsigned int nsol, nmax;
  double logW;
  char *fname;
  int fd;

  /* @nbmax: maximum number of branches at each tree level.
   * @eps: minimum discretization size for distance intervals.
   */
  unsigned int nbmax;
  double eps;

  /* @prune: (2d) array of pruning test function pointers.
   * @prune_sz: sizes of each inner array in @prune_test.
   * @prune_data: array of pruning data payloads.
   */
  enum_prune_test_fn **prune;
  unsigned int *prune_sz;
  void ***prune_data;

  /* @threads: array of enumerator threads.
   * @nthreads: number of enumerator threads.
   * @soln: vertices of the last solution accepted.
   */
  enum_node_t tree;
  vector_t *soln;

  /* pruning function control variables:
   *  @ddf_tol: error tolerance for ddf bounds checking.
   *  @rmsd_tol: minimum acceptable rmsd between solutions.
   *  @energy_tol: maximum acceptable energy for pruning.
   */
  double ddf_tol, rmsd_tol, energy_tol;
};

/* function declarations (enum.c): */

enum_t *enum_new (peptide_t *P, graph_t *G, opts_t *opts);

void enum_free (enum_t *E);

int enum_execute (enum_t *E);

#endif  /* !__IBPNG_ENUM_H__ */

