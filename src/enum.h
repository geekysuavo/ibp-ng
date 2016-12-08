
/* ensure once-only inclusion. */
#ifndef __IBPNG_ENUM_H__
#define __IBPNG_ENUM_H__

/* include the system headers. */
#include <sys/stat.h>
#include <unistd.h>
#include <fcntl.h>

/* include the pthread header. */
#ifdef __IBP_HAVE_PTHREAD
#include <pthread.h>
#endif

/* include the peptide, graph and options headers. */
#include "peptide.h"
#include "graph.h"
#include "opts.h"

/* include the vector header. */
#include "vector.h"

/* predeclare enum_t and enum_thread_t before defining them, in order
 * to allow the pruning function pointer specification below.
 */
typedef struct _enum_t enum_t;
typedef struct _enum_thread_t enum_thread_t;

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
 *  @th: pointer to the enumerator thread to check.
 *  @data: optional payload data pointer.
 *
 * returns:
 *  integer indicating whether (1) or not (0) the tree node should be
 *  pruned by the pruning function.
 */
typedef int (*enum_prune_test_fn) (struct _enum_t *E,
                                   struct _enum_thread_t *th,
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
 *  @th: pointer to the enumerator thread to access.
 *
 * returns:
 *  integer indicating whether (1) or not (0) the operation succeeded.
 */
typedef int (*enum_write_data_fn) (struct _enum_t *E,
                                   struct _enum_thread_t *th);

/* enum_write_close_fn: function pointer specification for closing
 * an enumerator data output system.
 *
 * arguments:
 *  @E: pointer to the enumerator data structure to utilize.
 */
typedef void (*enum_write_close_fn) (struct _enum_t *E);

/* enum_thread_node_t: data structure for holding the state of a single
 * node in an iDMDGP sub-tree. an array of these structures describes
 * the state of a single candidate solution.
 */
typedef struct {
  /* @idx: index [0..n-1] of the node at the current level.
   * @start: start index of the thread at the current level.
   * @end: end index of the thread at the current level.
   * @nb: number of nodes/branches at the current level.
   * @pos: position of the node for the candidate solution.
   * @prev: previous solution position of the node.
   * @energy: current energy at the node.
   */
  unsigned int idx, start, end, nb;
  vector_t pos, prev;
  double energy;
}
enum_thread_node_t;

/* enum_thread_t: data structure for holding the traversal state of a
 * single thread of an iDMDGP solution enumerator, which covers a
 * well-defined iDMDGP sub-tree.
 */
struct _enum_thread_t {
  /* @thread: system-level thread information.
   * @E: pointer back to the master enumerator.
   */
#ifdef __IBP_HAVE_PTHREAD
  pthread_t thread;
#endif
  enum_t *E;

  /* @state: array of states for each atom in the thread.
   * @level: current level of the thread in the tree.
   * @logW: logarithm of the thread size.
   */
  enum_thread_node_t *state;
  unsigned int level;
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
   * @write_mutex: mutual exclusion for multi-threaded data writes.
   * @write_open: function pointer for opening the output system.
   * @write_data: function pointer for writing output data.
   * @write_close: function pointer for closing the output system.
   */
  unsigned int *writeord;
#ifdef __IBP_HAVE_PTHREAD
  pthread_mutex_t write_mutex;
#endif
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
   */
  enum_thread_t *threads;
  unsigned int nthreads;

  /* @timer: unique thread for computing timing information.
   */
#ifdef __IBP_HAVE_PTHREAD
  pthread_t timer;
#endif

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

