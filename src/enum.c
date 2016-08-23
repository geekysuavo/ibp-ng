
/* include the enumerator headers. */
#include "enum.h"
#include "enum-thread.h"
#include "enum-write.h"
#include "enum-prune.h"

/* enum_format_map_t: structure for mapping between output format names
 * and enumerated type values.
 */
struct enum_format_map_t {
  /* @name: string name of the format.
   * @write_open, @write_data, @write_close: format function pointers.
   */
  char *name;
  enum_write_open_fn write_open;
  enum_write_data_fn write_data;
  enum_write_close_fn write_close;
};

/* enum_prune_map_t: structure for mapping between pruning method names
 * and function pointers.
 */
struct enum_prune_map_t {
  /* @name: string name of the pruning method.
   * @prune_init, @prune_test, @prune_report: pruning function pointers.
   */
  char *name;
  enum_prune_init_fn prune_init;
  enum_prune_test_fn prune_test;
  enum_prune_report_fn prune_report;
};

/* formats: mapping between name and type of all output formats
 * available to the ibp-ng enumerator.
 */
static const struct enum_format_map_t formats[] = {
  /* null output format. writes absolutely nothing. */
  { "null", NULL, NULL, NULL },

  /* dcd output format. creates a single dcd trajectory file. */
  { "dcd",
    enum_write_dcd_open,
    enum_write_dcd,
    enum_write_dcd_close
  },

  /* pdb output format. creates a directory full of pdb files. */
  { "pdb",
    enum_write_pdb_open,
    enum_write_pdb,
    NULL /* no close function required. */
  },

  /* null-terminator. */
  { NULL, NULL, NULL, NULL }
};

/* pruners: mapping between name and pointer of all pruning method
 * available to the ibp-ng enumerator.
 */
static const struct enum_prune_map_t pruners[] = {
  /* direct distance feasibility. */
  { "dist",
    enum_prune_ddf_init,
    enum_prune_ddf,
    enum_prune_ddf_report
  },

  /* dihedral torsion feasibility. */
  { "dihe",
    enum_prune_dihe_init,
    enum_prune_taf,
    enum_prune_dihe_report
  },

  /* improper torsion feasibility. */
  { "impr",
    enum_prune_impr_init,
    enum_prune_taf,
    enum_prune_impr_report
  },

  /* shortest path feasibility. */
  { "path",
    enum_prune_path_init,
    enum_prune_path,
    enum_prune_path_report
  },

  /* null-terminator. */
  { NULL, NULL, NULL, NULL }
};

/* enum_init_threads(): set up the thread array of an enumerator in
 * preparation for execution.
 *
 * arguments:
 *  @E: pointer to the enumerator structure to modify.
 *  @opts: pointer to an options data structure to access.
 *
 * returns:
 *  integer indicating whether (1) or not (0) the operation succeeded.
 */
static int enum_init_threads (enum_t *E, opts_t *opts) {
  /* declare required variables:
   *  @bytes: number of bytes to allocate for the thread array.
   *  @offset: offset for initializing thread states.
   *  @stride: stride for initializing thread states.
   */
  unsigned long bytes, offset, stride;
  char *stateptr;

  /* initialize the thread array and count. */
  E->threads = NULL;
  E->nthreads = opts->thread_num;

  /* compute the number of bytes to allocate. */
  bytes = E->G->n_order * sizeof(enum_thread_node_t);
  bytes = E->nthreads * (sizeof(enum_thread_t) + bytes);

  /* compute the thread offset and stride. */
  offset = E->nthreads * sizeof(enum_thread_t);
  stride = E->G->n_order * sizeof(enum_thread_node_t);

  /* allocate the array of threads. */
  E->threads = (enum_thread_t*) malloc(bytes);
  if (!E->threads)
    throw("unable to allocate array of %u threads (%lu bytes)",
          E->nthreads, bytes);

  /* initialize the thread contents. */
  for (unsigned int i = 0; i < E->nthreads; i++) {
    /* store the enumerator pointer and set the initial level. */
    E->threads[i].E = E;
    E->threads[i].level = 3;

    /* initialize the thread state pointer. */
    stateptr = ((char*) E->threads) + offset + i * stride;
    E->threads[i].state = (enum_thread_node_t*) stateptr;

    /* loop over the positions in the order. */
    for (unsigned int j = 0; j < E->G->n_order; j++) {
      /* set the node states. */
      E->threads[i].state[j].idx = 0;
      E->threads[i].state[j].end = 0;
      E->threads[i].state[j].nb = 0;

      /* set the node coordinates. */
      vector_set(&E->threads[i].state[j].pos, 0.0, 0.0, 0.0);
      vector_set(&E->threads[i].state[j].prev, 1.0e6, 1.0e6, 1.0e6);
    }
  }

  /* return success. */
  return 1;
}

/* enum_init_format(): set the output format of an enumerator by the string
 * name of the format.
 *
 * arguments:
 *  @E: pointer to the enumerator structure to modify.
 *  @opts: pointer to an options data structure to access.
 *
 * returns:
 *  integer indicating whether (1) or not (0) the operation succeeded.
 */
static int enum_init_format (enum_t *E, opts_t *opts) {
  /* declare required variables:
   *  @i: general-purpose loop counter.
   */
  unsigned int i;

  /* initialize with the default output system. */
  E->write_data = enum_write_dcd;
  E->write_open = enum_write_dcd_open;
  E->write_close = enum_write_dcd_close;

#ifdef __IBP_HAVE_PTHREAD
  /* initialize the write mutex. */
  pthread_mutex_init(&E->write_mutex, NULL);
#endif

  /* return if no output format was specified. */
  if (!opts->fmt_out)
    return 1;

  /* search for the format in the mapping. */
  for (i = 0; formats[i].name; i++) {
    /* check if the current format name matches. */
    if (strcmp(formats[i].name, opts->fmt_out) == 0) {
      /* match found. store the function pointers. */
      E->write_data = formats[i].write_data;
      E->write_open = formats[i].write_open;
      E->write_close = formats[i].write_close;

      /* return success. */
      return 1;
    }
  }

  /* unknown format name. */
  throw("unrecognized output format '%s'", opts->fmt_out);
}

/* enum_init_prune_add(): register a pruning device with an enumerator by
 * the string name of the pruning device.
 *
 * arguments:
 *  @E: pointer to the enumerator structure to modify.
 *  @name: pruning method name string.
 *
 * returns:
 *  integer indicating whether (1) or not (0) the operation succeeded.
 */
static int enum_init_prune_add (enum_t *E, const char *name) {
  /* declare required variables:
   *  @initfn: initialization function pointer.
   *  @i: general-purpose loop counter.
   */
  enum_prune_init_fn initfn;
  unsigned int i, lev;

  /* search for the pruning method in the mapping. */
  for (i = 0; pruners[i].name; i++) {
    /* check if the current pruning method name matches. */
    if (strcmp(pruners[i].name, name) == 0) {
      /* match found. get the pruning initialization function. */
      initfn = pruners[i].prune_init;

      /* loop over all levels of the graph order. */
      for (lev = 0; lev < E->G->n_order; lev++) {
        /* skip duplicate levels. */
        if (E->G->orig[lev])
          continue;

        /* initialize pruning methods for the current level. */
        if (!initfn(E, lev))
          throw("unable to initialize pruning method '%s'", name);
      }

      /* return success. */
      return 1;
    }
  }

  /* unknown pruning method. */
  throw("unrecognized pruning method '%s'", name);
}

/* initialize the pruning methods of an enumerator.
 *
 * arguments:
 *  @E: pointer to the enumerator structure to modify.
 *  @opts: pointer to an options data structure to access.
 *
 * returns:
 *  integer indicating whether (1) or not (0) the operation succeeded.
 */
static int enum_init_prune (enum_t *E, opts_t *opts) {
  /* declare required variables:
   *  @i: general-purpose loop counter.
   */
  unsigned int i;

  /* allocate the pruning method outer array. */
  E->prune = (enum_prune_test_fn**)
    malloc(E->G->n_order * sizeof(enum_prune_test_fn*));

  /* allocate the pruning method count array. */
  E->prune_sz = (unsigned int*)
    malloc(E->G->n_order * sizeof(unsigned int));

  /* initialize the pruning method data array. */
  E->prune_data = (void***) malloc(E->G->n_order * sizeof(void**));

  /* check if any allocation failed. */
  if (!E->prune || !E->prune_sz || !E->prune_data)
    throw("unable to allocate pruning arrays");

  /* initialize the inner arrays. */
  for (i = 0; i < E->G->n_order; i++) {
    E->prune[i] = NULL;
    E->prune_sz[i] = 0;
    E->prune_data[i] = NULL;
  }

  /* loop over the pruning method names in the options structure. */
  for (i = 0; i < opts->n_prune; i++) {
    /* add the pruning method to the enumerator. */
    if (!enum_init_prune_add(E, opts->prune[i]))
      return 0;
  }

  /* return success. */
  return 1;
}

/* enum_new(): allocate a new enumerator data structure.
 *
 * arguments:
 *  @P: pointer to a related peptide data structure to access.
 *  @G: pointer to a related graph data structure to access.
 *  @opts: pointer to an options data structure to access.
 *
 * returns:
 *  pointer to a newly allocated and initialized enumerator structure, or
 *  NULL on failure.
 */
enum_t *enum_new (peptide_t *P, graph_t *G, opts_t *opts) {
  /* declare required variables:
   *  @E: output structure pointer.
   */
  enum_t *E;

  /* allocate a new structure pointer. */
  E = (enum_t*) malloc(sizeof(enum_t));

  /* check if allocation failed. */
  if (!E) {
    /* raise an exception and return null. */
    raise("unable to allocate enumerator structure pointer");
    return NULL;
  }

  /* store the related data structures. */
  E->P = P;
  E->G = G;

  /* initialize the output system variables. */
  E->fd = -1;
  E->nsol = 0;
  E->ntree = 0.0;
  E->nmax = opts->nsol_limit;
  E->fname = strdup(opts->fname_out);

  /* initialize the termination variable. */
  E->term = 0;

  /* store the branching control variables. */
  E->nbmax = opts->branch_max / 2;
  E->eps = opts->branch_eps;

  /* initialize the threads. */
  if (!enum_init_threads(E, opts)) {
    /* raise an exception and return null. */
    raise("unable to initialize threads");
    enum_free(E);
    return NULL;
  }

  /* store the pruning control variables. */
  E->ddf_tol = opts->ddf_tol;
  E->rmsd_tol = (double) G->n_orig * pow(opts->rmsd_tol, 2.0);

  /* set the enumerator output format. */
  if (!enum_init_format(E, opts)) {
    /* raise an exception and return null. */
    raise("unable to set output format");
    enum_free(E);
    return NULL;
  }

  /* initialize the pruning methods. */
  if (!enum_init_prune(E, opts)) {
    /* raise an exception and return null. */
    raise("unable to initialize pruning methods");
    enum_free(E);
    return NULL;
  }

  /* return the new structure pointer. */
  return E;
}

/* enum_free(): free all allocated memory associated with an enumerator.
 *
 * arguments:
 *  @E: pointer to the enumerator structure to free.
 */
void enum_free (enum_t *E) {
  /* declare required variables:
   *  @i, @j: general-purpose loop counters.
   */
  unsigned int i, j;

  /* return if the structure pointer is null. */
  if (!E) return;

#ifdef __IBP_HAVE_PTHREAD
  /* destroy the write mutex. */
  pthread_mutex_destroy(&E->write_mutex);
#endif

  /* cleanup the output system. */
  if (E->write_close)
    E->write_close(E);

  /* free the directory name string. */
  free(E->fname);

  /* free the pruning function array. */
  if (E->prune) {
    /* free the inner array elements. */
    for (i = 0; i < E->G->n_order; i++)
      free(E->prune[i]);

    /* free the outer array. */
    free(E->prune);
  }

  /* free the pruning data array. */
  if (E->prune_data) {
    /* free the inner array elements. */
    for (i = 0; i < E->G->n_order; i++) {
      /* free the inner array elements. */
      for (j = 0; j < E->prune_sz[i]; j++)
        free(E->prune_data[i][j]);

      /* free the outer array. */
      free(E->prune_data[i]);
    }

    /* free the outer array. */
    free(E->prune_data);
  }

  /* free the pruning test sizes. */
  free(E->prune_sz);

  /* free the threads. */
  free(E->threads);

  /* finally, free the structure pointer. */
  free(E);
}

/* enum_report(): output a pruning report after enumeration.
 *
 * arguments:
 *  @E: pointer to the enumerator structure to access.
 */
static void enum_report (enum_t *E) {
  /* declare required variables:
   *  @i: pruning array index at each level.
   *  @m: pruning method index.
   *  @lev: reorder level index.
   *  @testfn: pruning test function pointer.
   *  @reportfn: pruning report function pointer.
   */
  unsigned int i, m, lev;
  enum_prune_test_fn testfn;
  enum_prune_report_fn reportfn;

  /* loop over the pruning methods. */
  for (m = 0; pruners[m].name; m++) {
    /* output an initial header. */
    printf("\nPruning results [%s]:\n", pruners[m].name);

    /* get the pruning function pointers. */
    testfn = pruners[m].prune_test;
    reportfn = pruners[m].prune_report;

    /* loop over the levels of the graph order. */
    for (lev = 0; lev < E->G->n_order; lev++) {
      /* skip duplicate levels. */
      if (E->G->orig[lev]) continue;

      /* loop over the registered pruning closures for the current level. */
      for (i = 0; i < E->prune_sz[lev]; i++) {
        /* skip closures that do not match the current pruning method. */
        if (E->prune[lev][i] != testfn) continue;

        /* print reporting information for the closure. */
        reportfn(E, lev, E->prune_data[lev][i]);
      }
    }
  }
}

/* enum_execute(): enumerate all solutions from an iDMDGP graph/peptide
 * structure pair using an enumerator.
 *
 * arguments:
 *  @E: pointer to the enumerator structure to utilize.
 *
 * returns:
 *  integer indicating the number of enumerated solutions, or zero on
 *  failure.
 */
int enum_execute (enum_t *E) {
  /* initialize the solution count and open the output system. */
  E->nsol = 0;
  if (E->write_open && !E->write_open(E))
    throw("unable to open enumerator output");

  /* initialize the threads for enumeration. */
  if (!enum_threads_init(E))
    throw("unable to initialize enumerator threads");

#if defined(__IBP_HAVE_PTHREAD)
#if defined(__IBP_HAVE_CUDA)

  /* FIXME: handle gpu thread kickoff here. */

#endif /* __IBP_HAVE_CUDA */

  /* execute the threads. */
  for (unsigned int i = 0; i < E->nthreads; i++) {
    /* create the thread. */
    int ret = pthread_create(&E->threads[i].thread, NULL,
                             enum_thread_execute,
                             (void*) (E->threads + i));

    /* check the thread creation result. */
    if (ret)
      throw("unable to create thread %u of %u", i + 1, E->nthreads);
  }

  /* wait for all threads to exit. */
  for (unsigned int i = 0; i < E->nthreads; i++)
    pthread_join(E->threads[i].thread, NULL);

#else /* __IBP_HAVE_PTHREAD */

  /* execute a single enumerator in the current thread. boring. */
  enum_thread_execute((void*) E->threads);\

#endif /* __IBP_HAVE_PTHREAD */

  /* close the output system. */
  if (E->write_close)
    E->write_close(E);

  /* write the pruning report. */
  enum_report(E);

  /* return success. */
  return 1;
}

