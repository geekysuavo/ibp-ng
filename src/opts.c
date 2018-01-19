
/* include the option parsing header. */
#include "opts.h"

/* define all accepted short options.
 */
#define OPTS_S_HELP        'h'
#define OPTS_S_VERBOSE     'v'
#define OPTS_S_INPUT       'i'
#define OPTS_S_CHAIN       'c'
#define OPTS_S_SEQID       'n'
#define OPTS_S_OUTPUT      'o'
#define OPTS_S_PSF         'p'
#define OPTS_S_DMDGP       'd'
#define OPTS_S_FORMAT      'f'
#define OPTS_S_RESTRAINT   'r'
#define OPTS_S_SIDECHAIN   's'
#define OPTS_S_TOPOLOGY    'T'
#define OPTS_S_PARAMS      'P'
#define OPTS_S_REORDER     'R'
#define OPTS_S_GPU         'g'
#define OPTS_S_THREADS     't'
#define OPTS_S_METHOD      'm'
#define OPTS_S_BRANCH_MAX  'b'
#define OPTS_S_BRANCH_EPS  'e'
#define OPTS_S_LIMIT       'l'
#define OPTS_S_COMPLETE   ('z'+1)
#define OPTS_S_VDW_SCALE  ('z'+2)
#define OPTS_S_DDF_TOL    ('z'+3)
#define OPTS_S_RMSD       ('z'+4)

/* define all accepted long options.
 */
#define OPTS_L_HELP       "help"
#define OPTS_L_VERBOSE    "verbose"
#define OPTS_L_INPUT      "input"
#define OPTS_L_CHAIN      "chain"
#define OPTS_L_SEQID      "seqid"
#define OPTS_L_OUTPUT     "output"
#define OPTS_L_PSF        "psf"
#define OPTS_L_DMDGP      "dmdgp"
#define OPTS_L_FORMAT     "format"
#define OPTS_L_RESTRAINT  "restraints"
#define OPTS_L_SIDECHAIN  "sidechain"
#define OPTS_L_TOPOLOGY   "topology"
#define OPTS_L_PARAMS     "params"
#define OPTS_L_REORDER    "reorder"
#define OPTS_L_GPU        "gpu"
#define OPTS_L_THREADS    "threads"
#define OPTS_L_METHOD     "method"
#define OPTS_L_BRANCH_MAX "branch-max"
#define OPTS_L_BRANCH_EPS "branch-eps"
#define OPTS_L_LIMIT      "limit"
#define OPTS_L_COMPLETE   "complete"
#define OPTS_L_VDW_SCALE  "vdw-scale"
#define OPTS_L_DDF_TOL    "ddf-tol"
#define OPTS_L_RMSD       "rmsd"

/* opts_config_t: option definition structure for informing opts_next()
 * about all supported command line options that the user may specify.
 */
typedef struct {
  /* @lname: long option name string.
   * @sname: short option name character.
   */
  const char *lname;
  int sname;

  /* @has_arg: whether or not this option requires an argument. */
  int has_arg;
}
opts_config_t;

/* declare the option definition array used by opts_next for reading
 * command line arguments specified by the user.
 */
static const opts_config_t conf[] = {
  { OPTS_L_HELP,       OPTS_S_HELP,       0 },
  { OPTS_L_VERBOSE,    OPTS_S_VERBOSE,    0 },
  { OPTS_L_INPUT,      OPTS_S_INPUT,      1 },
  { OPTS_L_CHAIN,      OPTS_S_CHAIN,      1 },
  { OPTS_L_SEQID,      OPTS_S_SEQID,      1 },
  { OPTS_L_OUTPUT,     OPTS_S_OUTPUT,     1 },
  { OPTS_L_PSF,        OPTS_S_PSF,        1 },
  { OPTS_L_DMDGP,      OPTS_S_DMDGP,      1 },
  { OPTS_L_FORMAT,     OPTS_S_FORMAT,     1 },
  { OPTS_L_RESTRAINT,  OPTS_S_RESTRAINT,  1 },
  { OPTS_L_SIDECHAIN,  OPTS_S_SIDECHAIN,  1 },
  { OPTS_L_TOPOLOGY,   OPTS_S_TOPOLOGY,   1 },
  { OPTS_L_PARAMS,     OPTS_S_PARAMS,     1 },
  { OPTS_L_REORDER,    OPTS_S_REORDER,    1 },
  { OPTS_L_GPU,        OPTS_S_GPU,        0 },
  { OPTS_L_THREADS,    OPTS_S_THREADS,    1 },
  { OPTS_L_METHOD,     OPTS_S_METHOD,     1 },
  { OPTS_L_BRANCH_MAX, OPTS_S_BRANCH_MAX, 1 },
  { OPTS_L_BRANCH_EPS, OPTS_S_BRANCH_EPS, 1 },
  { OPTS_L_LIMIT,      OPTS_S_LIMIT,      1 },
  { OPTS_L_COMPLETE,   OPTS_S_COMPLETE,   0 },
  { OPTS_L_VDW_SCALE,  OPTS_S_VDW_SCALE,  1 },
  { OPTS_L_DDF_TOL,    OPTS_S_DDF_TOL,    1 },
  { OPTS_L_RMSD,       OPTS_S_RMSD,       1 },

  /* null terminator. */
  { NULL,              '\0',              0 }
};

/* opts_next(): return the next parsed option in an argument array.
 *
 * this function is patterned after the getopt style of option parsing,
 * but does not allow for grouping single-character non-argument options
 * into a single chunk, like getopt does.
 *
 * arguments:
 *  @argc: argument count from main().
 *  @argv: argument array from main().
 *  @conf: array of accepted command line arguments.
 *  @argi: pointer to the current parsing position.
 *
 * returns:
 *  integer containing the short-form of the next parsed option in the
 *  @argv array, or 0 if an invalid option was passed, or -1 if no further
 *  options are available for parsing in @argv.
 */
int opts_next (int argc, char **argv, const opts_config_t *conf, int *argi) {
  /* declare required variables:
   *  @arg: shorthand storage location for the current argument string.
   *  @i: loop counter for searching the option definition array.
   */
  const char *arg;
  int i;

  /* ensure parsing begins at the correct index. */
  if (*argi < 1)
    *argi = 1;

  /* check that the argument index is in bounds. */
  if (*argi >= argc)
    return -1;

  /* get the currently indexed argument. */
  arg = argv[*argi];

  /* check the currently indexed option. */
  if (strlen(arg) == 2 && arg[0] == '-') {
    /* short option: search the option definition array. */
    for (i = 0; conf[i].lname; i++) {
      /* check if the current option is a match. */
      if (conf[i].sname == arg[1])
        break;
    }
  }
  else if (strlen(arg) >= 3 && arg[0] == '-' && arg[1] == '-') {
    /* long option: search the option definition array. */
    for (i = 0; conf[i].lname; i++) {
      /* check if the current option is a match. */
      if (strcmp(conf[i].lname, arg + 2) == 0)
        break;
    }
  }
  else {
    /* throw an exception. */
    throw("invalid argument '%s'", arg);
  }

  /* check if an option was identified. */
  if (conf[i].lname) {
    /* yes. check if an option argument is required. */
    if (conf[i].has_arg) {
      /* make sure that an extra argument exists. */
      if (*argi >= argc - 1)
        throw("option '--%s' requires an argument", conf[i].lname);
    }

    /* everything looks good. increment the argument index and return. */
    (*argi)++;
    return conf[i].sname;
  }
  else {
    /* throw an exception. */
    throw("invalid option '%s'", arg);
  }

  /* indicate that no arguments remain. */
  return -1;
}

/* opts_new(): allocate and initialize an empty options data structure.
 *
 * returns:
 *  newly allocated pointer to an empty options data structure.
 */
opts_t *opts_new (void) {
  /* declare required variables:
   *  @opts: pointer to the output options structure.
   */
  opts_t *opts;

  /* allocate a new structure pointer. */
  opts = (opts_t*) malloc(sizeof(opts_t));

  /* check if allocation failed. */
  if (!opts) {
    /* raise an exception and return null. */
    raise("unable to allocate options structure pointer");
    return NULL;
  }

  /* initialize the status fields of the data structure. */
  opts->help = 0;
  opts->verb = 0;

  /* initialize the filename fields. */
  opts->fname_in = NULL;
  opts->fname_out = NULL;
  opts->fname_top = NULL;
  opts->fname_par = NULL;
  opts->fname_ord = NULL;
  opts->fname_psf = NULL;
  opts->fname_dmdgp = NULL;

  /* initialize the file option fields. */
  opts->idx_in = NULL;
  opts->fmt_out = NULL;

  /* initialize the restraint filename fields. */
  opts->fname_restr = NULL;
  opts->n_restr = 0;

  /* initialize the sidechain fields. */
  opts->sidech = NULL;
  opts->n_sidech = 0;

  /* initialize the pruning fields. */
  opts->prune = NULL;
  opts->n_prune = 0;

  /* initialize branch control fields. */
  opts->thread_gpu = 0;
  opts->thread_num = 1;
  opts->branch_max = 20;
  opts->branch_eps = 0.05;

  /* initialize prune control fields. */
  opts->nsol_limit = 0;
  opts->complete = 0;
  opts->vdw_scale = 0.6;
  opts->ddf_tol = 0.001;
  opts->rmsd_tol = 0.0;

  /* return the new options data structure. */
  return opts;
}

/* opts_add_restraints(): append a new restraints filename to the
 * appropriate array of an options data structure.
 *
 * arguments:
 *  @opts: pointer to the options structure to modify.
 *  @fname: restraints filename to append.
 *
 * returns:
 *  integer indicating whether (1) or not (0) the operation succeeded.
 */
int opts_add_restraints (opts_t *opts, char *fname) {
  /* check that the structure pointer is valid. */
  if (!opts)
    throw("options structure pointer is invalid");

  /* increment the restraint filename count. */
  opts->n_restr++;

  /* reallocate the restraint filename array. */
  opts->fname_restr = (char**)
    realloc(opts->fname_restr, opts->n_restr * sizeof(char*));

  /* check if reallocation failed. */
  if (!opts->fname_restr)
    throw("unable to reallocate restraint filename array");

  /* store the new restraint filename. */
  opts->fname_restr[opts->n_restr - 1] = fname;

  /* return success. */
  return 1;
}

/* opts_add_sidechains(): append a new set of sidechain indices to the
 * appropriate array of an options data structure.
 *
 * arguments:
 *  @opts: pointer to the options structure to modify.
 *  @str: sidechain string to parse and append.
 *
 * returns:
 *  integer indicating whether (1) or not (0) the operation succeeded.
 */
int opts_add_sidechains (opts_t *opts, char *str) {
  /* declare required variables:
   *  @pa, @pb: start and end substring pointers.
   *  @idx: currently parsed sidechain index.
   *  @done: completion flag.
   */
  char *pa, *pb;
  int idx, done;

  /* check that the structure pointer is valid. */
  if (!opts)
    throw("options structure pointer is invalid");

  /* initialize the start pointer and completion flag. */
  pa = str;
  done = 0;

  /* loop until completion. */
  do {
    /* find the next integer. */
    pb = strstr(pa, ",");

    /* check if a new pointer was found. */
    if (!pb) {
      /* no. prepare the last string. */
      done = 1;
      pb = pa + strlen(pa);
    }

    /* skip empty strings. */
    if (pb - pa > 0) {
      /* parse the current integer. */
      idx = atoi(pa);

      /* check that the integer is in bounds. */
      if (idx < 1)
        throw("sidechain index %d out of bounds [1,inf)", idx);

      /* increment the sidechain count. */
      opts->n_sidech++;

      /* reallocate the sidechain index array. */
      opts->sidech = (unsigned int*)
        realloc(opts->sidech, opts->n_sidech * sizeof(unsigned int));

      /* check if reallocation failed. */
      if (!opts->sidech)
        throw("unable to reallocate sidechain array");

      /* store the new sidechain index. */
      opts->sidech[opts->n_sidech - 1] = idx;
    }

    /* move past the delimiter. */
    pa = pb + 1;
  }
  while (!done);

  /* return success. */
  return 1;
}

/* opts_add_pruner(): append a new set of pruning device names to the
 * appropriate array of an options data structure.
 *
 * arguments:
 *  @opts: pointer to the options structure to modify.
 *  @str: pruning device string to parse and append.
 *
 * returns:
 *  integer indicating whether (1) or not (0) the operation succeeded.
 */
int opts_add_pruner (opts_t *opts, char *str) {
  /* declare required variables:
   *  @pa, @pb: start and end substring pointers.
   *  @done: completion flag.
   */
  char *pa, *pb;
  int done;

  /* check that the structure pointer is valid. */
  if (!opts)
    throw("options structure pointer is invalid");

  /* initialize the start pointer and completion flag. */
  pa = str;
  done = 0;

  /* loop until completion. */
  do {
    /* find the next integer. */
    pb = strstr(pa, ",");

    /* check if a new pointer was found. */
    if (!pb) {
      /* no. prepare the last string. */
      done = 1;
      pb = pa + strlen(pa);
    }

    /* skip empty strings. */
    if (pb - pa > 0) {
      /* increment the pruner count. */
      opts->n_prune++;

      /* reallocate the pruner array. */
      opts->prune = (char**)
        realloc(opts->prune, opts->n_prune * sizeof(char*));

      /* check if reallocation failed. */
      if (!opts->prune)
        throw("unable to reallocate pruner name array");

      /* store the new pruner name. */
      opts->prune[opts->n_prune - 1] = strndup(pa, pb - pa);
    }

    /* move past the delimiter. */
    pa = pb + 1;
  }
  while (!done);

  /* return success. */
  return 1;
}

/* opts_new_from_strings(): allocate and fill a new options data structure
 * based on an array of command line argument strings.
 *
 * arguments:
 *  @argc: number of command line argument strings.
 *  @argv: array of command line argument strings.
 *
 * returns:
 *  newly allocated pointer to a filled options data structure, or null
 *  on failure.
 */
opts_t *opts_new_from_strings (int argc, char **argv) {
  /* declare variables for storing parsed options:
   *  @opts: pointer to the output options structure.
   */
  opts_t *opts;

  /* declare variables used for command line argument parsing:
   *  @argi: current argument array index in opts_next().
   *  @c: returned option character from opts_next().
   */
  int argi = 0;
  int c;

  /* initialize the options data structure. */
  opts = opts_new();
  if (!opts)
    return NULL;

  /* loop until the arguments are exhausted. */
  while ((c = opts_next(argc, argv, conf, &argi)) != -1) {
    /* determine which option was specified. */
    switch (c) {
      /* help mode. */
      case OPTS_S_HELP:
        /* increase the help message level. */
        opts->help++;
        break;

      /* verbose mode. */
      case OPTS_S_VERBOSE:
        /* increase the verbosity level. */
        opts->verb++;
        break;

      /* input filename. */
      case OPTS_S_INPUT:
        /* set the input filename. */
        opts->fname_in = argv[argi];
        argi++;
        break;

      /* output filename. */
      case OPTS_S_OUTPUT:
        /* set the output filename. */
        opts->fname_out = argv[argi];
        argi++;
        break;

      /* psf output filename. */
      case OPTS_S_PSF:
        opts->fname_psf = argv[argi];
        argi++;
        break;

      /* dmdgp output filename. */
      case OPTS_S_DMDGP:
        /* set the dmdgp filename. */
        opts->fname_dmdgp = argv[argi];
        argi++;
        break;

      /* chain identifier.
       * sequence number.
       */
      case OPTS_S_CHAIN:
      case OPTS_S_SEQID:
        /* set the input file option string. */
        opts->idx_in = argv[argi];
        argi++;
        break;

      /* output format. */
      case OPTS_S_FORMAT:
        /* set the output file option string. */
        opts->fmt_out = argv[argi];
        argi++;
        break;

      /* restraints filename. */
      case OPTS_S_RESTRAINT:
        /* add the new restraint filename. */
        if (!opts_add_restraints(opts, argv[argi])) {
          /* raise an exception and return null. */
          raise("unable to add restraint filename");
          opts_free(opts);
          return NULL;
        }

        /* increment the argument index and break. */
        argi++;
        break;

      /* sidechain index or index list. */
      case OPTS_S_SIDECHAIN:
        /* add the new sidechain index or index list. */
        if (!opts_add_sidechains(opts, argv[argi])) {
          /* raise an exception and return null. */
          raise("unable to add sidechain index/indices");
          opts_free(opts);
          return NULL;
        }

        /* increment the argument index and break. */
        argi++;
        break;

      /* topology filename. */
      case OPTS_S_TOPOLOGY:
        /* set the topology filename. */
        opts->fname_top = argv[argi];
        argi++;
        break;

      /* parameter filename. */
      case OPTS_S_PARAMS:
        /* set the parameter filename. */
        opts->fname_par = argv[argi];
        argi++;
        break;

      /* reorder filename. */
      case OPTS_S_REORDER:
        /* set the reorder filename. */
        opts->fname_ord = argv[argi];
        argi++;
        break;

      /* number of threads. */
      case OPTS_S_THREADS:
#if !defined(__IBP_HAVE_PTHREAD)
        /* raise an exception about thread support. */
        raise("this program was compiled without thread support");
        opts_free(opts);
        return NULL;
#else
        /* set the thread count. */
        opts->thread_num = atoi(argv[argi]);
        argi++;
        break;
#endif

      /* gpu usage flag. */
      case OPTS_S_GPU:
#if !defined(__IBP_HAVE_PTHREAD) || !defined(__IBP_HAVE_CUDA)
        /* raise an exception about thread support. */
        raise("this program was compiled without gpu support");
        opts_free(opts);
        return NULL;
#else
        /* set the gpu flag. */
        opts->thread_gpu++;
        break;
#endif

      /* pruning method. */
      case OPTS_S_METHOD:
        /* add the new pruning method or method list. */
        if (!opts_add_pruner(opts, argv[argi])) {
          /* raise an exception and return null. */
          raise("unable to add pruning method name/names");
          opts_free(opts);
          return NULL;
        }

        /* increment the argument index and break. */
        argi++;
        break;

      /* branch maximum. */
      case OPTS_S_BRANCH_MAX:
        opts->branch_max = atoi(argv[argi]);
        argi++;
        break;

      /* branch epsilon. */
      case OPTS_S_BRANCH_EPS:
        opts->branch_eps = atof(argv[argi]);
        argi++;
        break;

      /* solution limit. */
      case OPTS_S_LIMIT:
        opts->nsol_limit = atoi(argv[argi]);
        argi++;
        break;

      /* graph completion flag. */
      case OPTS_S_COMPLETE:
        opts->complete++;
        break;

      /* vdw scale factor. */
      case OPTS_S_VDW_SCALE:
        opts->vdw_scale = atof(argv[argi]);
        argi++;
        break;

      /* ddf error tolerance. */
      case OPTS_S_DDF_TOL:
        opts->ddf_tol = atof(argv[argi]);
        argi++;
        break;

      /* rmsd tolerance. */
      case OPTS_S_RMSD:
        opts->rmsd_tol = atof(argv[argi]);
        argi++;
        break;

      /* unknown option or argument. */
      default:
        raise("failed to parse arguments");
        opts_free(opts);
        return NULL;
    }
  }

  /* set the verbosity level. */
  verbosity_set(opts->verb);

  /* return the newly initialized options structure pointer. */
  return opts;
}

/* opts_free(): free all allocated memory belonging to an options data
 * structure.
 *
 * arguments:
 *  @opts: pointer to the options data structure to free.
 */
void opts_free (opts_t *opts) {
  /* declare required variables:
   *  @i: general-purpose loop counter.
   */
  unsigned int i;

  /* return if the structure pointer is null. */
  if (!opts) return;

  /* free the array of restraint filenames. */
  if (opts->n_restr)
    free(opts->fname_restr);

  /* free the array of sidechain indices. */
  if (opts->n_sidech)
    free(opts->sidech);

  /* free the array of pruning device names. */
  if (opts->n_prune) {
    /* free the array elements. */
    for (i = 0; i < opts->n_prune; i++) {
      if (opts->prune[i])
        free(opts->prune[i]);
    }

    /* free the array. */
    free(opts->prune);
  }

  /* free the structure pointer. */
  free(opts);
}

/* opts_validate(): performs a few simple checks on the contents of
 * an options data structure to make sure that later operations will
 * have sane inputs.
 *
 * arguments:
 *  @opts: pointer to the options structure to validate.
 *
 * returns:
 *  integer indicating whether (1) or not (0) validation succeeded.
 */
int opts_validate (opts_t *opts) {
  /* check that an input filename was specified. */
  if (!opts->fname_in)
    throw("expected input filename not specified");

  /* check that an output filename was specified. */
  if (!opts->fname_out) {
    /* if not, build a directory from the input filename. */
    opts->fname_out = (char*)
      malloc((strlen(opts->fname_in) + 16) * sizeof(char));

    /* if the allocation failed, just return failure. */
    if (!opts->fname_out)
      throw("expected output filename not specified");

    /* construct the output filename. */
    sprintf(opts->fname_out, "%s.d", opts->fname_in);

    /* but make sure to output a warning as well. */
    warn("output filename unspecified, defaulting to '%s'",
         opts->fname_out);
  }

  /* check that a topology filename was specified. */
  if (!opts->fname_top)
    raise("expected topology filename not specified");

  /* check that a parameter filename was specified. */
  if (!opts->fname_par)
    raise("expected parameter filename not specified");

  /* check that a reorder filename was specified. */
  if (!opts->fname_ord)
    raise("expected reorder filename not specified");

  /* validate the thread count. */
  if (opts->thread_num == 0)
    raise("thread count must be non-zero");

  /* validate the branch count. */
  if (opts->branch_max == 0)
    raise("maximum branch count must be positive");

  /* validate the branch epsilon. */
  if (opts->branch_eps <= 0.0)
    raise("branch epsilon out must be positive");

  /* validate the vdw scaling factor. */
  if (opts->vdw_scale <= 0.0)
    raise("vdw scale factor must be positive");

  /* validate the ddf tolerance. */
  if (opts->ddf_tol < 0.0)
    raise("DDF: error tolerance must be non-negative");

  /* return valid. */
  return (traceback_length() == 0);
}

