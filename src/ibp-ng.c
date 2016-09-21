
/* include the main header. */
#include "ibp-ng.h"

/* IBPNG_HELPSTR: short string that is displayed when the user specifies
 * no arguments or passes the --help argument.
 */
#define IBPNG_HELPSTR "\
 ibp-ng: Interval Branch and Prune utility for protein structures.\n\
 Copyright (C) 2016 Institut Pasteur. All rights reserved.\n\
\n\
 Usage:\n\
  ibp-ng [OPTIONS]\n\
\n\
 Program output options:\n\
  -h, --help              Display this help message                   [off]\n\
  -v, --verbose           Enable verbose messages                     [off]\n\
\n\
 Input file options:\n\
  -i, --input FIN         Input filename\n\
  -c, --chain CID         PDB chain identifier\n\
  -n, --seqid SID         FASTA sequence index\n\
\n\
 Output file options:\n\
  -p, --psf FPSF          Output PSF file                            [none]\n\
  -d, --dmdgp FDMD        Output DMDGP file                          [none]\n\
  -o, --output FOUT       Output filename                            [auto]\n\
  -f, --format FMT        Output format                               [dcd]\n\
  -r, --restraints RES    Input restraints filename                  [none]\n\
  -s, --sidechain SL      Sidechain(s) to make explicit              [none]\n\
\n\
 Configuration file options:\n\
  -T, --topology TOP      Topology file for protein geometry\n\
  -P, --params PAR        Parameter file for protein geometry\n\
  -R, --reorder ORD       Repetition order definitions filename\n\
\n\
 Enumeration options:\n\
  -m, --method ML         Pruning method(s) to use                   [none]\n\
  -b, --branch-max NB     Maximum number of branches per node          [20]\n\
  -e, --branch-eps EPS    Minimum interval discretization            [0.05]\n\
  -l, --limit NSOL        Maximum number of solutions                 [off]\n\
      --vdw-scale VF      Atomic radius scaling factor                [0.6]\n\
      --ddf-tol TOL       DDF error tolerance                       [0.001]\n\
\n\
 Parallel execution options:\n\
  -g, --gpu               Flag to execute on the GPU                  [off]\n\
  -t, --threads NT        Number of threads to execute                  [1]\n\
\n\
 The ibp-ng utility enumerates all feasible solutions to a given Interval\n\
 Discretizable Molecular Distance Geometry Problem (iDMDGP) instance, or\n\
 more simply, produces all protein structures that agree with a combined\n\
 set of local and non-local distance and angle constraints.\n\
\n\
"

/* declare a set of custom structures:
 *  @opts: options data structure for storing command line arguments.
 *  @ord: repetition order structure for constructing graphs.
 *  @top: molecular topology structure for constructing graphs.
 *  @par: molecular parameters structure for constructing graphs.
 *  @P: peptide structure for initializing molecular information.
 *  @G: graph structure for summarizing a problem instance.
 *  @E: enum structure for enumerating all feasible solutions.
 */
opts_t *opts = NULL;
reorder_t *ord = NULL;
topol_t *top = NULL;
param_t *par = NULL;
peptide_t *P = NULL;
graph_t *G = NULL;
enum_t *E = NULL;

/* main_handler(): handle interrupt signals in order to properly clean
 * up any running enumeration threads.
 *
 * arguments:
 *  @sig: code of the caught signal.
 */
void main_handler (int sig) {
  /* write an informational message and raise the termination flag. */
  info("attempting a clean getaway...");
  E->term = 1;
}

/* main(): application entry point.
 *
 * arguments:
 *  @argc: number of command line arguments.
 *  @argv: array of command line arguments.
 *
 * returns:
 *  integer indicating application success (0) or failure (!0).
 */
int main (int argc, char **argv) {
  /* declare other required variables:
   *  @i: general purpose loop counter.
   */
  unsigned int i;

  /* check if exactly zero arguments were provided. */
  if (argc == 1) {
    /* print the help message string and return without an error. */
    fprintf(stdout, IBPNG_HELPSTR);
    return 0;
  }

  /* parse the command line arguments. */
  opts = opts_new_from_strings(argc, argv);

  /* check that command line argument parsing succeeded. */
  if (!opts)
    die("unable to parse command line arguments");

  /* check if the help flag was raised. */
  if (opts->help) {
    /* print the help message string and return without an error. */
    fprintf(stdout, IBPNG_HELPSTR);
    return 0;
  }

  /* validate the option values. */
  if (!opts_validate(opts))
    die("one or more invalid program arguments");

  /* create a new topology structure from the specified file. */
  top = topol_new_from_file(opts->fname_top);

  /* check that the topology structure was successfully created. */
  if (!top)
    die("unable to read topology from '%s'", opts->fname_top);

  /* create a new parameter structure from the specified file. */
  par = param_new_from_file(opts->fname_par, opts->vdw_scale);

  /* check that the parameter structure was successfully created. */
  if (!par)
    die("unable to read parameters from '%s'", opts->fname_par);

  /* create a new reorder structure from the specified file. */
  ord = reorder_new_from_file(opts->fname_ord);

  /* check that the reorder structure was successfully created. */
  if (!ord)
    die("unable to read reorder from '%s'", opts->fname_ord);

  /* create a new peptide structure from an input file. */
  P = peptide_new_from_file(opts->fname_in, opts->idx_in);

  /* check that the peptide structure was successfully created. */
  if (!P)
    die("unable to read initial sequence from '%s'", opts->fname_in);

  /* make all requested peptide residue sidechains explicit. */
  for (i = 0; i < opts->n_sidech; i++) {
    /* make the current sidechain explicit. */
    if (!peptide_add_sidechain(P, opts->sidech[i] - 1))
      die("unable to make sidechain %u explicit", opts->sidech[i]);
  }

  /* apply topology information to the peptide structure. */
  if (!topol_apply_all(top, P))
    die("unable to apply peptide topology");

  /* apply parameter information to the peptide structure. */
  if (!param_apply_all(par, P))
    die("unable to apply peptide parameters");

  /* loop over the restraint filenames. */
  for (i = 0; i < opts->n_restr; i++) {
    /* feed the restraint file to the graph. */
    if (!assign_set_from_file(P, opts->fname_restr[i]))
      die("unable to add restraints from '%s'", opts->fname_restr[i]);
  }

  /* back-calculate the peptide force-field / probability parameters. */
  if (!peptide_field(P, opts->ddf_tol))
    die("unable to recompute force field parameters");

  /* create a graph structure from the peptide information. */
  G = peptide_graph(P, ord);

  /* check if graph creation failed. */
  if (!G)
    die("failed to build peptide graph");

  /* check if a psf output file was requested. */
  if (opts->fname_psf) {
    /* attempt to write the intermediate output file. */
    if (!psf_write(opts->fname_psf, P, G))
      die("unable to write peptide data to '%s'", opts->fname_psf);
  }

  /* check if a dmdgp output file was requested. */
  if (opts->fname_dmdgp) {
    /* attempt to write the intermediate output file. */
    if (!dmdgp_write(opts->fname_dmdgp, P, G))
      die("unable to write graph data to '%s'", opts->fname_dmdgp)
  }

  /* create an enumerator structure for computing all feasible solutions. */
  E = enum_new(P, G, opts);

  /* check if enumerator creation failed. */
  if (!E)
    die("failed to build graph enumerator");

  /* catch interrupt signals. */
  signal(SIGINT, main_handler);

  /* enumerate all solutions from the graph. */
  if (!enum_execute(E))
    die("failed to enumerate graph solutions");

/* death: label used by all die() macro functions to cleanly
 * terminate application execution without leaving allocated
 * memory on the heap.
 */
death:
  /* free the options structure. */
  opts_free(opts);

  /* free the molecular information structures. */
  reorder_free(ord);
  topol_free(top);
  param_free(par);

  /* free the enumerator structure. */
  enum_free(E);

  /* free the peptide and graph structures. */
  peptide_free(P);
  graph_free(G);

  /* check if the traceback contains entries. */
  if (traceback_length()) {
    /* clean up the traceback array and return failure. */
    traceback_clear();
    return 1;
  }

  /* return success. */
  return 0;
}

