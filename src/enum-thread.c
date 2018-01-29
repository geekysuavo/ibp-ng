
/* include the enumerator headers. */
#include "enum.h"
#include "enum-thread.h"
#include "enum-reduce.h"
#include "enum-write.h"

/* state_div(): divide the size of a thread state into uniform pieces,
 * storing the result into the index of the thread state.
 *
 * arguments:
 *  @state: thread state to operate on.
 *  @len: number of elements in the state.
 *  @n: value to divide the state size by.
 */
static void state_div (enum_thread_node_t *state,
                       unsigned int len,
                       unsigned int n) {
  /* @bN: base-2 log of the tree size/leaf count.
   * @bn: base-2 log of the division value.
   * @bs: base-2 log of the quotient.
   * @i: general-purpose loop counter.
   */
  double bN, bn, bs;
  unsigned int i;

  /* compute the tree size. */
  for (i = 0, bN = 0.0; i < len; i++)
    bN += log2((double) state[i].nb);

  /* compute the quotient. */
  bn = bN - log2((double) n);

  /* loop to compute the division result. */
  for (i = len - 1, bs = 0.0; i < len; i--) {
    /* sum the current node size into the "stride". */
    bs += log2((double) state[i].nb);

    /* check if the "stride" is sufficiently large. */
    if (bs == bn) {
      /* equal: amounts to one step of the parent node index. */
      state[i - 1].end = 1;
      break;
    }
    else if (bs > bn) {
      /* larger: compute the closest integer step amount based
       * on the required amount of bits in the quotient.
       */
      bs -= log2((double) state[i].nb);
      state[i].end = 1 << (unsigned int) floor(bn - bs);
      break;
    }
  }
}

/* state_add(): add the indices of two states, storing the result into
 * the index of the destination (first) state.
 *
 * arguments:
 *  @dest: destination state in the sum.
 *  @src1: source state in the sum.
 *  @len: number of elements in each state.
 */
static void state_add (enum_thread_node_t *dest,
                       enum_thread_node_t *src,
                       unsigned int len) {
  /* initialize the output values. */
  for (unsigned int i = 0; i < len; i++)
    dest[i].end = 0;

  /* loop over each index of the result. */
  for (unsigned int i = len - 1; i < len; i--) {
    /* compute the current term of the sum. */
    dest[i].end += dest[i].idx + src[i].end;

    /* perform any required carries. */
    while (i && dest[i].nb && dest[i].end >= dest[i].nb) {
      dest[i].end -= dest[i].nb;
      dest[i - 1].end++;
    }
  }
}

/* state_valid(): check whether a state is "valid", meaning that its index
 * has not yet reached its end.
 *
 * arguments:
 *  @state: state to check for validity.
 *  @len: number of nodes in the state.
 *
 * returns:
 *  integer indicating state validity (1) or non-validity (0).
 */
static inline int state_valid (enum_thread_node_t *state,
                               unsigned int len) {
  /* if any one index has passed the end value, return false. */
  for (unsigned int i = 0; i < len; i++) {
    if (state[i].idx > state[i].end)
      return 0;
  }

  /* return true. */
  return 1;
}

/* state_increment(): increment the index of a state.
 *
 * arguments:
 *  @state: state to increment.
 *  @len: number of nodes in the state.
 *  @lev: level at which to begin incrementation.
 *
 * returns:
 *  index of the most upstream node that was modified by the increment.
 */
static inline unsigned int state_increment (enum_thread_node_t *state,
                                            const unsigned int len,
                                            const unsigned int lev) {
  /* declare required variables:
   *  @imod: index of the highest modified parent node.
   */
  unsigned int imod = lev;

  /* loop over the skipped state indices (when lev < len - 1). */
  for (unsigned int i = lev + 1; i < len; i++)
    state[i].idx = 0;

  /* loop over the state indices. */
  for (unsigned int i = lev; i < len; i--) {
    /* increment the current index. */
    state[i].idx++;
    imod = i;

    /* check for overflow. */
    if (i && state[i].idx >= state[i].nb)
      state[i].idx = 0;
    else
      break;
  }

  /* return the highest modified node index. */
  return imod;
}

/* state_widths(): compute the number of leaf nodes in the implicit
 * tree that are to the left and right of the tree path defined by
 * the current state.
 *
 * arguments:
 *  @state: state to compute widths.
 *  @len: number of nodes in the state.
 *  @wl: base-10 log of the number of leaves to the left.
 *  @wr: base-10 log of the number of leaves to the right.
 */
static inline void state_widths (enum_thread_node_t *state,
                                 const unsigned int len,
                                 double *wl, double *wr) {
  /* initialize the output values. */
  long double L = 0.0, R = 0.0;

  /* loop over each level of the tree. */
  for (unsigned int j = 0; j < len; j++) {
    /* compute the size of each subtree under this level. */
    long double nj = 1.0;
    for (unsigned int k = j + 1; k < len; k++)
      nj *= (long double) state[k].nb;

    /* sum in the number of leaves to the left and right
     * of the state at this level.
     */
    L += nj * (long double) state[j].idx;
    R += nj * (long double) (state[j].nb - state[j].idx - 1);
  }

  /* return the base-10 logarithm of the leaf counts. */
  *wl = (double) log10l(L);
  *wr = (double) log10l(R);
}

/* state_rmsd(): determine the rmsd-step of a given state.
 *
 * this function utilizes the quaternion characteristic polynomial method,
 * referenced here:
 *
 *  D. L. Theobald, "Rapid calculation of RMSD using a quaternion-based
 *  characteristic polynomial." Acta Crysta. A, 2005, 61(4): 478-480.
 *
 * arguments:
 *  @state: state from which to compute the rmsd.
 *  @len: number of nodes in the state.
 *
 * returns:
 *  root-mean-square deviation (rmsd) from the previous solution to
 *  the current solution.
 */
static inline double state_rmsd (enum_thread_node_t *state,
                                 const unsigned int len) {
  /* declare required variables:
   */
  double Sxx, Sxy, Sxz, Syx, Syy, Syz, Szx, Szy, Szz, G1, G2;
  double x1, x2, y1, y2, z1, z2;
  unsigned int i;

  /* initialize the inner product matrix. */
  Sxx = Sxy = Sxz = 0.0;
  Syx = Syy = Syz = 0.0;
  Szx = Szy = Szz = 0.0;

  /* loop to compute the inner product matrix. */
  for (i = 0, G1 = 0.0, G2 = 0.0; i < len; i++) {
    /* get the first coordinate. */
    x1 = state[i].prev.x;
    y1 = state[i].prev.y;
    z1 = state[i].prev.z;

    /* get the second coordinate. */
    x2 = state[i].pos.x;
    y2 = state[i].pos.y;
    z2 = state[i].pos.z;

    /* compute the inner product terms. */
    G1 += x1 * x1 + y1 * y1 + z1 * z1;
    G2 += x2 * x2 + y2 * y2 + z2 * z2;

    /* update the first-row matrix elements. */
    Sxx += x1 * x2;
    Sxy += x1 * y2;
    Sxz += x1 * z2;

    /* update the second-row matrix elements. */
    Syx += y1 * x2;
    Syy += y1 * y2;
    Syz += y1 * z2;

    /* update the third-row matrix elements. */
    Szx += z1 * x2;
    Szy += z1 * y2;
    Szz += z1 * z2;
  }

  /* compute the final inner product. */
  const double E0 = (G1 + G2) * 0.5;
  double E = E0;

  /* compute squared matrix elements (diagonal). */
  const double Sxx2 = Sxx * Sxx;
  const double Syy2 = Syy * Syy;
  const double Szz2 = Szz * Szz;

  /* compute squared matrix elements (off-diagonals). */
  const double Sxy2 = Sxy * Sxy;
  const double Syz2 = Syz * Syz;
  const double Sxz2 = Sxz * Sxz;

  /* compute squared matrix elements (off-diagonals). */
  const double Syx2 = Syx * Syx;
  const double Szy2 = Szy * Szy;
  const double Szx2 = Szx * Szx;

  /* compute some temporaries... */
  const double SyzSzymSyySzz2 = 2.0 * (Syz * Szy - Syy * Szz);
  const double Sxx2Syy2Szz2Syz2Szy2 = Syy2 + Szz2 - Sxx2 + Syz2 + Szy2;

  /* ... and another temporary. */
  const double C2 = -2.0 * (Sxx2 + Syy2 + Szz2 + Sxy2 + Syx2 +
                            Sxz2 + Szx2 + Syz2 + Szy2);

  /* ... and yet another temporary. */
  const double C1 = 8.0 * (Sxx*Syz*Szy + Syy*Szx*Sxz + Szz*Sxy*Syx -
                           Sxx*Syy*Szz - Syz*Szx*Sxy - Szy*Syx*Sxz);

  /* ... and even more temporaries. */
  const double SxzpSzx = Sxz + Szx;
  const double SyzpSzy = Syz + Szy;
  const double SxypSyx = Sxy + Syx;
  const double SyzmSzy = Syz - Szy;
  const double SxzmSzx = Sxz - Szx;
  const double SxymSyx = Sxy - Syx;
  const double SxxpSyy = Sxx + Syy;
  const double SxxmSyy = Sxx - Syy;
  const double Sxy2Sxz2Syx2Szx2 = Sxy2 + Sxz2 - Syx2 - Szx2;

  /* ... and one final temporary. */
  const double C0 = Sxy2Sxz2Syx2Szx2 * Sxy2Sxz2Syx2Szx2 +
    (Sxx2Syy2Szz2Syz2Szy2 + SyzSzymSyySzz2) *
    (Sxx2Syy2Szz2Syz2Szy2 - SyzSzymSyySzz2) +
    (-SxzpSzx*SyzmSzy + SxymSyx*(SxxmSyy-Szz)) *
    (-SxzmSzx*SyzpSzy + SxymSyx*(SxxmSyy+Szz)) +
    (-SxzpSzx*SyzpSzy - SxypSyx*(SxxpSyy-Szz)) *
    (-SxzmSzx*SyzmSzy - SxypSyx*(SxxpSyy+Szz)) +
    (+SxypSyx*SyzpSzy + SxzpSzx*(SxxmSyy+Szz)) *
    (-SxymSyx*SyzmSzy + SxzpSzx*(SxxpSyy+Szz)) +
    (+SxypSyx*SyzmSzy + SxzmSzx*(SxxmSyy-Szz)) *
    (-SxymSyx*SyzpSzy + SxzmSzx*(SxxpSyy-Szz));

  /* newton-rhapson. */
  for (i = 0; i < 50; i++) {
    double Eprev = E;
    x2 = E * E;
    const double b = (x2 + C2) * E;
    const double a = b + C1;
    const double delta = (a * E + C0) / (2.0 * x2 * E + b + a);
    E -= delta;

    /* break on convergence. */
    if (fabs(E - Eprev) < fabs(1.0e-11 * E))
      break;
  }

  /* compute and return the rmsd. */
  return 2.0 * (E0 - E);
}

/* enum_threads_init(): initialize a set of threads prior to use.
 *
 * arguments:
 *  @E: pointer to the enumerator structure to modify.
 *
 * returns:
 *  integer indicating whether (1) or not (0) the operation succeeded.
 */
int enum_threads_init (enum_t *E) {
  /* compute the branch counts at each level. */
  for (unsigned int i = 0; i < E->G->n_order; i++) {
    /* check the special case of single or duplicate nodes. */
    if (i < 3 || E->G->orig[i]) {
      E->threads[0].state[i].nb = 1;
      continue;
    }

    /* determine whether the current vertex has an exact edge back to
     * a friend vertex, which will reduce the set of possible solutions
     * to two points, one point, or null.
     */
    unsigned int has_exact_edge = 0;
    for (unsigned int k = 0; k < E->G->n_friends[i]; k++) {
      value_t dik = graph_get_edge(E->G, E->G->order[i],
                                   E->G->friends[i][k]);

      if (value_is_scalar(dik)) {
        has_exact_edge = 1;
        break;
      }
    }

    /* set the branch count based on whether the current vertex has
     * an exact edge to a friend vertex.
     */
    if (has_exact_edge) {
      /* scalar edges produce two branches. */
      E->threads[0].state[i].nb = 2;
    }
    else {
      /* interval edges produce multiple branches. */
      const unsigned int i0 = E->G->order[i];
      const unsigned int i3 = E->G->order[i - 3];
      value_t d03 = graph_get_edge(E->G, i0, i3);
      unsigned int nb = (d03.u - d03.l) / E->eps;
      if (nb > E->nbmax)
        nb = E->nbmax;
      else if (nb == 0)
        nb = 1;

      /* set the branch count. */
      E->threads[0].state[i].nb = 2 * nb;
    }
  }

  /* copy the branch counts from the first thread
   * into every remaining thread.
   */
  for (unsigned int t = 1; t < E->nthreads; t++)
    for (unsigned int i = 0; i < E->G->n_order; i++)
      E->threads[t].state[i].nb = E->threads[0].state[i].nb;

  /* check if more than one thread will be executed. */
  if (E->nthreads > 1) {
    /* determine the division of labor. */
    state_div(E->threads[0].state, E->G->n_order, E->nthreads);

    /* set up the threads based on the computed division. */
    for (unsigned int t = 1; t < E->nthreads; t++) {
      /* compute the state start. */
      for (unsigned int i = 0; i < E->G->n_order; i++)
        E->threads[t].state[i].idx = E->threads[t - 1].state[i].end;

      /* compute the state end. */
      state_add(E->threads[t].state,
                E->threads[0].state,
                E->G->n_order);
    }

    /* set the end state. */
    for (unsigned int t = 1; t < E->nthreads; t++)
      for (unsigned int i = 0; i < E->G->n_order; i++)
        E->threads[t - 1].state[i].end = E->threads[t].state[i].idx;
  }

  /* set the end state to the problem size. */
  for (unsigned int i = 0; i < E->G->n_order; i++) {
    E->threads[E->nthreads - 1].state[i].end =
    E->threads[E->nthreads - 1].state[i].nb - 1;
  }

  /* set the start state to the initial index. */
  for (unsigned int t = 0; t < E->nthreads; t++)
    for (unsigned int i = 0; i < E->G->n_order; i++)
      E->threads[t].state[i].start = E->threads[t].state[i].idx;

  /* initialize the energies. */
  for (unsigned int t = 0; t < E->nthreads; t++)
    for (unsigned int i = 0; i < E->G->n_order; i++)
      E->threads[t].state[i].energy = 0.0;

  /* initialize the dihedral lists. */
  for (unsigned int t = 0; t < E->nthreads; t++) {
    for (unsigned int i = 0; i < E->G->n_order; i++) {
      E->threads[t].state[i].omega =
        malloc(E->threads[t].state[i].nb * sizeof(double));
      if (!E->threads[t].state[i].omega)
        throw("unable to allocate dihedrals array, level %u", i);
    }
  }

  /* compute the tree size. */
  for (unsigned int i = 0; i < E->G->n_order; i++)
    E->logW += log10((double) E->threads[0].state[i].nb);

  /* output an informational message about the tree size. */
  info("dense tree size: 10^%.3lf leaves", E->logW);

  /* return success. */
  return 1;
}

/* enum_thread_feasible(): determine whether an enumeration tree node
 * is feasible or not.
 *
 * arguments:
 *  @th: pointer to the thread to feasibility-check.
 *
 * returns:
 *  integer indicating whether (1) or not (0) the atom is feasible.
 */
static inline int enum_thread_feasible (enum_thread_t *th) {
  /* get local references to the pruning data. */
  enum_t *E = th->E;
  void **data = E->prune_data[th->level];
  enum_prune_test_fn *func = E->prune[th->level];
  const unsigned int n = E->prune_sz[th->level];

  /* loop over the array of pruning function pointers. */
  for (unsigned int i = 0; i < n; i++) {
    /* return infeasible if any function returns a prune. */
    if ((func[i])(E, th, data[i]))
      return 0;
  }

  /* return feasible. */
  return 1;
}

/* enum_thread_timer(): timer thread function for enumeration timing
 * information.
 *
 * arguments:
 *  @pdata: pointer to the enumerator data structure.
 */
void *enum_thread_timer (void *pdata) {
  /* get a reference to the current enumerator. */
  enum_t *E = (enum_t*) pdata;

  /* get references to commonly used sizes. */
  const unsigned int nt = E->nthreads;
  const unsigned int len = E->G->n_order;

  /* initialize the timer. */
  unsigned int last_nsol = 0;
  double t = 0.0, dt = 1.0;

  /* loop until the master thread cancels us. */
  while (1) {
    /* sleep for the required duration. */
    sleep((unsigned int) (dt * 60.0));
    t += dt;

    /* start outputting timing information. */
    info("elapsed time: %.0lf min.", t);

    /* compute the log-width of each thread. */
    for (unsigned int tid = 0; tid < nt; tid++) {
      /* compute the thread width. */
      double L1 = 0.0;
      double L2 = 0.0;
      state_widths(E->threads[tid].state, len, &L1, &L2);

      /* write the current thread information. */
      fprintf(stderr, "   #%-3u: ~10^%.3lf min. remaining\n",
              tid + 1, L2 - L1 + log(t));
    }

    /* if the number of solutions has not increased since our
     * last timing report, increase the delay time, up to a maximum
     * of two hours between reports. when a solution is logged, set
     * the delay time back to the minimum.
     */
    if (E->nsol <= last_nsol)
      dt = (dt * 5.0 > 120.0 ? 120.0 : dt * 5.0);
    else
      dt = 1.0;

    /* store the new solution count. */
    last_nsol = E->nsol;
  }

  /* end thread execution. */
  return NULL;
}

/* enum_thread_execute(): core thread function for enumerator threads.
 *
 * arguments:
 *  @pdata: pointer to the current enumerator thread data structure.
 */
void *enum_thread_execute (void *pdata) {
  /* get references to the current thread data. */
  enum_thread_t *thread = (enum_thread_t*) pdata;
  enum_thread_node_t *state = thread->state;
  graph_t *G = thread->E->G;
  enum_t *E = thread->E;

  /* get a reference to the length of the order. */
  const unsigned int len = G->n_order;
  const unsigned int *dup = G->orig;
  unsigned int lev = thread->level;

  /* get the distances required to embed the first three atoms. */
  double d01 = graph_get_edge_exact(G, G->order[0], G->order[1]);
  double d02 = graph_get_edge_exact(G, G->order[0], G->order[2]);
  double d12 = graph_get_edge_exact(G, G->order[1], G->order[2]);

  /* define distances to each newly embedded atom. */
  double d13, d23;
  value_t val03;

  /* compute the cosine and sine of the angle formed by the atoms. */
  double ct = distances_to_angle(d01, d02, d12);
  double st = sqrt(1.0 - ct * ct);

  /* define dihedral angular quantities for embedding atoms. */
  double cw, sw;

  /* define vector quantities and extra scalars for embedding atoms. */
  vector_t x0, x1, x2, x3, r01, r02, r12, rv, p1, p2, p3;
  double fp, fv, fd;

  /* initialize the first three atom positions. */
  vector_set(&state[0].pos, 0.0, 0.0, 0.0);
  vector_set(&state[1].pos, -d01, 0.0, 0.0);
  vector_set(&state[2].pos, d12 * ct - d01, d12 * st, 0.0);

  /* loop over the set of states apportioned to the thread. */
  while (state_valid(state, len)) {
    /* check if we should terminate enumeration. */
    if (E->term)
      return NULL;

    /* check if we've computed enough solutions. */
    if (E->nmax && E->nsol >= E->nmax)
      return NULL;

    /* embed all modified atoms in the state. */
    while (lev < len) {
      /* do not attempt to embed the first three atoms. */
      if (lev < 3)
        lev = 3;

      /* check for duplicate atoms. */
      if (dup[lev]) {
        /* store the previously computed position, and move on. */
        state[lev].pos = state[lev - dup[lev]].pos;
        state[lev].energy = state[lev - dup[lev]].energy;
        lev++; continue;
      }

      /* pull some embedded atom positions into local variables. */
      x0 = state[lev - 3].pos;
      x1 = state[lev - 2].pos;
      x2 = state[lev - 1].pos;

      /* r01 = x1 - x0 == x_{i-2} - x_{i-3} */
      r01.x = x1.x - x0.x;
      r01.y = x1.y - x0.y;
      r01.z = x1.z - x0.z;

      /* r02 = x2 - x0 == x_{i-1} - x_{i-3} */
      r02.x = x2.x - x0.x;
      r02.y = x2.y - x0.y;
      r02.z = x2.z - x0.z;

      /* r12 = x2 - x1 == x_{i-1} - x_{i-2} */
      r12.x = x2.x - x1.x;
      r12.y = x2.y - x1.y;
      r12.z = x2.z - x1.z;

      /* rv = cross(r12, r01) */
      rv.x = r12.y * r01.z - r12.z * r01.y;
      rv.y = r12.z * r01.x - r12.x * r01.z;
      rv.z = r12.x * r01.y - r12.y * r01.x;

      /* fd = dot(r12, r01) */
      fd = r12.x * r01.x + r12.y * r01.y + r12.z * r01.z;

      /* compute distances between the previously embedded atoms. */
      d01 = sqrt(r01.x * r01.x + r01.y * r01.y + r01.z * r01.z);
      d02 = sqrt(r02.x * r02.x + r02.y * r02.y + r02.z * r02.z);
      d12 = sqrt(r12.x * r12.x + r12.y * r12.y + r12.z * r12.z);

      /* obtain distances to the atom to be embedded. */
      val03 = graph_get_edge(G, G->order[lev - 3], G->order[lev]);
      d13 = graph_get_edge_exact(G, G->order[lev - 2], G->order[lev]);
      d23 = graph_get_edge_exact(G, G->order[lev - 1], G->order[lev]);

      /* compute the cosine and sine of theta. */
      ct = distances_to_angle(d12, d13, d23);
      st = sqrt(1.0 - ct * ct);

/*FIXME*/
printf("\n[***    LEVEL %3u    ***]\n",lev);
printf("[idx,   nb,  start,  end]\n");
printf(" %3u   %3u     %3u   %3u\n",
state[lev].idx,state[lev].nb,state[lev].start,state[lev].end);
/*FIXME*/
      /* upon entering a level at the starting index, compute
       * the set of dihedral points to use for embedding.
       */
      if (state[lev].idx == state[lev].start) {
        /* check if the (i,i-3) dihedral is directly available. */
        if (value_is_dihedral(val03)) {
          /* dihedral case: bound the initial dihedrals. */
          val03 = value_bound(value_scal(*val03.src, M_PI / 180.0),
                              value_interval(-M_PI, M_PI));
        }
        else {
          /* distance case: do not bound the initial dihedrals. */
          val03 = value_interval(-M_PI, M_PI);
        }

        /* compute the discretized dihedral values. */
        enum_reduce(thread, lev, val03.l, val03.u);
      }
/*FIXME*/
printf("omega (%u) = ( ",state[lev].n_omega);
for (unsigned int kk = 0; kk < state[lev].n_omega; kk++)
  printf("%lf%s", state[lev].omega[kk],
         kk + 1 == state[lev].n_omega ? "" : ", ");
printf(" )\n");
/*FIXME*/

      /* prune if we have exhausted our set of feasible points
       * at the current level.
       */
      if (state[lev].idx >= state[lev].n_omega) {
        lev = state_increment(state, len, lev - 1);
        goto infeasible;
      }

      /* compute the cosine and sine of omega. */
      /* FIXME: THIS IS ABSOLUTELY NOT MULTITHREAD COMPATIBLE!!!
       *  -- see notes in src/opts.c on future remedies.
       */
      cw = cos(state[lev].omega[state[lev].idx]);
      sw = sin(state[lev].omega[state[lev].idx]);
/*FIXME*/printf("cos(omega) = %lf\n",cw);
/*FIXME*/printf("sin(omega) = %lf\n",sw);

      /* compute the scale factor for all p-vectors. */
      fv = st / sqrt(rv.x * rv.x + rv.y * rv.y + rv.z * rv.z);
      fp = -d23 / d12;

      /* compute the first anchor position. */
      p1.x = fp * ((ct + 1.0 / fp) * x2.x - ct * x1.x);
      p1.y = fp * ((ct + 1.0 / fp) * x2.y - ct * x1.y);
      p1.z = fp * ((ct + 1.0 / fp) * x2.z - ct * x1.z);

      /* compute the second anchor position. */
      fp *= fv;
      p2.x = fp * (d12 * d12 * r01.x - fd * r12.x);
      p2.y = fp * (d12 * d12 * r01.y - fd * r12.y);
      p2.z = fp * (d12 * d12 * r01.z - fd * r12.z);

      /* compute the third anchor position. */
      p3.x = fp * d12 * rv.x;
      p3.y = fp * d12 * rv.y;
      p3.z = fp * d12 * rv.z;

      /* compute and store the newly embedded atom position. */
      x3.x = p1.x + cw * p2.x + sw * p3.x;
      x3.y = p1.y + cw * p2.y + sw * p3.y;
      x3.z = p1.z + cw * p2.z + sw * p3.z;
      state[lev].pos = x3;

      /* check feasibility of the newly embedded atom. */
      thread->level = lev;
      if (!enum_thread_feasible(thread)) {
        /* infeasible:
         *  1. skip all sub-trees of the infeasible atom/node.
         *  2. move back into the loop without incrementing.
         */
        lev = state_increment(state, len, lev);
        goto infeasible;
      }

      /* check if the atom is feasible and terminal. */
      if (lev == len - 1) {
        /* compute the center of the structure. */
        x0.x = x0.y = x0.z = 0.0;
        for (unsigned int i = 0; i < len; i++) {
          x0.x += state[i].pos.x;
          x0.y += state[i].pos.y;
          x0.z += state[i].pos.z;
        }

        /* scale the computed mean value. */
        fp = 1.0 / ((double) len);
        x0.x *= fp;
        x0.y *= fp;
        x0.z *= fp;

        /* center the coordinates of the current candidate solution. */
        for (unsigned int i = 0; i < len; i++) {
          state[i].pos.x -= x0.x;
          state[i].pos.y -= x0.y;
          state[i].pos.z -= x0.z;
        }

        /* break if:
         *  1. the rmsd-step of the candidate solution is too low.
         *  2. the energy of the candidate solution is too high.
         */
        if (state_rmsd(state, len) < E->rmsd_tol ||
            state[len - 1].energy > E->energy_tol) {
          E->nrej++;
          break;
        }

        /* store the new solution into the "previous" slot. */
        for (unsigned int i = 0; i < len; i++)
          state[i].prev = state[i].pos;

#ifdef __IBP_HAVE_PTHREAD
        /* obtain a lock on the write mutex. */
        pthread_mutex_lock(&E->write_mutex);
#endif

        /* increment the solution count. */
        E->nsol++;
        E->energy_tol = state[len - 1].energy;

        /* write some output. */
        info("solution %u found, U = %.32le",
             E->nsol, E->energy_tol);

        /* write the solution. */
        if (E->write_data && !E->write_data(E, thread)) {
          /* raise an exception and end thread execution. */
          raise("failed to write solution %u", E->nsol);
          return NULL;
        }

#ifdef __IBP_HAVE_PTHREAD
        /* unlock the write mutex. */
        pthread_mutex_unlock(&E->write_mutex);
#endif
      }

      /* move down a level. */
      lev++;
    }

    /* increment the state. */
    lev = state_increment(state, len, len - 1);

/* causes the thread to re-enter the level loop without an increment. */
infeasible:;
  }

  /* end thread execution. */
  return NULL;
}

