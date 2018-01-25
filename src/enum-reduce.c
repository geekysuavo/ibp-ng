
/* include the enumerator headers. */
#include "enum.h"
#include "enum-thread.h"

/* _det3x3(): macro function for computing the determinant of
 * a three-by-three matrix, stored linearly in row-major order.
 */
#define _det3x3(A) \
  (A[0] * A[4] * A[8] \
 - A[0] * A[5] * A[7] \
 - A[1] * A[3] * A[8] \
 + A[1] * A[5] * A[6] \
 + A[2] * A[3] * A[7] \
 - A[2] * A[4] * A[6])

/* _solve_det3x3_0(): macro function to solve for the *first* unknown
 * of the system A*x=b, where A is an invertible 3x3 matrix.
 */
#define _solve_det3x3_0(A, b) \
  (b[0] * A[4] * A[8] \
 - b[0] * A[5] * A[7] \
 - A[1] * b[1] * A[8] \
 + A[1] * A[5] * b[2] \
 + A[2] * b[1] * A[7] \
 - A[2] * A[4] * b[2])

/* _solve_det3x3_1(): macro function to solve for the *second* unknown
 * of the system A*x=b, where A is an invertible 3x3 matrix.
 */
#define _solve_det3x3_1(A, b) \
  (A[0] * b[1] * A[8] \
 - A[0] * A[5] * b[2] \
 - b[0] * A[3] * A[8] \
 + b[0] * A[5] * A[6] \
 + A[2] * A[3] * b[2] \
 - A[2] * b[1] * A[6])

/* _solve_det3x3_2(): macro function to solve for the *third* unknown
 * of the system A*x=b, where A is an invertible 3x3 matrix.
 */
#define _solve_det3x3_2(A, b) \
  (A[0] * A[4] * b[2] \
 - A[0] * b[1] * A[7] \
 - A[1] * A[3] * b[2] \
 + A[1] * b[1] * A[6] \
 + b[0] * A[3] * A[7] \
 - b[0] * A[4] * A[6])

/* add_circular_interval(): union an interval into an interval set
 * while respecting the fact that all intervals lie on the circle,
 * i.e. \theta \in [-pi, pi).
 *
 * arguments:
 *  @I: pointer to the interval set to modify.
 *  @a: lower bound of the new interval.
 *  @b: upper bound of the new interval.
 */
void add_circular_interval (intervals_t *I, double a, double b) {
  /* if the interval bounds are reversed, swap them. */
  if (a > b) {
    const double c = a;
    a = b;
    b = c;
  }

  /* compute the interval length. */
  const double len = b - a;

  /* check if the interval crosses the point pi = -pi. */
  if (len <= (M_PI + 1.0e-8)) {
    /* no. add the entire interval. */
    intervals_union(I, a, b);
  }
  else {
    /* yes. break the interval in two and add each separately. */
    intervals_union(I, -M_PI, a);
    intervals_union(I, b, M_PI);
  }
}

/* solve_3x3(): solve the linear system of equations: A * p = b
 *
 * where A is an invertible 3x3 matrix, stored linearly in row-major order,
 * b is a 3-element array, and p is the output vector to store the result
 * into.
 *
 * arguments:
 *  @A: array of nine (9) matrix elements.
 *  @b: array of three (3) vector elements.
 *  @p: vector structure pointer to store the result into.
 */
void solve_3x3 (const double *A, const double *b, vector_t *p) {
  /* compute the matrix determinant. */
  const double detA = _det3x3(A);

  /* compute each unknown vector element. */
  p->x = _solve_det3x3_0(A, b) / detA;
  p->y = _solve_det3x3_1(A, b) / detA;
  p->z = _solve_det3x3_2(A, b) / detA;
}

/* solve_tsi(): compute the intersection of three spheres as a linear
 * system with two solutions.
 *
 * arguments:
 *  @a, @b, @c: centers of the three spheres.
 *  @ra, @rb, @rc: radii of the three spheres.
 *  @n: unit normal vector to the plane (a,b,c).
 *  @p: projection of the solution point-pair onto the plane (a,b,c).
 *  @alpha: distance between @p and each point in the solution point-pair.
 */
void solve_tsi (vector_t *a, vector_t *b, vector_t *c,
                double ra,   double rb,   double rc,
                vector_t *n, vector_t *p, double *alpha) {
  /* declare a pair of temporary vectors. */
  vector_t u, v;

  /* u <- b - a */
  u.x = b->x - a->x;
  u.y = b->y - a->y;
  u.z = b->z - a->z;

  /* v <- c - a */
  v.x = c->x - a->x;
  v.y = c->y - a->y;
  v.z = c->z - a->z;

  /* n \propto cross(u, v). */
  vector_cross(&u, &v, n);
  vector_normalize(n);

  /* define the left-hand matrix. */
  const double A[] = {
    u.x,  u.y,  u.z,
    v.x,  v.y,  v.z,
    n->x, n->y, n->z
  };

  /* compute some intermediate values. */
  const double a2 = a->x * a->x + a->y * a->y + a->z * a->z;
  const double b2 = b->x * b->x + b->y * b->y + b->z * b->z;
  const double c2 = c->x * c->x + c->y * c->y + c->z * c->z;
  const double n_dot_a = n->x * a->x + n->y * a->y + n->z * a->z;

  /* define the right-hand vector. */
  const double y[] = { 
    (b2 - a2 - rb * rb + ra * ra) / 2.0,
    (c2 - a2 - rc * rc + ra * ra) / 2.0,
    n_dot_a
  };

  /* solve the linear system. */
  solve_3x3(A, y, p);

  /* u <- p - a */
  u.x = p->x - a->x;
  u.y = p->y - a->y;
  u.z = p->z - a->z;
  /* FIXME: issues with solving the systems in this way...
   * in some cases, |p-a|^2 > ra^2, which results in
   * alpha becoming nan.
   */

  /* compute the plane-to-solution length. */
  *alpha = sqrt((ra * ra) - (u.x * u.x + u.y * u.y + u.z * u.z));
}

/* solve_iomega_k(): compute the two dihedral interval arcs belonging
 * to the set of vertices (xk,x2,x1,x0), computed in the reference
 * frame of the vertices (x3,x2,x1,x0), for unknown x0.
 *
 * arguments:
 *  @x1: once-removed prior vertex in the order.
 *  @x2: twice-removed prior vertex in the order.
 *  @x3: thrice-removed prior vertex in the order.
 *  @xk: another prior vertex in the order (k <= i-3).
 *  @d01: distance between x[i] and x[i-1].
 *  @d02: distance between x[i] and x[i-2].
 *  @lk: lower bound distance between x[i] and x[k].
 *  @uk: upper bound distance between x[i] and x[k].
 *  @iomega_k: output interval set to hold the dihedral arcs.
 */
void solve_iomega_k (vector_t *x1, vector_t *x2, vector_t *x3, vector_t *xk,
                     double d01, double d02, double lk, double uk,
                     intervals_t *iomega_k) {
  /* declare required variables:
   *  @p: midpoint between the two solutions.
   *  @n: unit normal vector to plane (x0,x1,xk).
   *  @alpha: distance from midpoint to each solution.
   */
  vector_t p, n;
  double alpha;

  /* compute the arc points closest to x[k], and their omega values. */
  solve_tsi(x1, x2, xk, d01, d02, lk, &n, &p, &alpha);
  vector_axpy(&p, alpha, &n); /* p <- p + alpha n */
  const double a0 = vector_dihedral(x3, x2, x1, &p);
  vector_axpy(&p, -2.0 * alpha, &n); /* p <- p - 2 alpha n */
  const double b0 = vector_dihedral(x3, x2, x1, &p);

  /* compute the arc points furthest from x[k], and their omega values. */
  solve_tsi(x1, x2, xk, d01, d02, uk, &n, &p, &alpha);
  vector_axpy(&p, alpha, &n); /* p <- p + alpha n */
  const double a1 = vector_dihedral(x3, x2, x1, &p);
  vector_axpy(&p, -2.0 * alpha, &n); /* p <- p - 2 alpha n */
  const double b1 = vector_dihedral(x3, x2, x1, &p);

  /* store the results into the interval set. */
  iomega_k->size = 0;
  add_circular_interval(iomega_k, a0, a1);
  add_circular_interval(iomega_k, b0, b1);
}

/* enum_reduce(): at a given level within an enumerator thread, compute
 * a discretized set of omega dihedral angles using interval reduction,
 * described in the 'torsion Branch-and-Prune' paper.
 *
 * arguments:
 *  @th: enumerator thread structure pointer.
 *  @lev: level in the order at which to compute the omega values.
 *  @l0: initial lower bound for the reduced interval set.
 *  @u0: initial upper bound for the reduced interval set.
 *
 * returns:
 *  integer indicating whether (1) or not (0) dihedral interval reduction
 *  yielded a feasible (i.e. non-empty) set of intervals.
 */
int enum_reduce (enum_thread_t *th, unsigned int lev,
                 double l0, double u0) {
  /* get some required struct pointers:
   *  @E: master enumerator.
   *  @G: distance graph.
   *  @state: enumerator state.
   */
  enum_t *E = th->E;
  graph_t *G = E->G;
  enum_thread_node_t *state = th->state;

  /* get the positions of the three preceeding vertices in the order. */
  vector_t x3 = state[lev - 3].pos;
  vector_t x2 = state[lev - 2].pos;
  vector_t x1 = state[lev - 1].pos;

  /* declare variables for friends:
   *  @xk: position of the friend vertex.
   *  @d0k: distance to the friend vertex.
   */
  vector_t xk;
  value_t d0k;

  /* get some more required variables:
   *  @n_omega: number of discretization points.
   *  @isa, @isb, @isk: temporary interval sets.
   */
  const unsigned int n_omega = state[lev].nb;
  intervals_t *isa, *isb, *isk;
  isk = state[lev].isk;

  /* get some required uint arrays:
   *  @order: graph repetition order array.
   *  @ordrev: graph inverse re-order array.
   *  @friends: friend list for the current vertex.
   *  @n_friends: number of friends of the current vertex.
   */
  const unsigned int *order = G->order;
  const unsigned int *ordrev = G->ordrev;
  const unsigned int *friends = G->friends[lev];
  const unsigned int n_friends = G->n_friends[lev];

  /* ugh, and some more required variables:
   *  @v: graph vertex index at the current level in the order.
   *  @d01: distance to the once-removed vertex.
   *  @d02: distance to the twice-removed vertex.
   */
  const unsigned int v = order[lev];
  const double d01 = graph_get_edge(G, v, order[lev - 1]).l;
  const double d02 = graph_get_edge(G, v, order[lev - 2]).l;

  /* initialize the interval set to the entire circle. */
  state[lev].isa->size = 0;
  state[lev].isb->size = 0;
  intervals_union(state[lev].isa, l0, u0);

  /* loop over all friend vertices (adjacent non-repeat vertices
   * prior to the current vertex in the order), but do not use
   * the last two friends, because they are a part of the
   * 3-clique (they are x[i-1] and x[i-2]).
   */
  for (unsigned int k = 0; k < n_friends; k++) {
    /* alternate between the two interval sets for intersections. */
    isa = k % 2 ? state[lev].isb : state[lev].isa;
    isb = k % 2 ? state[lev].isa : state[lev].isb;

    /* get the vertex index and position of the current friend. */
    const unsigned int vk = friends[k];
    xk = state[ordrev[vk]].pos;

    /* get the graph edge connecting us to the current friend. */
    d0k = graph_get_edge(G, v, vk);

    /* solve for the two interval arcs related to the current friend. */
    solve_iomega_k(&x1, &x2, &x3, &xk, d01, d02,
    /* FIXME: always adding tolerances to d(i,k) can hurt us bad!
     * sometimes, d(i,k) is exact, so tolerances will cause pruning.
     * however, removing tolerances will require a re-implementation
     * of intervals_grid() to support zero-length intervals.
     */
                   d0k.l - E->ddf_tol / 2.0,
                   d0k.u + E->ddf_tol / 2.0,
                   isk);

    /* intersect these interval arcs with the current interval set. */
    intervals_intersect(isa, isk, isb);

    /* the interval set is empty! return false. */
    if (isb->size == 0)
      return 0;
  }

  /* finally, discretize the reduced interval set. */
  intervals_grid(isb, state[lev].omega, n_omega);

  /* feasible intervals are available, return true. */
  return 1;
}

