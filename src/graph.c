
/* include the graph header. */
#include "graph.h"

/* graph_new(): allocate a new empty graph data structure.
 *
 * this function creates a pointer to a fully initialized graph_t data type,
 * which may then be filled with information using other graph_*() functions
 * and/or freed using graph_free().
 *
 * arguments:
 *  @n_vertices: number of vertices of the graph.
 *
 * returns:
 *  pointer to an allocated and initialized graph_t structure. the pointer
 *  must be freed after use by graph_free().
 */
graph_t *graph_new (unsigned int n_vertices) {
  /* declare required variables:
   *  @G: pointer to a newly allocated graph.
   */
  graph_t *G;

  /* allocate the graph structure pointer. */
  G = (graph_t*) malloc(sizeof(graph_t));

  /* check if the pointer was allocated successfully. */
  if (!G) {
    /* raise an exception and return null. */
    raise("unable to allocate graph structure pointer");
    return NULL;
  }

  /* initialize the vertex count. */
  G->nv = n_vertices;

  /* allocate the edge array. */
  G->E = (value_t*) malloc(G->nv * G->nv * sizeof(value_t));

  /* check if the array was allocated successfully. */
  if (!G->E) {
    /* free the structure pointer. */
    free(G);

    /* raise an exception and return null. */
    raise("unable to allocate graph edge array");
    return NULL;
  }

  /* allocate the vertex reverse-lookup array. */
  G->ordrev = (unsigned int*) malloc(G->nv * sizeof(unsigned int));

  /* check if the array was allocated successfully. */
  if (!G->ordrev) {
    /* free the structure pointer. */
    free(G);

    /* raise an exception and return null. */
    raise("unable to allocate graph reverse lookup array");
    return NULL;
  }

  /* initialize the edge array. */
  for (unsigned int i = 0; i < G->nv; i++)
    for (unsigned int j = 0; j < G->nv; j++)
      G->E[i + G->nv * j] = value_undefined();

  /* initialize the reverse-lookup array. */
  for (unsigned int i = 0; i < G->nv; i++)
    G->ordrev[i] = UINT_MAX;

  /* initialize the re-order array. */
  G->rmsd = NULL;
  G->order = G->orig = NULL;
  G->n_order = G->n_orig = 0;

  /* return the newly allocated graph. */
  return G;
}

/* graph_free(): free all allocated memory associated with an iDMDGP graph.
 *
 * arguments:
 *  @G: pointer to the graph structure to free.
 */
void graph_free (graph_t *G) {
  /* return if the structure pointer is null. */
  if (!G) return;

  /* free the edge array. */
  free(G->E);
  G->E = NULL;

  /* free the reverse lookup array. */
  free(G->ordrev);
  G->ordrev = NULL;

  /* reset the vertex count. */
  G->nv = 0;

  /* check if the re-order arrays are allocated. */
  if (G->n_order) {
    /* yes. free the arrays and reset the array size. */
    free(G->order);
    free(G->orig);
    free(G->rmsd);
    G->n_order = 0;
    G->n_orig = 0;
  }

  /* finally, free the structure pointer. */
  free(G);
}

/* graph_has_edge(): test whether two vertices of a graph are adjacent.
 *
 * arguments:
 *  @G: pointer to the graph structure to access.
 *  @va: index of the first graph vertex to test.
 *  @vb: index of the second graph vertex to test.
 *
 * returns:
 *  value type of the graph edge, which will evaluate to zero for
 *  undefined values (i.e. missing edges).
 */
value_type_t graph_has_edge (graph_t *G, unsigned int va, unsigned int vb) {
  /* return true if the graph edge is either scalar or interval. */
  return G->E[va + G->nv * vb].type;
}

/* graph_get_edge(): get the generalized value of an edge between
 * two indexed vertices of a graph, if such an edge exists.
 *
 * arguments:
 *  @G: pointer to the graph structure to access.
 *  @va, @vb: vertices of the query edge.
 *
 * returns:
 *  value of the edge, or an undefined value if no such edge was found.
 *
 * note:
 *  this function performs no bounds checking. the calling function must
 *  ensure that the specified indices @va and @vb are in bounds.
 */
value_t graph_get_edge (graph_t *G, unsigned int va, unsigned int vb) {
  /* directly return the edge value. */
  return (va < G->nv && vb < G->nv
            ? G->E[va + G->nv * vb]
            : value_undefined());
}

/* graph_get_edge_exact(): get the value of an exact edge between two
 * vertices of a graph, if such an edge exists.
 *
 * arguments:
 *  @G: pointer to the graph structure to access.
 *  @va, @vb: vertices of the query edge.
 *
 * returns:
 *  value of the exact edge, or NAN if no such edge was found.
 *
 * note:
 *  this function performs no bounds checking. the calling function must
 *  ensure that the specified indices @va and @vb are in bounds.
 */
double graph_get_edge_exact (graph_t *G, unsigned int va, unsigned int vb) {
  /* directly return the edge lower bound. */
  return G->E[va + G->nv * vb].l;
}

/* graph_set_edge(): set the generalized value of an edge between two
 * vertices of a graph.
 *
 * arguments:
 *  @G: pointer to the graph structure to access.
 *  @va, @vb: vertices of the edge to set.
 *  @w: edge weight to set.
 */
void graph_set_edge (graph_t *G, unsigned int va, unsigned int vb,
                     value_t w) {
  /* return if the edge indices are out of bounds or equal. */
  if (va == vb || va >= G->nv || vb >= G->nv)
    return;

  /* directly assign the edge value. */
  G->E[va + G->nv * vb] = G->E[vb + G->nv * va] = w;
}

/* graph_refine_edge(): high-level function for refining the weight of
 * a graph edge.
 *
 * if no edge exists between the two vertices, then an edge is added to
 * the appropriate array (depending on whether the weight is a scalar
 * or an interval).
 *
 * however, if an edge exists, then the weight of the existing edge is
 * refined when possible:
 *
 * exists:        refinement:        result:
 *  scalar         scalar             error
 *  interval       scalar             replacement
 *  interval       interval           refinement
 *
 * - replacements remove the existing inexact edge and add a new exact edge.
 * - refinements update the existing inexact edge by the new inexact edge,
 *   when the lower and/or upper bound of the new edge is an improvement
 *   over the previous bound.
 *
 * arguments:
 *  @G: pointer to the graph structure to modify.
 *  @va: index of the first graph vertex.
 *  @vb: index of the second graph vertex.
 *  @w: refinement weight of the graph edge.
 *
 * returns:
 *  integer indicating whether (1) or not (0) the operation succeeded.
 */
int graph_refine_edge (graph_t *G, unsigned int va, unsigned int vb,
                       value_t w) {
  /* declare required variables:
   *  @wcur: current edge weight.
   */
  value_t wcur;

  /* check that the graph pointer is valid. */
  if (!G)
    throw("graph structure pointer is invalid");

  /* check if an exact edge already exists. */
  if (graph_has_edge(G, va, vb) == VALUE_TYPE_SCALAR)
    return 1;

  /* determine the type of the refinement weight. */
  if (value_is_scalar(w)) {
    /* add the new exact edge to the graph. */
    graph_set_edge(G, va, vb, w);
  }
  else if (value_is_interval(w)) {
    /* check if an interval edge already exists. */
    if (graph_has_edge(G, va, vb) == VALUE_TYPE_INTERVAL) {
      /* an edge exists. get its value. */
      wcur = graph_get_edge(G, va, vb);

      /* refine its lower bound, if possible. */
      if (w.l > wcur.l)
        wcur.l = w.l;

      /* and refine its upper bound, if possible. */
      if (w.u < wcur.u)
        wcur.u = w.u;

      /* validate the new edge bounds. */
      if (wcur.l > wcur.u)
        throw("invalid refined interval edge [%.3lf,%.3lf]",
              wcur.l, wcur.u);

      /* update the interval edge in the graph. */
      graph_set_edge(G, va, vb, wcur);
    }
    else {
      /* no edge exists. add a new interval edge. */
      graph_set_edge(G, va, vb, w);
    }
  }
  else {
    /* undefined weight. throw an exception. */
    throw("edge weight is undefined");
  }

  /* return success. */
  return 1;
}

/* graph_count_edges(): count the number of edges in a graph.
 *
 * arguments:
 *  @G: pointer to the graph structure to access.
 *  @ne: pointer to the output variable for the exact edge count.
 *  @ni: pointer to the output variable for the interval edge count.
 */
void graph_count_edges (graph_t *G, unsigned int *ne, unsigned int *ni) {
  /* initialize the output values. */
  unsigned int le = 0, li = 0;

  /* loop over the edge array. */
  for (unsigned int i = 0; i < G->nv; i++) {
    for (unsigned int j = i + 1; j < G->nv; j++) {
      /* get the type of the current edge. */
      const value_type_t et = graph_has_edge(G, i, j);

      /* update the edge counts. */
      switch (et) {
        case VALUE_TYPE_SCALAR:   le++; break;
        case VALUE_TYPE_INTERVAL: li++; break;
        default: break;
      }
    }
  }

  /* dirty trick to get rid of warnings. */
  *ne = le;
  *ni = li;
}

/* graph_extend_order(): append a new vertex into the order array
 * of a graph.
 *
 * arguments:
 *  @G: pointer to the graph structure to modify.
 *  @v: vertex index to add to the order.
 *
 * returns:
 *  integer indicating whether (1) or not (0) the operation was successful.
 */
int graph_extend_order (graph_t *G, unsigned int v) {
  /* declare required variables:
   *  @i: order array index.
   *  @j: array loop counter.
   */
  unsigned int i, j;

  /* check that the graph pointer is valid. */
  if (!G || G->nv == 0)
    throw("graph structure pointer is invalid");

  /* check that the vertex index is in bounds. */
  if (v >= G->nv)
    throw("vertex index %u out of bounds [0,%u]", v, G->nv - 1);

  /* increment the array length. */
  i = G->n_order;
  G->n_order++;

  /* reallocate the order array. */
  G->order = (unsigned int*)
    realloc(G->order, G->n_order * sizeof(unsigned int));

  /* reallocate the originality array. */
  G->orig = (unsigned int*)
    realloc(G->orig, G->n_order * sizeof(unsigned int));

  /* reallocate the rmsd array. */
  G->rmsd = (double*)
    realloc(G->rmsd, G->n_order * sizeof(double));

  /* check if reallocation failed. */
  if (!G->order || !G->orig || !G->rmsd)
    throw("unable to reallocate re-order arrays");

  /* store the new array element. */
  G->order[i] = v;
  G->orig[i] = 0;
  G->rmsd[i] = 0.0;

  /* loop over the previous elements to determine originality. */
  for (j = 0; j < i; j++) {
    /* check if the current node is a prior visit. */
    if (G->order[j] == v) {
      /* yes, this is the original visit. */
      G->orig[i] = i - j;
      break;
    }
  }

  /* if this is an original node... */
  if (G->orig[i] == 0) {
    /* store the reverse-lookup index and increment the
     * number of original vertices.
     */
    G->ordrev[v] = i;
    G->n_orig++;
  }

  /* return success. */
  return 1;
}

