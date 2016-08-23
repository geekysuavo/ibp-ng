
/* ensure once-only inclusion. */
#ifndef __IBPNG_GRAPH_H__
#define __IBPNG_GRAPH_H__

/* include the traceback and value headers. */
#include "trace.h"
#include "value.h"

/* graph_t: structure for holding a single interval discretizable molecular
 * distance geometry problem (iDMDGP) instance in graph form.
 *
 * vertices:
 *  no metadata is stored per-vertex, in order to keep the graph structure
 *  as simple, abstract and lightweight as possible. the vertex count of
 *  a graph is set at allocation time.
 *
 * edges:
 *  graph edges are stored within a flat array of value_t structures,
 *  treated internally as a hollow symmetric matrix.
 */
typedef struct {
  /* @E: two-dimensional array of edges in the graph.
   * @nv: number of vertices in the graph.
   */
  value_t *E;
  unsigned int nv;

  /* @order: re-order array for graph traversal.
   * @orig: re-order originality offset array.
   * @rmsd: deviation contribution array.
   * @n_order: length of the re-order array.
   * @n_orig: number of original atoms in the order.
   */
  unsigned int *order, *orig, n_order, n_orig;
  double *rmsd;
}
graph_t;

/* function declarations (graph.c): */

graph_t *graph_new (unsigned int n_vertices);

void graph_free (graph_t *G);

value_type_t graph_has_edge (graph_t *G, unsigned int va, unsigned int vb);

value_t graph_get_edge (graph_t *G, unsigned int va, unsigned int vb);

double graph_get_edge_exact (graph_t *G, unsigned int va, unsigned int vb);

void graph_set_edge (graph_t *G, unsigned int va, unsigned int vb,
                     value_t w);

int graph_refine_edge (graph_t *G, unsigned int va, unsigned int vb,
                       value_t w);

#define graph_remove_edge(G,va,vb) \
  graph_set_edge(G, va, vb, value_undefined())

void graph_count_edges (graph_t *G, unsigned int *ne, unsigned int *ni);

int graph_extend_order (graph_t *G, unsigned int v);

#endif  /* !__IBPNG_GRAPH_H__ */

