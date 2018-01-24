
/* ensure once-only inclusion. */
#pragma once

/* include the integer limits header. */
#include <limits.h>

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
  /* core graph elements:
   *
   *  @E: two-dimensional array of edges in the graph.
   *  @nv: number of vertices in the graph.
   */
  value_t *E;
  unsigned int nv;

  /* graph properties related to the repetition order:
   *
   *  @order: re-order array for graph traversal. at a given level 'i'
   *          in the order, order[i] gives the vertex index.
   *  @ordrev: re-order reverse-lookup array. for a given vertex 'v'
   *           in the graph, ordrev[v] gives the first occurrence of
   *           that vertex in the order.
   *  @orig: re-order originality offset array. at a given level 'i'
   *         in the order, orig[i] is nonzero if 'i' is a repetition,
   *         in which case it holds the offset from 'i' back to the
   *         original level of the vertex in the order.
   *  @n_order: length of the re-order array.
   *  @n_orig: number of original atoms in the order.
   *  @friends: array of adjacent vertices that are earlier in the order
   *            (referred to as "friends") for every vertex in the order.
   *  @n_friends: array of friends for every vertex in the order.
   */
  unsigned int *order, *ordrev, *orig, n_order, n_orig;
  unsigned int **friends, *n_friends;

  /* miscellaneous graph properties:
   *
   *  @rmsd: deviation contribution array.
   */
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
                       value_t w, value_t *psrc, const value_semantic_t sem);

#define graph_remove_edge(G,va,vb) \
  graph_set_edge(G, va, vb, value_undefined())

void graph_count_edges (graph_t *G, unsigned int *ne, unsigned int *ni);

int graph_extend_order (graph_t *G, unsigned int v);

