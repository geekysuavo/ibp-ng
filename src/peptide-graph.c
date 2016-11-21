
/* include all peptide headers. */
#include "peptide.h"
#include "peptide-atoms.h"
#include "peptide-bonds.h"
#include "peptide-angles.h"
#include "peptide-torsions.h"
#include "peptide-impropers.h"

/* peptide_graph_complete(): transform a sparsely connected graph into a
 * completely connected one by adding interval edges were no edges
 * previously existed.
 *
 * interval edge lower bounds are computed from atomic radii, and upper
 * bounds are determined using the Floyd-Warshall shortest-path algorithm.
 *
 * arguments:
 *  @P: pointer to the peptide structure to access.
 *  @G: pointer to the graph structure to modify.
 *
 * returns:
 *  integer indicating whether (1) or not (0) the operation succeeded.
 */
int peptide_graph_complete (peptide_t *P, graph_t *G) {
  /* declare required variables:
   *  @i, @j, @k: loop indices for updating bounds.
   *  @W: temporary duplicate matrix of bounds.
   */
  unsigned int i, j, k;
  value_t *W;

  /* store the graph vertex count in a local variable. */
  const unsigned int n = G->nv;

  /* allocate the edge matrix. */
  W = (value_t*) malloc(n * n * sizeof(value_t));
  if (!W)
    throw("unable to allocate edge matrix");

  /* initialize the edge matrix. */
  memcpy(W, G->E, n * n * sizeof(value_t));

  /* initialize the unknown edges. */
  for (i = 0; i < n; i++) {
    for (j = i + 1; j < n; j++) {
      /* skip existing edges. */
      if (W[i + n * j].type)
        continue;

      /* initialize the lower bound. */
      W[i + n * j].l = W[j + n * i].l =
        P->atoms[i].radius + P->atoms[j].radius;

      /* initialize the upper bound. */
      W[i + n * j].u = W[j + n * i].u = 1.0e+6;
    }
  }

  /* refine the upper bounds using shortest paths. */
  for (k = 0; k < n; k++) {
    for (i = 0; i < n; i++) {
      for (j = 0; j < n; j++) {
        /* refine the upper bound. */
        if (W[i + n * j].u > W[i + n * k].u + W[k + n * j].u)
          W[i + n * j].u = W[i + n * k].u + W[k + n * j].u;
      }
    }
  }

  /* install the computed bounds into the graph. */
  for (i = 0; i < n; i++) {
    for (j = i + 1; j < n; j++) {
      /* skip self-loops and existing edges. */
      if (i == j || graph_has_edge(G, i, j))
        continue;

      /* add the edge. */
      W[i + n * j].type = VALUE_TYPE_INTERVAL;
      graph_set_edge(G, i, j, W[i + n * j]);
    }
  }

  /* free the edge matrix. */
  free(W);

  /* return success. */
  return 1;
}

/* peptide_graph_order(): use the combined information present in a peptide
 * structure and a reorder structure to construct a repetition ordering for
 * an iDMDGP graph.
 *
 * arguments:
 *  @P: pointer to the peptide structure to access.
 *  @G: pointer to the graph structure to modify.
 *  @ord: pointer to the reorder structure to use.
 *
 * returns:
 *  integer indicating whether (1) or not (0) the operation succeeded.
 */
int peptide_graph_order (peptide_t *P, graph_t *G, reorder_t *ord) {
  /* declare required variables:
   *  @i: peptide residue index.
   *  @igrp: reorder atom index.
   *  @j: (i) atom index for reorder edge validation.
   *  @j1: (i-1) atom index for reorder edge validation.
   *  @j2: (i-2) atom index for reorder edge validation.
   *  @j3: (i-3) atom index for reorder edge validation.
   *  @visits: array for counting atom occurences in the order.
   *  @grp: reorder atom group pointer for each residue.
   *  @iatom: peptide atom index.
   */
  unsigned int i, igrp, *visits;
  unsigned int j, j1, j2, j3;
  unsigned int ires;
  reorder_t *grp;
  int iatom;

  /* allocate the visit counters. */
  visits = (unsigned int*) malloc(P->n_atoms * sizeof(unsigned int));
  if (!visits)
    throw("unable to allocate visit counter array");

  /* initialize the visit counter. */
  memset(visits, 0, P->n_atoms * sizeof(unsigned int));

  /* loop over the residues of the peptide. */
  for (i = 0; i < P->n_res; i++) {
    /* locate the necessary reorder atom group. */
    grp = reorder_get_residue(ord, peptide_get_restype(P, i));

    /* check that the reorder group was found. */
    if (!grp)
      throw("residue #%u (%s) has no re-order entry",
            i + 1, resid_get_code3(P->res[i]));

    /* loop over the atoms of the reorder group. */
    for (igrp = 0; igrp < grp->n_atoms; igrp++) {
      /* validate the residue index. */
      if (i == 0 && grp->atoms[igrp].off < 0)
        throw("invalid re-order offset for %s", grp->atoms[igrp].name);

      /* lookup the specified atom. */
      ires = i + grp->atoms[igrp].off;
      iatom = peptide_atom_find(P, ires, grp->atoms[igrp].name);

      /* check if the atom was found. */
      if (iatom < 0 && !grp->atoms[igrp].opt)
        throw("no atom '%s' found in %s%u (%s) of peptide",
              grp->atoms[igrp].name,
              peptide_get_resname(P, ires), ires + 1,
              peptide_get_restype(P, ires));

      /* add the atom into the graph order. */
      if (iatom >= 0 && !graph_extend_order(G, iatom))
        throw("unable to add atom %s%u %s to graph order",
              peptide_get_resname(P, ires), ires + 1,
              grp->atoms[igrp].name);
    }
  }

  /* loop over the repetition order entries. */
  for (i = 0; i < G->n_order; i++) {
    /* increment the visit counter. */
    j = G->order[i];
    visits[j]++;

    /* check edges to previous vertices in the order. */
    if (i >= 3) {
      /* get the three preceeding vertices in the order. */
      j1 = G->order[i - 1];
      j2 = G->order[i - 2];
      j3 = G->order[i - 3];

      /* check if this is the first visit to the vertex. */
      if (visits[j] == 1) {
        /* search for the required i,i-1 edge into the vertex. */
        if (graph_has_edge(G, j1, j) != VALUE_TYPE_SCALAR)
          raise("missing i-1,i edge from %u.%s to %u.%s",
                P->atoms[j1].res_id + 1, P->atoms[j1].name,
                P->atoms[j].res_id + 1, P->atoms[j].name);

        /* search for the required i,i-2 edge into the vertex. */
        if (graph_has_edge(G, j2, j) != VALUE_TYPE_SCALAR)
          raise("missing i-2,i edge from %u.%s to %u.%s",
                P->atoms[j2].res_id + 1, P->atoms[j2].name,
                P->atoms[j].res_id + 1, P->atoms[j].name);

        /* search for the required i,i-3 edge into the vertex. */
        if (!graph_has_edge(G, j3, j))
          raise("missing i-3,i edge from %u.%s to %u.%s",
                P->atoms[j3].res_id + 1, P->atoms[j3].name,
                P->atoms[j].res_id + 1, P->atoms[j].name);
      }
    }
  }

  /* loop over the atoms in the peptide. */
  for (i = 0; i < P->n_atoms; i++) {
    /* warn if the atom is not visited in the order. */
    if (visits[i] == 0)
      warn("atom '%s' in %s%u (%s) is not in the graph order",
           P->atoms[i].name,
           peptide_get_resname(P, P->atoms[i].res_id),
           P->atoms[i].res_id + 1,
           peptide_get_restype(P, P->atoms[i].res_id));
  }

  /* free the visit counter. */
  free(visits);

  /* check if errors were raised during graph re-order validation. */
  if (traceback_length())
    throw("graph is missing required edges");

  /* return success. */
  return 1;
}

/* peptide_graph(): use the information contained in a peptide structure
 * to add an initial set of vertices and edges to a graph.
 *
 * arguments:
 *  @P: pointer to the peptide structure to access.
 *  @ord: pointer to the reorder structure to access.
 *
 * returns:
 *  pointer to a newly allocated, initialized and filled graph structure,
 *  or NULL on failure.
 */
graph_t *peptide_graph (peptide_t *P, reorder_t *ord) {
  /* declare required variables:
   *  @G: pointer to a new graph structure.
   */
  graph_t *G;

  /* check that the structure pointers are valid. */
  if (!P || !ord) {
    /* nope. raise an exception and return null. */
    raise("input structure pointers are null");
    return NULL;
  }

  /* check that the peptide structure is completed. */
  if (P->n_atoms == 0) {
    /* nope. raise an exception and return null. */
    raise("peptide structure is incomplete");
    return NULL;
  }

  /* allocate a new graph structure for holding peptide information. */
  G = graph_new(P->n_atoms);

  /* check if graph allocation failed. */
  if (!G) {
    /* yes. raise an exception and return null. */
    raise("unable to allocate new graph");
    return NULL;
  }

  /* convert bonds to graph edges. */
  if (!peptide_graph_bonds(P, G)) {
    /* raise an exception and return null. */
    raise("unable to update bond-derived graph edges");
    graph_free(G);
    return NULL;
  }

  /* convert angles to graph edges. */
  if (!peptide_graph_angles(P, G)) {
    /* raise an exception and return null. */
    raise("unable to update angle-derived graph edges");
    graph_free(G);
    return NULL;
  }

  /* convert torsions to graph edges. */
  if (!peptide_graph_torsions(P, G)) {
    /* raise an exception and return null. */
    raise("unable to update dihedral-derived graph edges");
    graph_free(G);
    return NULL;
  }

  /* convert impropers to graph edges. */
  if (!peptide_graph_impropers(P, G)) {
    /* raise an exception and return null. */
    raise("unable to update improper-derived graph edges");
    graph_free(G);
    return NULL;
  }

  /* construct the graph order. */
  if (!peptide_graph_order(P, G, ord)) {
    /* raise an exception and return null. */
    raise("unable to build graph re-order");
    graph_free(G);
    return NULL;
  }

  /* complete the graph edge set. */
  if (!peptide_graph_complete(P, G)) {
    /* raise an exception and return null. */
    raise("unable to complete graph edge set");
    graph_free(G);
    return NULL;
  }

  /* return the final graph. */
  return G;
}

