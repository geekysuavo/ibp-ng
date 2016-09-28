
# To-do list for ibp-ng

Just a small document for recording the tasks that need tasking and the
dos that need doing.

## Small tasks

 * Possible to more intelligently precompute discretization points?
 * Compute write orderings from the re-order that preserve residue order.
 * Beneficial to keep sign of dihedrals when known, best way to implement?

## Big tasks

 * Energetic pruning for true global optimization.
    - Conversion of force fields and restraints into parameters of
      appropriate energy functionals.
    - Pruning and best-first search based on evaluated energy.
 * Support for GPU-based tree traversal. [back burner]

## Misc. notes

 * Moving from literal trees to indices enables multi-threading.
 * Flattening vertices and matrices into the nodes has no measurable
   effect on performance.
 * The search order in graph_get_edge() is *CRITICAL* to performance.
 * Binary search is definitely faster: 4m vs 18m for the test case.
 * Inlining graph_cmp_{edge,vertex}() makes absolutely no difference.
 * Interleaving (I1,I2) and (E1,E2) into Iv and Ev, resp., is slightly
   slower: 82.04% vs. 81.63% usage by graph_get_edge().
 * Ugh, facepalm time: moved from sparse to dense storage of graph edges.
   Edge lookups have become array accesses (fast!) and since the VDW/DSP
   edges will *always* make the edge set dense, there was no advantage to
   using sparse storage.

