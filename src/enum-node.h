
/* ensure once-only inclusion. */
#ifndef __IBPNG_ENUM_NODE_H__
#define __IBPNG_ENUM_NODE_H__

/* function declarations (enum-node.c): */

void enum_node_init (enum_node_t *node, enum_node_t *prev);

void enum_node_free (enum_t *E, enum_node_t *node);

int enum_node_set_branches (enum_node_t *node, unsigned int nb);

int enum_node_compute (enum_t *E, enum_node_t *node);

#endif /* !__IBPNG_ENUM_NODE_H__ */

