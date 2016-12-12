
/* ensure once-only inclusion. */
#ifndef __IBPNG_ENUM_QUEUE_H__
#define __IBPNG_ENUM_QUEUE_H__

/* function declarations (enum-queue.c): */

int enum_queue_insert (enum_t *E, enum_node_t *node);

enum_node_t *enum_queue_pull (enum_t *E);

#endif /* !__IBPNG_ENUM_QUEUE_H__ */

