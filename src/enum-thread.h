
/* ensure once-only inclusion. */
#ifndef __IBPNG_ENUM_THREAD_H__
#define __IBPNG_ENUM_THREAD_H__

/* function declarations (enum-thread.c): */

int enum_threads_init (enum_t *E);

void *enum_thread_execute (void *pdata);

#endif /* !__IBPNG_ENUM_THREAD_H__ */

