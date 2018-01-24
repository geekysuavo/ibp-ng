
/* ensure once-only inclusion. */
#pragma once

/* function declarations (enum-thread.c): */

int enum_threads_init (enum_t *E);

void *enum_thread_timer (void *pdata);

void *enum_thread_execute (void *pdata);

