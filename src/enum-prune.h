
/* ensure once-only inclusion. */
#ifndef __IBPNG_ENUM_PRUNE_H__
#define __IBPNG_ENUM_PRUNE_H__

/* function declarations (enum-prune.c): */

int enum_prune_add_closure (enum_t *E, unsigned int lev,
                            enum_prune_test_fn func,
                            void *data);

unsigned int enum_prune_get_level (unsigned int *order,
                                   unsigned int lev,
                                   unsigned int id);

/* function declarations (enum-prune-ddf.c): */

int enum_prune_ddf_init (enum_t *E, unsigned int lev);

int enum_prune_ddf (enum_t *E, enum_thread_t *th, void *data);

void enum_prune_ddf_report (enum_t *E, unsigned int lev, void *data);

/* function declarations (enum-prune-taf.c): */

int enum_prune_dihe_init (enum_t *E, unsigned int lev);

int enum_prune_impr_init (enum_t *E, unsigned int lev);

int enum_prune_taf (enum_t *E, enum_thread_t *th, void *data);

void enum_prune_dihe_report (enum_t *E, unsigned int lev, void *data);

void enum_prune_impr_report (enum_t *E, unsigned int lev, void *data);

/* function declarations (enum-prune-path.c): */

int enum_prune_path_init (enum_t *E, unsigned int lev);

int enum_prune_path (enum_t *E, enum_thread_t *th, void *data);

void enum_prune_path_report (enum_t *E, unsigned int lev, void *data);

/* function declarations (enum-prune-future.c): */

int enum_prune_future_init (enum_t *E, unsigned int lev);

int enum_prune_future (enum_t *E, enum_thread_t *th, void *data);

void enum_prune_future_report (enum_t *E, unsigned int lev, void *data);

/* function declarations (enum-prune-energy.c): */

int enum_prune_energy_init (enum_t *E, unsigned int lev);

int enum_prune_energy (enum_t *E, enum_thread_t *th, void *data);

void enum_prune_energy_report (enum_t *E, unsigned int lev, void *data);

#endif /* !__IBPNG_ENUM_PRUNE_H__ */

