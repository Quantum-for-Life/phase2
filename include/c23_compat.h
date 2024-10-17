#ifndef C23_COMPAT_H
#define C23_COMPAT_H

/* --------------------------------------------------------------------------*/
/* Remove this section when C23 arrives.                                     */
#include <stdbool.h>
#define nullptr (void *)0
#define stdc_count_ones_ul(v) (__builtin_popcountl(v))
#define unreachable() (__builtin_unreachable())
/* --------------------------------------------------------------------------*/
/* Uncomment this section when C23 arrives.                                  */
/*
#include <stdbit.h>
*/
/* --------------------------------------------------------------------------*/

#endif /* C23_COMPAT_H */
