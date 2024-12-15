#ifndef BUG_H
#define BUG_H

#ifndef DEBUG
#define NDEBUG
#else
#undef NDEBUG
#endif

#include <assert.h>

#define BUG_ON(x) assert(!(x))

#endif /* BUG_H */
