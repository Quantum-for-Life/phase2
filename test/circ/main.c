#include <stdarg.h>
#include <stdio.h>

#include "test.h"

TEST(main, void)

        TEST_ASSERT(0 == 1, "Not good enough")

        TEST_FIN
TEST_END
