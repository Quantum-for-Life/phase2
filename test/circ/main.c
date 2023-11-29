#include <stdarg.h>
#include <stdio.h>

#define TEST_INIT(desc)                                                       \
        const char* __test_desc = desc;                                       \
        int __test_res = 0;

#define TEST_FAIL(...)                                                        \
        fprintf(stderr, "FAIL: %s:%d, ", __FILE__, __LINE__);                 \
        fprintf(stderr, __VA_ARGS__);                                         \
        fprintf(stderr, "\n");                                                \
        goto __error;

#define TEST_END                                                              \
        fprintf(stderr, "PASS: %s\n", __test_desc);                           \
        goto __exit;                                                          \
        __error:  __test_res = -1;                                            \
        __exit:

#define TEST_RES __test_res

int test_case(int c) {
        TEST_INIT("test case")

        if (c < 0) {
                TEST_FAIL("something's wrong")
        }

        TEST_END
        /* Your cleanup code here */

        return TEST_RES;
}


int main() {
        int rc = 0;

        rc |= test_case(0);
        rc |= test_case(1);

        return rc;
}
