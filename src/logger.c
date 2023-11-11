#include "logger.h"
#include <stdio.h>
#include <time.h>

void logger(const char *tag, const char *message) {
    char buff[100];
    time_t now;
    time(&now);
    strftime(buff, 100, "%Y-%m-%d %H:%M:%S.000", localtime(&now));
    fprintf(stderr, "%s [%s]: %s\n", buff, tag, message);
}