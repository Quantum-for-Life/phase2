#include <stdio.h>
#include <time.h>
#include <stdarg.h>

void log_format(const char *tag, const char *message, va_list args) {
    char date[100];
    time_t now;
    time(&now);
    strftime(date, 100, "%Y-%m-%d %H:%M:%S.000", localtime(&now));
    fprintf(stderr, "%s [%s] ", date, tag);
    vfprintf(stderr, message, args);
    fprintf(stderr, "\n");
}

void log_error(const char *message, ...) {
    va_list args;
    va_start(args, message);
    log_format("error", message, args);
    va_end(args);
}

void log_info(const char *message, ...) {
    va_list args;
    va_start(args, message);
    log_format("info", message, args);
    va_end(args);
}

void log_debug(const char *message, ...) {
    va_list args;
    va_start(args, message);
    log_format("debug", message, args);
    va_end(args);
}
