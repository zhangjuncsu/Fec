#include "util.h"

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>

void
assertion_fail_handler(const char* expr, const char* file, const int line, const char* fmt, ...)
{
    fprintf(stderr, "Assertion Fail At '%s: %d'\n", file, line);
    fprintf(stderr, "Expression: '%s'\n", expr);
    if (!fmt) return;
    fprintf(stderr, "Context Information: '");
    va_list ap;
    va_start(ap, fmt);
    vfprintf(stderr, fmt, ap);
    va_end(ap);
    fprintf(stderr, "'\n");
    exit(1);
}

void
oc_dump_info(FILE* out, 
        const char* level, 
        const char* file, 
        const int line, 
        const char* fmt, 
        ...)
{
    fprintf(out, "%s: ", level);
    va_list ap;
    va_start(ap, fmt);
    vfprintf(out, fmt, ap);
    va_end(ap);
    if (file) fprintf(out, " (%s, %d)", file, line);
    fprintf(out, "\n");
}