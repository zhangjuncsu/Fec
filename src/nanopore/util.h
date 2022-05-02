#pragma once

#include <stdio.h>

void assertion_fail_handler(const char* expr, const char* file, const int line, const char* fmt, ...);

#define __oc_assert(expr, ...) \
    do { \
        if (!(expr)) { \
            assertion_fail_handler(#expr, __VA_ARGS__, NULL); \
            exit(1); \
        } \
    } while(0)

#define oc_assert(expr, args...) __oc_assert(expr, __FILE__, __LINE__, ##args);

void 
oc_dump_info(FILE* out, 
            const char* level, 
            const char* file, 
            const int line, 
            const char* fmt, 
            ...);
#define OC_ERROR(fmt, args...) do { oc_dump_info(stderr, "ERROR", __FILE__, __LINE__, fmt, ##args); exit(1); } while(0)

#define OC_MIN(a, b) ((a) < (b) ? (a) : (b))
#define OC_MAX(a, b) ((a) > (b) ? (a) : (b))