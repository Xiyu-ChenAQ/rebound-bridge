/* SPDX-License-Identifier: GPL-3.0-or-later */
#include "val_cli.h"

#include <errno.h>
#include <stdlib.h>
#include <string.h>

int val_parse_double(const char* text, double* out) {
    char* end = NULL;
    errno = 0;
    const double value = strtod(text, &end);
    if (errno != 0 || end == text || *end != '\0') return -1;
    *out = value;
    return 0;
}

int val_parse_int(const char* text, int* out) {
    char* end = NULL;
    errno = 0;
    const long value = strtol(text, &end, 10);
    if (errno != 0 || end == text || *end != '\0' || value <= 0 || value > 2147483647L) return -1;
    *out = (int)value;
    return 0;
}

int val_parse_nonnegative_int(const char* text, int* out) {
    char* end = NULL;
    errno = 0;
    const long value = strtol(text, &end, 10);
    if (errno != 0 || end == text || *end != '\0' || value < 0 || value > 2147483647L) return -1;
    *out = (int)value;
    return 0;
}

int val_parse_bridge_method(const char* text, enum reb_bridge_method* out) {
    if (strcmp(text, "order2")       == 0) { *out = REB_BRIDGE_METHOD_ORDER2;      return 0; }
    if (strcmp(text, "eos-lf4-2")   == 0) { *out = REB_BRIDGE_METHOD_EOS_LF4_2;   return 0; }
    if (strcmp(text, "eos-lf8-6-4") == 0) { *out = REB_BRIDGE_METHOD_EOS_LF8_6_4; return 0; }
    if (strcmp(text, "forward-fg4") == 0) { *out = REB_BRIDGE_METHOD_FORWARD_FG4;  return 0; }
    return -1;
}

const char* val_bridge_method_name(enum reb_bridge_method method) {
    switch (method) {
        case REB_BRIDGE_METHOD_ORDER2:      return "order2";
        case REB_BRIDGE_METHOD_EOS_LF4_2:   return "eos-lf4-2";
        case REB_BRIDGE_METHOD_EOS_LF8_6_4: return "eos-lf8-6-4";
        case REB_BRIDGE_METHOD_FORWARD_FG4:  return "forward-fg4";
        case REB_BRIDGE_METHOD_INVALID:      return "invalid";
    }
    return "unknown";
}
