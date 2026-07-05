/* SPDX-License-Identifier: GPL-3.0-or-later */
#ifndef VAL_CLI_H
#define VAL_CLI_H

/**
 * @file  val_cli.h
 * @brief Shared CLI parsing utilities for validation programs.
 */

#include "rebound_bridge.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Parse a positive double from text.  Returns 0 on success, -1 on error. */
int val_parse_double(const char* text, double* out);

/* Parse a positive int from text.  Returns 0 on success, -1 on error. */
int val_parse_int(const char* text, int* out);

/* Parse a non-negative int from text.  Returns 0 on success, -1 on error. */
int val_parse_nonnegative_int(const char* text, int* out);

/*
 * Parse a bridge method name ("order2", "eos-lf4-2", "eos-lf8-6-4",
 * "forward-fg4") into the corresponding enum value.
 * Returns 0 on success, -1 on unknown name.
 */
int val_parse_bridge_method(const char* text, enum reb_bridge_method* out);

/* Return the canonical CLI name for a bridge method. */
const char* val_bridge_method_name(enum reb_bridge_method method);

#ifdef __cplusplus
}
#endif

#endif /* VAL_CLI_H */
