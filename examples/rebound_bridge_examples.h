/* SPDX-License-Identifier: GPL-3.0-or-later */

#ifndef REBOUND_BRIDGE_EXAMPLES_H
#define REBOUND_BRIDGE_EXAMPLES_H

/**
 * @file    rebound_bridge_examples.h
 * @brief   Convenience constructors used by the bundled examples.
 *
 * This header is NOT installed and NOT part of the public API.
 * It is only compiled when REBOUND_BRIDGE_BUILD_EXAMPLES=ON.
 */

#include "rebound_bridge.h"

#ifdef __cplusplus
extern "C" {
#endif

struct reb_simulation* reb_bridge_make_main_sun_emb(void);
struct reb_simulation* reb_bridge_make_sub_earth_moon(void);
struct reb_bridge*     reb_bridge_create_earth_moon(double dt_outer, int sub_ratio);

#ifdef __cplusplus
}
#endif

#endif /* REBOUND_BRIDGE_EXAMPLES_H */
