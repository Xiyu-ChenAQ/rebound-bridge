/* Standalone GR-potential probe — Mercury perihelion precession.
 * GR formula extracted from REBOUNDx 4.6.2 gr_potential.c (pure arithmetic,
 * no REBOUNDx dependency). Validates that the GR physics can be bridged into
 * a REBOUND 5 simulation via operator-split kicks, the same way the bridge
 * applies its cross-kicks. Expected: ~42.98 arcsec/century. */

#include <stdio.h>
#include <math.h>
#include "rebound.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define G_AU_YR_MSUN (4.0 * M_PI * M_PI)
#define C_AU_YR 63239.7263            /* speed of light, AU/yr */

/* gr_potential acceleration kick, applied for time dt_kick.
 * particles[0] is the dominant central body (the Sun). */
static void gr_potential_kick(struct reb_simulation* sim, double dt_kick) {
    struct reb_particle* p = sim->particles;
    const int N = sim->N;
    const double G = sim->G;
    const double C2 = C_AU_YR * C_AU_YR;
    const struct reb_particle source = p[0];
    const double prefac1 = 6.0 * (G * source.m) * (G * source.m) / C2;
    for (int i = 1; i < N; i++) {
        const double dx = p[i].x - source.x;
        const double dy = p[i].y - source.y;
        const double dz = p[i].z - source.z;
        const double r2 = dx*dx + dy*dy + dz*dz;
        const double prefac = prefac1 / (r2 * r2);
        p[i].vx -= prefac * dx * dt_kick;
        p[i].vy -= prefac * dy * dt_kick;
        p[i].vz -= prefac * dz * dt_kick;
        const double q = p[i].m / source.m;
        p[0].vx += q * prefac * dx * dt_kick;
        p[0].vy += q * prefac * dy * dt_kick;
        p[0].vz += q * prefac * dz * dt_kick;
    }
}

int main(void) {
    struct reb_simulation* sim = reb_simulation_create();
    if (!sim) { fprintf(stderr, "sim create failed\n"); return 1; }
    sim->G = G_AU_YR_MSUN;

    struct reb_integrator_whfast_state* wh =
        (struct reb_integrator_whfast_state*)reb_simulation_set_integrator(sim, "whfast");
    if (!wh) { fprintf(stderr, "set integrator failed\n"); return 1; }
    wh->safe_mode = 1;

    reb_simulation_add_fmt(sim, "m", 1.0);                                  /* Sun */
    reb_simulation_add_fmt(sim, "m a e", 1.66012e-7, 0.387098, 0.205630);   /* Mercury */
    reb_simulation_move_to_com(sim);

    const double dt = 0.387098 / 100.0;   /* ~1% of Mercury's period */
    const double years = 100.0;
    const long nsteps = (long)(years / dt);

    struct reb_orbit o0 = reb_orbit_from_particle(sim->G, sim->particles[1], sim->particles[0]);
    const double pomega0 = o0.pomega;

    /* Operator split: half GR kick, WHFast drift step, half GR kick. */
    for (long n = 0; n < nsteps; n++) {
        gr_potential_kick(sim, 0.5 * dt);
        reb_simulation_integrate(sim, sim->t + dt);
        gr_potential_kick(sim, 0.5 * dt);
    }

    struct reb_orbit o1 = reb_orbit_from_particle(sim->G, sim->particles[1], sim->particles[0]);
    double dpomega = o1.pomega - pomega0;
    while (dpomega >  M_PI) dpomega -= 2.0 * M_PI;
    while (dpomega < -M_PI) dpomega += 2.0 * M_PI;

    const double total_yr = sim->t;
    const double arcsec_per_century = dpomega * (180.0 / M_PI) * 3600.0 * (100.0 / total_yr);
    printf("GR-potential probe: Mercury perihelion precession\n");
    printf("  integrated %.2f yr, dt = %.5f yr\n", total_yr, dt);
    printf("  d(pomega)  = %.6e rad\n", dpomega);
    printf("  precession = %.3f arcsec/century   (GR predicts ~42.98)\n", arcsec_per_century);

    reb_simulation_free(sim);
    return 0;
}
