/**
 * solar_system_bridge.c - Multi-subsystem rebound_bridge validation example.
 *
 * Main simulation: Sun, eight planets, Earth-Moon host, and Jupiter-system host.
 * Subsystems: resolved Earth-Moon and Jupiter-Io-Europa-Ganymede systems.
 *
 * Initial vectors are fixed JPL Horizons geometric states for
 * 2026-Jan-01 00:00 TDB.  They make this example deterministic and suitable
 * for validating bridge plumbing; this is not a full ephemeris-grade model.
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "rebound_bridge.h"

enum main_body {
    MAIN_SUN = 0,
    MAIN_MERCURY = 1,
    MAIN_VENUS = 2,
    MAIN_EARTH_MOON_HOST = 3,
    MAIN_MARS = 4,
    MAIN_JUPITER_HOST = 5,
    MAIN_SATURN = 6,
    MAIN_URANUS = 7,
    MAIN_NEPTUNE = 8
};

enum jovian_body {
    JOVIAN_JUPITER = 0,
    JOVIAN_IO = 1,
    JOVIAN_EUROPA = 2,
    JOVIAN_GANYMEDE = 3
};

struct range_stats {
    double min;
    double max;
};

static struct range_stats range_init(double value) {
    struct range_stats r;
    r.min = value;
    r.max = value;
    return r;
}

static void range_add(struct range_stats* r, double value) {
    if (value < r->min) r->min = value;
    if (value > r->max) r->max = value;
}

static double wrap_degrees(double angle) {
    angle = fmod(angle, 360.0);
    if (angle < 0.0) angle += 360.0;
    return angle;
}

static double centered_delta_degrees(double angle, double center) {
    double delta = angle - center;
    while (delta <= -180.0) delta += 360.0;
    while (delta > 180.0) delta -= 360.0;
    return delta;
}

static double distance_particles(const struct reb_particle* a, const struct reb_particle* b) {
    const double dx = a->x - b->x;
    const double dy = a->y - b->y;
    const double dz = a->z - b->z;
    return sqrt(dx * dx + dy * dy + dz * dz);
}

static void add_cartesian(
    struct reb_simulation* sim,
    double mass,
    double x,
    double y,
    double z,
    double vx,
    double vy,
    double vz
) {
    reb_simulation_add_fmt(sim, "m x y z vx vy vz", mass, x, y, z, vx, vy, vz);
}

static struct reb_simulation* make_main_solar_system(void) {
    struct reb_simulation* sim = reb_simulation_create();
    if (!sim) return NULL;

    sim->G = REB_BRIDGE_G_AU_YR_MSUN;
    sim->integrator = REB_INTEGRATOR_WHFAST;
    sim->ri_whfast.safe_mode = 1;

    const double m_mercury = 1.660120e-7;
    const double m_venus = 2.4478383e-6;
    const double m_earth_moon = 3.0034896e-6 + 3.6940e-8;
    const double m_mars = 3.227151e-7;
    const double m_jupiter_system = 9.547919e-4 + 2.388e-8 + 2.527e-8 + 7.804e-8;
    const double m_saturn = 2.8588567e-4;
    const double m_uranus = 4.366244e-5;
    const double m_neptune = 5.151389e-5;

    reb_simulation_add_fmt(sim, "m", 1.0);
    add_cartesian(sim, m_mercury, -2.15201304377667713e-01, -4.09207615700272087e-01, -1.37032695919236600e-02, 7.02443404629697010e+00, -4.28716799788450942e+00, -9.94632701667513053e-01);
    add_cartesian(sim, m_venus, 8.88772496163565545e-02, -7.21762382317080675e-01, -1.50440565750619305e-02, 7.28255910865529810e+00, 8.76733642322772511e-01, -4.08156573658943866e-01);
    add_cartesian(sim, m_earth_moon, -1.74269758548305093e-01, 9.67785674324191270e-01, -5.68660813624550973e-05, -6.28653514426191418e+00, -1.13726197010495311e+00, 1.00945156455726687e-04);
    add_cartesian(sim, m_mars, 3.40579676862215075e-01, -1.38700201594525407e+00, -3.74172267877010828e-02, 5.15711521620981550e+00, 1.65830921706158052e+00, -9.17089273335817273e-02);
    add_cartesian(sim, m_jupiter_system, -1.69400303018584597e+00, 4.92888271463623262e+00, 1.74262399033241998e-02, -2.64098245722763814e+00, -7.67587434003239921e-01, 6.22754635347801122e-02);
    add_cartesian(sim, m_saturn, 9.50734336522187107e+00, 2.57738457802152210e-01, -3.82935528507753087e-01, -1.69141631745304527e-01, 2.03249650346815258e+00, -2.86682018818476998e-02);
    add_cartesian(sim, m_uranus, 9.88031797519427535e+00, 1.68000154098179486e+01, -6.57197991373143220e-02, -1.25149298434788658e+00, 6.61240423623731699e-01, 1.86840216329132620e-02);
    add_cartesian(sim, m_neptune, 2.98721200408889800e+01, 5.18939463241125942e-01, -6.99067175203629043e-01, -2.98764379233138963e-02, 1.15313346010349438e+00, -2.30704057616862752e-02);

    reb_simulation_move_to_com(sim);
    return sim;
}

static struct reb_simulation* make_sub_earth_moon_horizons(void) {
    struct reb_simulation* sim = reb_simulation_create();
    if (!sim) return NULL;

    sim->G = REB_BRIDGE_G_AU_YR_MSUN;
    sim->integrator = REB_INTEGRATOR_WHFAST;
    sim->ri_whfast.safe_mode = 1;

    add_cartesian(sim, 3.0034896e-6, -1.17223722782580293e-05, -2.67540695880473792e-05, -2.57902881285356581e-06, 2.57421785107276755e-03, -1.07874487935336202e-03, -1.43157903344392045e-05);
    add_cartesian(sim, 3.6940e-8, 9.53035527126304600e-04, 2.17512105974572998e-03, 2.09676507944608496e-04, -2.09285374018137460e-01, 8.77025716574588027e-02, 1.16388188872971468e-03);

    reb_simulation_move_to_com(sim);
    return sim;
}

static struct reb_simulation* make_sub_jovian_horizons(void) {
    struct reb_simulation* sim = reb_simulation_create();
    if (!sim) return NULL;

    sim->G = REB_BRIDGE_G_AU_YR_MSUN;
    sim->integrator = REB_INTEGRATOR_WHFAST;
    sim->ri_whfast.safe_mode = 1;

    reb_simulation_add_fmt(sim, "m", 9.547919e-4);
    add_cartesian(sim, 2.388e-8, 2.48120485765584205e-03, -1.34556438159649208e-03, -1.31630007347439604e-05, 1.75341752994026812e+00, 3.20163492363371871e+00, 1.39542430489110164e-01);
    add_cartesian(sim, 2.527e-8, 5.51813635549570501e-04, -4.47551321186325738e-03, -1.27610904592036212e-04, 2.85755306406471643e+00, 3.72439324002732552e-01, 7.16871253077594472e-02);
    add_cartesian(sim, 7.804e-8, 6.76397567656043287e-03, -2.29852901298470890e-03, 1.01873701583559195e-05, 7.40898140660848692e-01, 2.17447712709783314e+00, 9.37799108604124437e-02);

    reb_simulation_move_to_com(sim);
    return sim;
}

static struct reb_bridge* create_solar_system_bridge(
    double dt_outer,
    double dt_earth_moon,
    double dt_jovian
) {
    struct reb_simulation* main_sim = make_main_solar_system();
    struct reb_simulation* earth_moon = make_sub_earth_moon_horizons();
    struct reb_simulation* jovian = make_sub_jovian_horizons();
    if (!main_sim || !earth_moon || !jovian) {
        if (main_sim) reb_simulation_free(main_sim);
        if (earth_moon) reb_simulation_free(earth_moon);
        if (jovian) reb_simulation_free(jovian);
        return NULL;
    }

    struct reb_bridge* bridge = reb_bridge_create(main_sim, dt_outer);
    if (!bridge) {
        reb_simulation_free(main_sim);
        reb_simulation_free(earth_moon);
        reb_simulation_free(jovian);
        return NULL;
    }

    reb_bridge_set_owns_main(bridge, 1);
    if (reb_bridge_add_subsystem(bridge, MAIN_EARTH_MOON_HOST, earth_moon, dt_earth_moon, 1) != 0) {
        reb_bridge_free(bridge);
        return NULL;
    }
    if (reb_bridge_add_subsystem(bridge, MAIN_JUPITER_HOST, jovian, dt_jovian, 1) != 0) {
        reb_bridge_free(bridge);
        return NULL;
    }
    return bridge;
}

static double host_distance(const struct reb_simulation* main_sim, int host_index) {
    return distance_particles(&main_sim->particles[MAIN_SUN], &main_sim->particles[host_index]);
}

static struct reb_orbit jovian_orbit(const struct reb_simulation* jovian, int moon_index) {
    return reb_orbit_from_particle(
        jovian->G,
        jovian->particles[moon_index],
        jovian->particles[JOVIAN_JUPITER]
    );
}

static double jovian_laplace_angle_degrees(const struct reb_simulation* jovian) {
    const struct reb_orbit io = jovian_orbit(jovian, JOVIAN_IO);
    const struct reb_orbit europa = jovian_orbit(jovian, JOVIAN_EUROPA);
    const struct reb_orbit ganymede = jovian_orbit(jovian, JOVIAN_GANYMEDE);
    return wrap_degrees((io.l - 3.0 * europa.l + 2.0 * ganymede.l) * 180.0 / M_PI);
}

int main(void) {
    const double dt_outer = 1.0 / 365.25;
    const double dt_earth_moon = dt_outer / 20.0;
    const double dt_jovian = dt_outer / 50.0;
    const double years = 100.0;
    const int samples = 1000;

    struct reb_bridge* bridge = create_solar_system_bridge(dt_outer, dt_earth_moon, dt_jovian);
    if (!bridge) {
        fprintf(stderr, "failed to create solar-system bridge\n");
        return EXIT_FAILURE;
    }

    if (reb_bridge_subsystem_count(bridge) != 2) {
        fprintf(stderr, "expected two subsystems\n");
        reb_bridge_free(bridge);
        return EXIT_FAILURE;
    }

    struct reb_simulation* main_sim = reb_bridge_main_sim(bridge);
    struct reb_simulation* earth_moon = reb_bridge_sub_sim(bridge, 0);
    struct reb_simulation* jovian = reb_bridge_sub_sim(bridge, 1);

    double earth_moon_sum = 0.0;
    double earth_moon_sum2 = 0.0;
    struct range_stats earth_host = range_init(host_distance(main_sim, MAIN_EARTH_MOON_HOST));
    struct range_stats jupiter_host = range_init(host_distance(main_sim, MAIN_JUPITER_HOST));
    struct range_stats laplace_delta = range_init(
        centered_delta_degrees(jovian_laplace_angle_degrees(jovian), 180.0)
    );

    for (int i = 0; i < samples; i++) {
        const double target = years * (double)(i + 1) / (double)samples;
        if (reb_bridge_advance(bridge, target - main_sim->t) != 0) {
            fprintf(stderr, "advance failed at sample %d\n", i + 1);
            reb_bridge_free(bridge);
            return EXIT_FAILURE;
        }

        const double em_distance = distance_particles(&earth_moon->particles[0], &earth_moon->particles[1]);
        earth_moon_sum += em_distance;
        earth_moon_sum2 += em_distance * em_distance;
        range_add(&earth_host, host_distance(main_sim, MAIN_EARTH_MOON_HOST));
        range_add(&jupiter_host, host_distance(main_sim, MAIN_JUPITER_HOST));
        range_add(&laplace_delta, centered_delta_degrees(jovian_laplace_angle_degrees(jovian), 180.0));
    }

    const double mean_em = earth_moon_sum / (double)samples;
    const double var_em = earth_moon_sum2 / (double)samples - mean_em * mean_em;
    const double std_em = sqrt(var_em > 0.0 ? var_em : 0.0);
    const double max_abs_laplace = fmax(fabs(laplace_delta.min), fabs(laplace_delta.max));

    printf("rebound_bridge - solar-system multi-subsystem demo\n");
    printf("==================================================\n");
    printf("  main bodies       = %d\n", main_sim->N);
    printf("  subsystems        = %d\n", reb_bridge_subsystem_count(bridge));
    printf("  integration span  = %.1f yr\n", years);
    printf("  dt_outer          = %.3f days\n", dt_outer * 365.25);
    printf("  dt_earth_moon     = %.3f days\n", dt_earth_moon * 365.25);
    printf("  dt_jovian         = %.3f days\n", dt_jovian * 365.25);
    printf("  Earth-Moon d_mean = %.8f AU\n", mean_em);
    printf("  Earth-Moon d_std  = %.3e AU\n", std_em);
    printf("  Earth host r      = %.6f .. %.6f AU\n", earth_host.min, earth_host.max);
    printf("  Jupiter host r    = %.6f .. %.6f AU\n", jupiter_host.min, jupiter_host.max);
    printf("  Jovian phi_L-180  = %.6f .. %.6f deg (max_abs %.6f)\n",
           laplace_delta.min,
           laplace_delta.max,
           max_abs_laplace);

    reb_bridge_free(bridge);
    return max_abs_laplace < 90.0 ? EXIT_SUCCESS : EXIT_FAILURE;
}
