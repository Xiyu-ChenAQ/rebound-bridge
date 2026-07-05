/* SPDX-License-Identifier: GPL-3.0-or-later */
#include "solar_system_fixture.h"
#include "val_particle.h"

#include <stdlib.h>

static void add_cartesian(struct reb_simulation* sim,
    double m, double x, double y, double z,
    double vx, double vy, double vz)
{
    reb_simulation_add_fmt(sim, "m x y z vx vy vz", m, x, y, z, vx, vy, vz);
}

void ss_populate_bodies(struct reb_particle bodies[SS_BODY_COUNT]) {
    bodies[SS_BODY_SUN]    = val_particle_from_state(1.0, 0,0,0, 0,0,0);
    bodies[SS_BODY_MERCURY]= val_particle_from_state(1.660120e-7,
        -2.15201304377667713e-01,-4.09207615700272087e-01,-1.37032695919236600e-02,
         7.02443404629697010e+00,-4.28716799788450942e+00,-9.94632701667513053e-01);
    bodies[SS_BODY_VENUS]  = val_particle_from_state(2.4478383e-6,
         8.88772496163565545e-02,-7.21762382317080675e-01,-1.50440565750619305e-02,
         7.28255910865529810e+00, 8.76733642322772511e-01,-4.08156573658943866e-01);
    bodies[SS_BODY_EARTH]  = val_particle_from_state(3.0034896e-6,
        -1.74269758548305093e-01, 9.67785674324191270e-01,-5.68660813624550973e-05,
        -6.28653514426191418e+00,-1.13726197010495311e+00, 1.00945156455726687e-04);
    bodies[SS_BODY_MOON]   = val_particle_from_state(3.6940e-8,
        -1.74269758548305093e-01+9.53035527126304600e-04,
         9.67785674324191270e-01+2.17512105974572998e-03,
        -5.68660813624550973e-05+2.09676507944608496e-04,
        -6.28653514426191418e+00-2.09285374018137460e-01,
        -1.13726197010495311e+00+8.77025716574588027e-02,
         1.00945156455726687e-04+1.16388188872971468e-03);
    bodies[SS_BODY_MARS]   = val_particle_from_state(3.227151e-7,
         3.40579676862215075e-01,-1.38700201594525407e+00,-3.74172267877010828e-02,
         5.15711521620981550e+00, 1.65830921706158052e+00,-9.17089273335817273e-02);
    bodies[SS_BODY_JUPITER]= val_particle_from_state(9.547919e-4,
        -1.69400303018584597e+00, 4.92888271463623262e+00, 1.74262399033241998e-02,
        -2.64098245722763814e+00,-7.67587434003239921e-01, 6.22754635347801122e-02);
    bodies[SS_BODY_IO]     = val_particle_from_state(2.388e-8,
        -1.69400303018584597e+00+2.48120485765584205e-03,
         4.92888271463623262e+00-1.34556438159649208e-03,
         1.74262399033241998e-02-1.31630007347439604e-05,
        -2.64098245722763814e+00+1.75341752994026812e+00,
        -7.67587434003239921e-01+3.20163492363371871e+00,
         6.22754635347801122e-02+1.39542430489110164e-01);
    bodies[SS_BODY_EUROPA] = val_particle_from_state(2.527e-8,
        -1.69400303018584597e+00+5.51813635549570501e-04,
         4.92888271463623262e+00-4.47551321186325738e-03,
         1.74262399033241998e-02-1.27610904592036212e-04,
        -2.64098245722763814e+00+2.85755306406471643e+00,
        -7.67587434003239921e-01+3.72439324002732552e-01,
         6.22754635347801122e-02+7.16871253077594472e-02);
    bodies[SS_BODY_GANYMEDE]=val_particle_from_state(7.804e-8,
        -1.69400303018584597e+00+6.76397567656043287e-03,
         4.92888271463623262e+00-2.29852901298470890e-03,
         1.74262399033241998e-02+1.01873701583559195e-05,
        -2.64098245722763814e+00+7.40898140660848692e-01,
        -7.67587434003239921e-01+2.17447712709783314e+00,
         6.22754635347801122e-02+9.37799108604124437e-02);
    bodies[SS_BODY_SATURN] = val_particle_from_state(2.8588567e-4,
         9.50734336522187107e+00, 2.57738457802152210e-01,-3.82935528507753087e-01,
        -1.69141631745304527e-01, 2.03249650346815258e+00,-2.86682018818476998e-02);
    bodies[SS_BODY_URANUS] = val_particle_from_state(4.366244e-5,
         9.88031797519427535e+00, 1.68000154098179486e+01,-6.57197991373143220e-02,
        -1.25149298434788658e+00, 6.61240423623731699e-01, 1.86840216329132620e-02);
    bodies[SS_BODY_NEPTUNE]= val_particle_from_state(5.151389e-5,
         2.98721200408889800e+01, 5.18939463241125942e-01,-6.99067175203629043e-01,
        -2.98764379233138963e-02, 1.15313346010349438e+00,-2.30704057616862752e-02);
}

void ss_add_all_bodies(struct reb_simulation* sim) {
    struct reb_particle b[SS_BODY_COUNT];
    ss_populate_bodies(b);
    for (int i = 0; i < SS_BODY_COUNT; i++) {
        add_cartesian(sim, b[i].m, b[i].x, b[i].y, b[i].z,
                           b[i].vx, b[i].vy, b[i].vz);
    }
}

struct reb_simulation* ss_make_whfast(double dt_days) {
    struct reb_simulation* sim = reb_simulation_create();
    sim->G = REB_BRIDGE_G_AU_YR_MSUN;
    if (reb_bridge_set_integrator_whfast(sim) != 0) {
        reb_simulation_free(sim);
        return NULL;
    }
    sim->dt = dt_days / 365.25;
    sim->exact_finish_time = 1;
    ss_add_all_bodies(sim);
    reb_simulation_move_to_com(sim);
    return sim;
}

struct reb_simulation* ss_make_ias15(void) {
    struct reb_simulation* sim = reb_simulation_create();
    sim->G = REB_BRIDGE_G_AU_YR_MSUN;
    sim->integrator = REB_INTEGRATOR_IAS15;
    sim->exact_finish_time = 1;
    ss_add_all_bodies(sim);
    reb_simulation_move_to_com(sim);
    return sim;
}
