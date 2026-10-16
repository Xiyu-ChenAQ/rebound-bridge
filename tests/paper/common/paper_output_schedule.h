#ifndef PAPER_OUTPUT_SCHEDULE_H
#define PAPER_OUTPUT_SCHEDULE_H

#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>

struct paper_output_schedule {
    uint64_t total_steps;
    uint64_t steps_per_output;
    double dt_years;
    double output_duration_years;
};

/* Accept only a whole number of outer steps and evenly spaced outputs on
 * that grid. The tolerance accounts for decimal-input arithmetic, not a
 * physical partial step; no final step or output interval is rounded down. */
static int paper_make_output_schedule(double years, double dt_days, int samples,
                                      struct paper_output_schedule* schedule) {
    const double steps = years * 365.25 / dt_days;
    const double rounded = round(steps);
    if (!isfinite(years) || !isfinite(dt_days) || !(years > 0.0) ||
        !(dt_days > 0.0) || samples < 1 || !isfinite(steps) ||
        rounded < 1.0 || rounded > 9007199254740991.0 ||
        fabs(steps - rounded) > fmin(0.25, 64.0 * DBL_EPSILON * fmax(1.0, steps))) {
        fprintf(stderr, "requested duration must be an integer number of outer steps\n");
        return -1;
    }
    schedule->total_steps = (uint64_t)rounded;
    if (schedule->total_steps % (uint64_t)samples != 0) {
        fprintf(stderr, "requested samples must divide the total number of outer steps\n");
        return -1;
    }
    schedule->steps_per_output = schedule->total_steps / (uint64_t)samples;
    schedule->dt_years = dt_days / 365.25;
    schedule->output_duration_years =
        (double)schedule->steps_per_output * schedule->dt_years;
    if (!(schedule->dt_years > 0.0) || !isfinite(schedule->output_duration_years)) {
        fprintf(stderr, "output schedule is not representable\n");
        return -1;
    }
    return 0;
}

#endif
