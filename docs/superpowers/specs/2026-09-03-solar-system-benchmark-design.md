# Solar System Bridge Benchmark Design

The benchmark will compare Bridge EOS LF(4,2), global WHFast, and IAS15 on a fixed, reproducible Solar System fixture. The fixture contains the Sun, eight planets, five named dwarf planets, the Moon, Phobos/Deimos, four Galilean moons, and four major Saturnian moons. Initial Cartesian states are deterministic analytic/JPL-inspired values stored in the benchmark source; no runtime network access is required.

Each method integrates 100 years and writes sampled Cartesian states and wall-clock timings to CSV. IAS15 is the reference trajectory. Bridge uses a 7-day outer step and subsystem inner steps selected from each satellite group; WHFast uses a conservative 0.25-day global step. Python then computes per-body relative position error (normalised by the IAS15 heliocentric distance), maximum/RMS summaries, and speedups, producing PNG figures and a JSON summary.

The benchmark is intentionally not a cost-equivalent comparison: Bridge's multirate outer step is fixed at 7 days as requested, while WHFast uses its required small global step.
