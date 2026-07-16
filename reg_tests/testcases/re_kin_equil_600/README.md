# re_kin_equil_600: kinetic RE drift-surface equilibrium test

Tests the `re_kinetic_equilibrium` feature end to end with model 600 and the
`kinetic_main` workflow:

1. `jorek_model600` (restart=.f.) computes the pure-RE fixed-boundary
   equilibrium in which the current of the three RE classes of
   `re_distribution.dat` lies on drift (constant canonical toroidal momentum)
   surfaces, with the common profile function iterated until q matches
   `qprofile_target.dat`. Outputs: the standard fluid restart
   (`jorek00000.h5`), the per-class handoff file `re_equilibrium.dat`, and
   `re_eq_convergence.log`.
2. `kinetic_main` (restart=.t., restart_particles=.f.) initializes the
   full-orbit RE markers from the stationary per-class density
   `n_s ~ w_s Nprof(Ahat_s)/R` (init_function='equilibrium') and runs a few
   coupled steps; the regression run restarts from `part_restart.h5`.

Physics checks when preparing/updating the reference data:

* `re_eq_convergence.log` must end with max|q/q_t - 1| < 1e-3 (typically
  ~10 outer iterations).
* `qprofile.dat` of the equilibrium run must agree with
  `qprofile_target.dat`.
* The per-class drift-axis positions in `re_equilibrium.dat` (columns
  R_axis) must shift outboard with class energy: at B0 = 0.3 T the 1, 5 and
  10 MeV classes have drift parameters d_s of about 1.2e-2, 5.4e-2 and
  1.1e-1, i.e. axis shifts of order q0*d_s*a ~ 1.5, 7 and 14 cm.
* During the coupled `kinetic_main` phase, the projected RE current and the
  q profile must stay stationary at marker-noise level (the whole point of
  initializing in invariant space); compare the binned j_phi at t=0 and
  after the run.

See models/mod_re_kinetic_equilibrium.f90 for the physics documentation and
util/re_equilibrium_prototype/ for the standalone reference implementation.
