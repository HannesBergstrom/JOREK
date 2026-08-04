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
   The marker weights are normalized so the marker-carried toroidal current
   equals the equilibrium I_RE exactly (`num_re` is ignored): a mismatch
   leaves eta*(j - j_RE) finite and the current channel decays resistively,
   drifting inboard.

Physics checks when preparing/updating the reference data:

* `re_eq_convergence.log` must end with max|q/q_t - 1| < 1e-3 (typically
  ~10 outer iterations).
* Current on drift orbits leaving the domain is REMOVED (linear taper of
  width `re_eq_edge_taper` in the label beyond the last closed surface):
  RE orbits crossing the wall are lost, and the marker loader applies the
  identical truncation. The last confined surface is `Ahat = 1`, whose
  `A_edge` is the extremum of `alpha_s R - psi` over the BOUNDARY nodes,
  i.e. the last drift surface that does not reach the wall. On a boundary
  that is itself a flux surface (as here) this equals the previous
  `alpha_s * max(R) - psi_bnd` exactly, so this case is unchanged; the two
  differ only on a boundary with varying psi (a vessel-shaped contour),
  where the old form truncated the beam well inside the wall by an amount
  growing linearly with class energy.
* The convergence log reports `max_edge_fraction`: the fraction of the
  CARRIED current sitting in the outermost 5% of the label range, i.e. how
  hard the beam edge is against the loss boundary. (It replaces the former
  `max_lost_fraction`, which evaluated the clipped `Nprof(1)` over every
  point outside the beam and was dominated by vacuum volume on a domain
  larger than the plasma.) When a large current fraction sits at the beam
  edge, exactly matching q_t is outside the range of the single-Nprof
  ansatz; the solver then stops at the achievable optimum and accepts it if
  below `re_eq_tol_q_soft` (with a warning), so choose the class energies of
  a regression case such that the edge fraction stays modest (as the shipped
  distribution does).
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
