# re_kin_equil_xpoint_600: kinetic RE drift-surface equilibrium, DIVERTED (X-point)

**Status: ready to attempt — UNVERIFIED on X-point.** The separatrix-based
label normalization is now implemented (edge radius from the LCFS, open-region
nodes excluded, q matching capped below the separatrix, the `xpoint2` guard
relaxed). This is the first diverted validation case and has NOT yet been run;
the key unknown is whether a clean X-point forms from a pure-RE current (check
`ES%ifail_xpoint` / `LCFS_is_lost` in the equilibrium log first). Free-boundary
equilibria remain unsupported.

## What it is

JET-like lower single-null (R0 = 3 m, a = 1 m, B0 ~ 3.3 T) pure-RE
fixed-boundary equilibrium: a mono-energetic 10 MeV RE beam whose current
lies on drift (constant canonical toroidal momentum) surfaces, with the
common profile function iterated until q matches `qprofile_target.dat`. The
geometry (elongation, lower squareness, and the X-point block `xpoint=.t.`,
`xampl`, and the diverted grid `n_open`/`n_leg`/`n_private`) follows
`fixedbound_equil_xpoint`; the pure-RE profile block (`FF_0=FF_1=0` with
finite shape widths, negligible T/rho) follows `re_kin_equil_600`.

## Why these choices for a FIRST diverted test

* **`re_eq_l_beam = 0.9`**: the beam current is confined well inside the
  separatrix, leaving a current-free annulus to the X-point/wall. This
  exercises the separatrix-based label normalization (outboard edge and
  `psi_bnd` taken from the LCFS, open-region nodes excluded) **without** the
  harder physics of current on drift orbits that cross the separatrix. Raise
  `l_beam` towards 1 only after the confined case is validated.
* **Modest energy (10 MeV) at reactor-scale field**: drift parameter
  `d_s ~ 1e-2`, drift shift ~1.5 cm — large enough to be meaningful, small
  enough that the beam edge stays clear of the separatrix.
* **`re_eq_map_mode = 'contour'`**: diverted plasmas are strongly shaped, so
  the geometry-agnostic nodal contour label map is the appropriate default.
* **`re_eq_tol_q = 1e-2`, `alpha_out = 0.2`**: the shaping lesson — ~1% is a
  physically adequate q-match and a milder outer under-relaxation avoids the
  limit cycle seen in strongly shaped cases.

## Physics checks (once it runs)

* `re_eq_convergence.log` converges (hard or soft) with `edge_fraction` small
  (beam confined inside the separatrix).
* The drift axis and beam sit inside the LCFS; per-class beam-edge `psihat`
  reported at iteration 1 should be < 0.98 (clear of the separatrix).
* `qprofile.dat` agrees with `qprofile_target.dat` below the beam edge; q
  between the beam edge and the separatrix is an outcome (steepening toward
  q -> infinity at the separatrix), not a target.
* Coupled `kinetic_main` phase: the projected RE current and q stay
  stationary at marker-noise level — the same stationarity payoff verified
  for circular and shaped limiter cases, now with an X-point.

See `models/mod_re_kinetic_equilibrium.f90` and
`util/re_equilibrium_prototype/doc/re_kinetic_equilibrium.tex` (esp. the
limitations/outlook section on X-point support).
