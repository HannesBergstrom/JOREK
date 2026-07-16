#!/usr/bin/env python3
"""
M0 validation suite for the kinetic RE drift-surface equilibrium prototype.

Test 1 (low-energy limit): at 100 keV the drift shift must be negligible,
        the solution must coincide with a flux-surface-parametrized solve,
        and the q-transplant update must converge monotonically (this pins
        the exponent sign and all sign conventions).
Test 2 (Bandaru benchmark):  ITER-like circular plasma R0=6.2 m, a=2 m,
        B0=5.3 T, I_RE=10 MA, mono-energetic 40/80 MeV: drift-surface
        shifts of order ~8 cm / ~15 cm and scaling ~ gamma between them
        (Bandaru & Hoelzl 2023, Fig. 1).
Test 3 (q matching): mono 20 MeV, mono 60 MeV, and an 8-node two-decade
        spectrum must converge to max|q/q_t - 1| < 1e-3 in < 30 outer
        iterations, in both full_q and q_shape modes.

Run:  python3 run_validation.py
"""

import sys
import time
import numpy as np

from re_eq_prototype import (REClasses, PolarGSSolver, REEquilibrium,
                             Nprof, QTarget, read_q_target)

R0, A_MIN, B0 = 6.2, 2.0, 5.3
F0 = B0 * R0

PASS, FAIL = "PASS", "FAIL"
results = []


def check(name, ok, detail):
    results.append((name, ok, detail))
    print(f"[{PASS if ok else FAIL}] {name}: {detail}")


def q_target_default():
    return read_q_target('qprofile_target.dat')


# ---------------------------------------------------------------------------
# Test 1: low-energy limit
# ---------------------------------------------------------------------------
def test_low_energy():
    print("\n=== Test 1: low-energy limit (100 keV) ===")
    cl = REClasses(E_kin=[1.0e5], xi=[-0.99], w=[1.0])
    gs = PolarGSSolver(R0, A_MIN, Nr=96, Nt=192)
    nprof = Nprof((1.0 - np.linspace(0, 1, 101)**1.5)**2 * 1.0e16)

    eq = REEquilibrium(cl, gs, F0, match_mode='q_shape', I_RE=10.0e6,
                       verbose=False)
    n_in, res = eq.solve_fixed_nprof(Nprof(nprof.N.copy(), nprof.l))
    shifts = eq.drift_shifts()[0]

    # 1a: drift shift negligible (analytic estimate ~ q * gamma m v/(e B) ~ 0.4 mm)
    check("low-E axis shift negligible", shifts['axis_shift'] < 5.0e-3,
          f"axis shift = {shifts['axis_shift'] * 1e3:.3f} mm (inner {n_in} it, res {res:.1e})")

    # 1b: identical to a flux-surface solve (alpha_s -> 0 exactly)
    cl0 = REClasses(E_kin=[1.0e5], xi=[-0.99], w=[1.0])
    cl0.alpha[:] = 0.0     # A_s = -psi: labels are exactly flux-surface labels
    eq0 = REEquilibrium(cl0, gs, F0, match_mode='q_shape', I_RE=10.0e6,
                        verbose=False)
    eq0.solve_fixed_nprof(Nprof(nprof.N.copy(), nprof.l))
    dpsi_rel = (np.abs(eq.psi - eq0.psi).max()
                / np.abs(eq0.psi - gs.psi_b).max())
    check("low-E flux-surface limit", dpsi_rel < 2.0e-3,
          f"max|psi - psi_fluxsurf|/|psi| = {dpsi_rel:.2e}")

    # 1c: q-transplant converges monotonically at low energy (sign check)
    qt = q_target_default()
    eqm = REEquilibrium(cl, gs, F0, q_target=qt, match_mode='full_q',
                        alpha_out=0.3, tol_q=1e-3, max_it_out=30,
                        verbose=False)
    ok = eqm.match_q()
    errs = [rec['q_err'] for rec in eqm.log]
    n_up = sum(1 for k in range(1, len(errs)) if errs[k] > errs[k - 1] * 1.05)
    check("low-E q-matching converges", ok and len(errs) <= 30,
          f"{len(errs)} outer iterations, final err = {errs[-1]:.2e}")
    check("low-E q-error decrease ~monotone", n_up == 0,
          f"error increased in {n_up}/{len(errs) - 1} steps "
          f"(exponent sign of the transplant update is correct)" if n_up == 0
          else f"error increased in {n_up} steps -- check transplant sign!")
    return eqm


# ---------------------------------------------------------------------------
# Test 2: mono-energetic Bandaru benchmark
# ---------------------------------------------------------------------------
def test_bandaru():
    """Drift-surface shifts of the ITER-like mono-energetic reference case.

    The exact comparison with Bandaru & Hoelzl Fig. 1 (5% tolerance) needs
    their exact current profile and belongs to the code-level test (M1/M2).
    Here the measured boundary shift is checked against the analytic
    passing-orbit estimate for the in-out asymmetry of a drift surface,
        Delta_bnd ~ 2 q_edge gamma m v_par / (e B_phi(R_edge)),
    (B_phi taken at the outboard edge, where the shift is measured) and
    against the literature ballpark of several-to-tens of cm."""
    print("\n=== Test 2: Bandaru benchmark (10 MA, 40 / 80 MeV) ===")
    from re_eq_prototype import MASS_ELECTRON, EL_CHG
    nprof0 = Nprof((1.0 - np.linspace(0, 1, 101)**1.5)**2 * 1.0e16)
    out = {}
    for E_MeV in (40.0, 80.0):
        cl = REClasses(E_kin=[E_MeV * 1e6], xi=[-0.99], w=[1.0])
        gs = PolarGSSolver(R0, A_MIN, Nr=128, Nt=256)
        eq = REEquilibrium(cl, gs, F0, match_mode='q_shape', I_RE=10.0e6,
                           verbose=False)
        t0 = time.time()
        n_in, res = eq.solve_fixed_nprof(Nprof(nprof0.N.copy(), nprof0.l))
        sh = eq.drift_shifts()[0]
        ph, qn = eq.q_profile()
        q_edge = qn[-1]
        B_edge = F0 / (R0 + A_MIN)
        rho_par = (cl.gamma[0] * MASS_ELECTRON * abs(cl.v_par[0])
                   / (EL_CHG * B_edge))
        sh['delta_est'] = 2.0 * q_edge * rho_par
        out[E_MeV] = sh
        print(f"  {E_MeV:5.1f} MeV: d_s = {sh['d_s']:+.4e}, q_edge = {q_edge:.3f}, "
              f"axis shift = {sh['axis_shift'] * 100:.2f} cm, "
              f"boundary shift = {sh['boundary_shift'] * 100:.2f} cm "
              f"(analytic est {sh['delta_est'] * 100:.2f} cm), "
              f"lost = {sh['lost_fraction']:.1e} "
              f"({n_in} it, res {res:.1e}, {time.time() - t0:.1f} s)")

    for E_MeV, lo, hi in ((40.0, 4.0, 16.0), (80.0, 8.0, 30.0)):
        b = out[E_MeV]['boundary_shift'] * 100.0
        est = out[E_MeV]['delta_est'] * 100.0
        ok = (abs(b / est - 1.0) < 0.2) and (lo <= b <= hi)
        check(f"boundary shift vs analytic estimate ({E_MeV:.0f} MeV)", ok,
              f"measured {b:.2f} cm vs estimate {est:.2f} cm "
              f"(ballpark [{lo:.0f}, {hi:.0f}] cm)")
    b40 = out[40.0]['boundary_shift']
    b80 = out[80.0]['boundary_shift']
    p_ratio = ((out[80.0]['d_s']) / (out[40.0]['d_s']))
    check("shift scaling ~ p (gamma ratio)",
          abs(b80 / b40 / p_ratio - 1.0) < 0.15,
          f"shift(80)/shift(40) = {b80 / b40:.2f} (p ratio {p_ratio:.2f})")


# ---------------------------------------------------------------------------
# Test 3: q matching
# ---------------------------------------------------------------------------
def test_q_matching():
    print("\n=== Test 3: q-profile matching ===")
    qt = q_target_default()

    cases = [
        ("mono 20 MeV", REClasses(E_kin=[2.0e7], xi=[-0.99], w=[1.0])),
        ("mono 60 MeV", REClasses(E_kin=[6.0e7], xi=[-0.99], w=[1.0])),
        ("8-node spectrum", REClasses.from_file('dist_spectrum_8nodes.dat')),
    ]
    for name, cl in cases:
        gs = PolarGSSolver(R0, A_MIN, Nr=96, Nt=192)
        eq = REEquilibrium(cl, gs, F0, q_target=qt, match_mode='full_q',
                           alpha_out=0.3, tol_q=1e-3, max_it_out=30,
                           verbose=False)
        t0 = time.time()
        ok = eq.match_q()
        n_out = len(eq.log)
        err = eq.log[-1]['q_err']
        I_MA = eq.log[-1]['I_RE'] / 1e6
        check(f"q-match ({name}, full_q)", ok and n_out < 30,
              f"{n_out} outer it, err = {err:.2e}, I_RE = {I_MA:.3f} MA, "
              f"{time.time() - t0:.1f} s")

    # q_shape mode: fix I_RE, match the shape, report amplitude
    cl = REClasses(E_kin=[2.0e7], xi=[-0.99], w=[1.0])
    gs = PolarGSSolver(R0, A_MIN, Nr=96, Nt=192)
    eq = REEquilibrium(cl, gs, F0, q_target=qt, match_mode='q_shape',
                       I_RE=8.0e6, alpha_out=0.3, tol_q=1e-3, max_it_out=30,
                       verbose=False)
    ok = eq.match_q()
    rec = eq.log[-1]
    check("q-match (mono 20 MeV, q_shape)", ok,
          f"{len(eq.log)} outer it, err = {rec['q_err']:.2e}, "
          f"I_RE = {rec['I_RE'] / 1e6:.3f} MA (prescribed 8), "
          f"q amplitude = {rec['q_amplitude']:.3f}")


# ---------------------------------------------------------------------------
if __name__ == '__main__':
    t_start = time.time()
    test_low_energy()
    test_bandaru()
    test_q_matching()

    print(f"\n=== Summary ({time.time() - t_start:.0f} s) ===")
    n_fail = sum(1 for _, ok, _ in results if not ok)
    for name, ok, _ in results:
        print(f"  [{PASS if ok else FAIL}] {name}")
    print(f"{len(results) - n_fail}/{len(results)} passed")
    sys.exit(1 if n_fail else 0)
