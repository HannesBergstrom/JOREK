# Kinetic RE drift-surface equilibrium prototype (M0)

Standalone Python reference implementation of the kinetic runaway-electron
equilibrium solver: given a discrete RE momentum/pitch distribution and a
target safety-factor profile `q_t(psihat_n)`, it computes the axisymmetric
fixed-boundary equilibrium in which the RE current lies on drift surfaces
(constant canonical toroidal momentum `A_s = gamma_s m_e v_par,s R - e psi`)
rather than flux surfaces.

Physics basis: V. Bandaru and M. Hoelzl, Phys. Plasmas 30, 092508 (2023),
Formulation B [Eqs. (3.23)-(3.27), Sec. IV], generalized from mono-energetic
to a set of momentum/pitch classes sharing one common profile function
`Nprof(Ahat)`, plus an outer profile-transplant loop matching `q_t`.

This tool validates the physics and the iteration scheme for the in-code
implementation (`models/mod_re_kinetic_equilibrium.f90`) and stays in the
repository as a documented reference. Everything is in SI units; physical
constants match `models/constants.f90`.

## Usage

```sh
# q-matched multi-class equilibrium (ITER-like circular, defaults R0=6.2 a=2 B0=5.3):
python3 re_eq_prototype.py --dist dist_spectrum_8nodes.dat \
                           --q-target qprofile_target.dat

# match only the q shape at prescribed RE current:
python3 re_eq_prototype.py --dist dist_mono_40MeV.dat \
                           --q-target qprofile_target.dat \
                           --match-mode q_shape --I-RE 10e6

# full validation suite (low-energy limit, Bandaru-like shifts, q-matching):
python3 run_validation.py
```

Outputs: `<output>_psi.npz` (psi map), `<output>_nprof.dat` (common profile
table), `<output>_classes.dat` (per-class gamma, v_par, weights, drift-axis
data, lost-current fractions), `<output>_convergence.log`.

## Input formats

* Distribution table (`--dist`, format `ekin_xi_w`): one class per line,
  columns `E_kin[eV]  xi  weight`, `#` comments. Readers are registered in
  `_DIST_READERS` in `re_eq_prototype.py`; adding a new table format means
  adding one reader function there (the in-code implementation mirrors this).
* Target q profile (`--q-target`): two columns `psihat_n  q`.

## Validation status

`run_validation.py`: 11/11 pass (2026-07-15).

1. **Low-energy limit (100 keV)**: drift shift 0.1 mm (negligible); psi field
   identical to a flux-surface-parametrized solve to 3e-4; q-matching
   converges monotonically (pins the transplant exponent sign and all sign
   conventions).
2. **Mono-energetic ITER-like benchmark (10 MA, 40/80 MeV)**: boundary
   drift-surface shifts 11.1 / 21.5 cm, within 15% of the analytic passing
   orbit estimate `2 q_edge gamma m v_par / (e B(R_edge))` and scaling
   linearly with momentum (cf. Bandaru & Hoelzl Fig. 1; the quantitative
   5%-tolerance figure comparison requires their exact current profile and is
   part of the in-code test suite).
3. **q-matching**: mono 20 MeV, mono 35 MeV, and an 8-node two-decade
   exponential spectrum all converge in both `full_q` and `q_shape` modes
   (alpha_out = 0.3): the 20 MeV and spectrum cases to
   `max|q/q_t - 1| < 1e-3` in **10 outer iterations**, the 35 MeV case to
   2.0e-3 in 27. The convergence metric covers the FULL evaluated psihat
   range; it was previously restricted to the beam-edge label range, which
   understated the error by ~5x once the drift smearing is wide (the former
   60 MeV case read 2.4e-3 but is really at 1.2e-2). That floor is set by the
   psihat span of a drift surface: the map Nprof(Ahat) -> I(psihat) smooths
   over ~|alpha| (R_out - R_in) / |dpsi|, so target structure finer than that
   cannot be represented, worst at the edge and growing with class energy.

## M0 findings that the in-code implementation MUST carry over

These were found the hard way; each one produced a stall or a limit cycle of
the outer iteration when absent:

1. **Cumulative transplant (default).** The literal pointwise update
   `Nprof *= (q/q_t)^alpha` converges but with a very slow tail (~150
   iterations) whenever the residual sits near the edge: the local Nprof
   there carries almost no current, while q responds to the *enclosed*
   current (lower-triangular response). Transplanting the cumulative profile
   `C(l) = int_0^l Nprof dl'` instead — `C *= ratio^alpha`, `Nprof = dC/dl`,
   clipped at 0 — is the cylindrical-exact Newton direction and converges in
   <= 10 iterations. Both variants share the same fixed point; the pointwise
   variant is kept selectable (`transplant='pointwise'`).
2. **Clamp the ratio** to `[1/2, 2]` per outer iteration before applying the
   exponent: a poor initial profile otherwise drives updates over several
   decades and the nested iteration limit-cycles.
3. **Evaluate `q_now` and `q_t` at the same clamped argument** where the
   label map leaves the computed q range (near axis/edge); otherwise the
   comparison has a bias floor `~ dq_t/dpsihat * clamp distance` that the
   iteration can never remove (observed: exactly 1.6e-3 with the default
   target profile).
4. **Measure convergence directly in q space** (`max|q(psihat)/q_t(psihat)-1|`
   over the q-diagnostic levels, restricted to the label-controllable range).
   Measuring through the label map adds interpolation noise of the map to the
   metric; the update cannot and need not remove it. Smooth interpolation of
   q(psihat) (monotone cubic) and a mild smoothing of the log-ratio before
   the transplant keep the update from chasing single-point q-evaluation
   artifacts. With JOREK's cubic finite elements the q diagnostic is much
   more accurate than this prototype's bilinear one, but the structural rules
   stand.
5. **Initial Nprof from q_t (cylindrical)**: use the analytic derivative
   `j = (B0/(mu0 R0)) (2/q - r q'/q^2)`; a finite-difference `I'(r)/r`
   blows up at the axis and poisons the first iterations.
6. **Per-class drift axes** must be relocated every Picard iteration (they
   sit 1-15 cm outboard of the psi axis for 0.6-80 MeV classes and move as
   psi converges).

Sign conventions validated here (electron charge q = -e, `A_s/e = alpha_s R -
psi` with `alpha_s = gamma m_e v_par / e`, `j_phi = -e sum_s v_par,s n_s`,
`Delta* psi = mu0 e sum_s v_par,s w_s Nprof(Ahat_s)`, psi maximum at axis for
j_phi > 0): the low-energy test reproduces the standard flux-surface
equilibrium and q > 0 with the correct monotone transplant response.

Known prototype-only simplifications: circular fixed boundary, small-pitch
form of A_s (the finite-pitch `R B_phi/B` variant is an in-code namelist
option), bilinear q diagnostic, no thermal pressure.
