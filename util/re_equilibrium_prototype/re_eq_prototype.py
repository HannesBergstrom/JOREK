#!/usr/bin/env python3
"""
Standalone prototype of the kinetic runaway-electron (RE) drift-surface
equilibrium solver (milestone M0 of the kinetic RE equilibrium project).

Physics
-------
Given a discrete set of RE momentum/pitch classes s = 1..N_s (read from a
plain-text table) and a target safety-factor profile q_t(psihat_n), compute
the axisymmetric fixed-boundary equilibrium in which the RE current lies on
drift surfaces (surfaces of constant canonical toroidal momentum) rather
than flux surfaces:

    A_s(R,Z)   = gamma_s m_e v_par,s R - e psi        (small-pitch form)
    Ahat_s     = (A_s - A_s,axis) / (A_s,edge - A_s,axis)
    n_s(R,Z)   = w_s Nprof(Ahat_s) / R
    j_phi      = -e sum_s v_par,s n_s
    Delta* psi = -mu0 R j_phi = mu0 e sum_s v_par,s w_s Nprof(Ahat_s)

with one common normalized profile function Nprof shared by all classes,
determined by an outer profile-transplant iteration that matches q(psihat_n)
to the target q_t.  See Bandaru & Hoelzl, Phys. Plasmas 30, 092508 (2023),
Formulation B [Eqs. (3.23)-(3.27) and Sec. IV], generalized from
mono-energetic to multiple classes.

Numerics
--------
Fixed-boundary circular plasma. Grad-Shafranov operator discretized with
second-order finite differences on a polar grid (r, theta) centred on the
geometric axis (R0, 0); r-nodes at half-integer positions avoid the
coordinate singularity at r = 0. The linear operator is factorized once
(sparse LU) and reused in all Picard iterations.

Everything is in SI units. Physical constants match models/constants.f90.

This tool is a documented physics reference for the in-code implementation;
it deliberately mirrors the structure planned for
models/mod_re_kinetic_equilibrium.f90.
"""

import numpy as np
from scipy import sparse
from scipy.sparse import linalg as sla

# --- Physical constants (values as in models/constants.f90) ---------------
MU_ZERO = 4.0e-7 * np.pi        # [Vs/Am]
EL_CHG = 1.602176565e-19        # [C]
MASS_ELECTRON = 9.10938291e-31  # [kg]
C_LIGHT = 299792458.0           # [m/s]


# ===========================================================================
# Input: RE momentum-space distribution
# ===========================================================================

def _read_table_ekin_xi_w(fname):
    """Reader for the 'ekin_xi_w' plain-text format: one class per line,
    columns E_kin [eV], pitch xi = p_par/p, relative weight. '#' comments."""
    rows = []
    with open(fname) as f:
        for line in f:
            line = line.split('#')[0].strip()
            if not line:
                continue
            cols = [float(c) for c in line.split()]
            if len(cols) != 3:
                raise ValueError(
                    f"{fname}: expected 3 columns (E_kin[eV] xi weight), got: {line}")
            rows.append(cols)
    if not rows:
        raise ValueError(f"{fname}: no distribution classes found")
    return np.array(rows)


# Registry of table formats: adding a new format = adding one reader here.
_DIST_READERS = {
    'ekin_xi_w': _read_table_ekin_xi_w,
}


class REClasses:
    """The discrete momentum/pitch classes of the RE distribution.

    Attributes (arrays of length n_s):
      E_kin  : kinetic energy [eV]
      xi     : pitch p_par/p (sign relative to the toroidal field direction)
      w      : weights, normalized to sum(w) = 1
      gamma  : Lorentz factor
      v      : speed [m/s]
      v_par  : parallel velocity [m/s]
      alpha  : gamma m_e v_par / e  [Wb/m], so that A_s/e = alpha_s R - psi
    """

    def __init__(self, E_kin, xi, w, xi_min=0.9):
        E_kin = np.atleast_1d(np.asarray(E_kin, dtype=float))
        xi = np.atleast_1d(np.asarray(xi, dtype=float))
        w = np.atleast_1d(np.asarray(w, dtype=float))
        if np.any(np.abs(xi) < xi_min):
            raise ValueError(
                f"Distribution contains classes with |xi| < xi_min = {xi_min}: "
                "trapped-particle invariants are out of scope for this solver.")
        if np.any(w < 0.0):
            raise ValueError("Distribution weights must be non-negative")
        self.E_kin = E_kin
        self.xi = xi
        self.w = w / w.sum()
        mec2_eV = MASS_ELECTRON * C_LIGHT**2 / EL_CHG
        self.gamma = 1.0 + E_kin / mec2_eV
        self.v = C_LIGHT * np.sqrt(1.0 - 1.0 / self.gamma**2)
        self.v_par = self.xi * self.v
        # A_s / e = alpha_s * R - psi
        self.alpha = self.gamma * MASS_ELECTRON * self.v_par / EL_CHG
        self.n_s = len(E_kin)

    @classmethod
    def from_file(cls, fname, fmt='ekin_xi_w', xi_min=0.9):
        rows = _DIST_READERS[fmt](fname)
        return cls(rows[:, 0], rows[:, 1], rows[:, 2], xi_min=xi_min)

    def alpha_eff(self):
        """Current-density-weighted effective alpha, used by the label map."""
        cw = self.w * self.v_par
        return np.sum(cw * self.alpha) / np.sum(cw)

    def d_shift(self, B0, a):
        """Per-class drift-shift parameter d_s = gamma m v_par/(e B0 a):
        the single most useful sanity number (Delta_s ~ d_s q / (r/a))."""
        return self.alpha / (B0 * a)


# ===========================================================================
# Target q profile
# ===========================================================================

def read_q_target(fname):
    """Plain-text target q profile: two columns psihat_n, q_t."""
    dat = np.loadtxt(fname)
    if dat.ndim != 2 or dat.shape[1] != 2:
        raise ValueError(f"{fname}: expected two columns psihat_n, q_t")
    return QTarget(dat[:, 0], dat[:, 1])


class QTarget:
    """Target q profile. Evaluated with MONOTONE-CUBIC interpolation: a
    piecewise-linear target has kinks that no smooth equilibrium q can
    match, which floors the achievable q error at
    ~ curvature * table-spacing^2 (observed: exactly 1e-3 with an 11-point
    table -- right at the default tolerance)."""

    def __init__(self, psihat, q):
        from scipy.interpolate import PchipInterpolator
        idx = np.argsort(psihat)
        self.psihat = np.asarray(psihat, dtype=float)[idx]
        self.q = np.asarray(q, dtype=float)[idx]
        self._interp = PchipInterpolator(self.psihat, self.q)

    def __call__(self, psihat):
        return self._interp(np.clip(psihat, self.psihat[0], self.psihat[-1]))


# ===========================================================================
# Grad-Shafranov solver on a polar grid (fixed circular boundary)
# ===========================================================================

class PolarGSSolver:
    """Fixed-boundary GS solver: Delta* psi = S(psi; R, Z) with psi = psi_b
    on the circle of radius a around (R0, 0).

    Polar grid: r_i = (i + 1/2) dr, i = 0..Nr-1 ; theta_j = j dtheta,
    j = 0..Nt-1 (periodic; Nt must be even for the across-centre stencil).
    """

    def __init__(self, R0, a, Nr=128, Nt=256, psi_b=0.0):
        assert Nt % 2 == 0, "Nt must be even (across-centre stencil at r=0)"
        self.R0, self.a = R0, a
        self.Nr, self.Nt = Nr, Nt
        self.psi_b = psi_b
        self.dr = a / Nr
        self.dth = 2.0 * np.pi / Nt
        self.r = (np.arange(Nr) + 0.5) * self.dr
        self.th = np.arange(Nt) * self.dth
        self.RR = self.R0 + np.outer(self.r, np.cos(self.th))   # (Nr, Nt)
        self.ZZ = np.outer(self.r, np.sin(self.th))
        self._factorize()

    # --- operator assembly --------------------------------------------------
    def _factorize(self):
        """Assemble Delta* in polar coordinates and LU-factorize it.

        Delta* psi = psi_rr + psi_r/r + psi_tt/r^2
                     - (1/R) (cos th * psi_r - sin th * psi_th / r)
        """
        Nr, Nt, dr, dth = self.Nr, self.Nt, self.dr, self.dth
        n = Nr * Nt

        def idx(i, j):
            return i * Nt + (j % Nt)

        rows, cols, vals = [], [], []
        self._bnd_rhs = np.zeros(n)  # contribution of the Dirichlet boundary

        for i in range(Nr):
            r = self.r[i]
            for j in range(Nt):
                th = self.th[j]
                R = self.R0 + r * np.cos(th)
                k = idx(i, j)

                # coefficients of psi_{i+1}, psi_{i-1}, psi_{i}, psi_{j+1}, psi_{j-1}
                c_rr = 1.0 / dr**2
                c_r = (1.0 / r - np.cos(th) / R) / (2.0 * dr)
                c_tt = 1.0 / (r * dth)**2
                c_t = (np.sin(th) / (R * r)) / (2.0 * dth)

                # centre
                rows.append(k); cols.append(k); vals.append(-2.0 * c_rr - 2.0 * c_tt)
                # theta neighbours (periodic)
                rows.append(k); cols.append(idx(i, j + 1)); vals.append(c_tt + c_t)
                rows.append(k); cols.append(idx(i, j - 1)); vals.append(c_tt - c_t)
                # r + dr
                cp = c_rr + c_r
                if i + 1 < Nr:
                    rows.append(k); cols.append(idx(i + 1, j)); vals.append(cp)
                else:
                    # ghost at r = a + dr/2: psi_ghost = 2 psi_b - psi_{Nr-1}
                    rows.append(k); cols.append(k); vals.append(-cp)
                    self._bnd_rhs[k] += -2.0 * self.psi_b * cp
                # r - dr
                cm = c_rr - c_r
                if i - 1 >= 0:
                    rows.append(k); cols.append(idx(i - 1, j)); vals.append(cm)
                else:
                    # across the centre: (r=-dr/2, th) == (dr/2, th+pi)
                    rows.append(k); cols.append(idx(0, j + Nt // 2)); vals.append(cm)

        A = sparse.csr_matrix((vals, (rows, cols)), shape=(n, n))
        self._lu = sla.splu(A.tocsc())

    def solve(self, S):
        """Solve Delta* psi = S for the given source field S (Nr, Nt)."""
        rhs = S.reshape(-1) + self._bnd_rhs
        return self._lu.solve(rhs).reshape(self.Nr, self.Nt)

    # --- interpolation helpers ---------------------------------------------
    def interp(self, psi, R, Z):
        """Bilinear interpolation of a polar-grid field at points (R, Z)."""
        R = np.asarray(R, dtype=float)
        Z = np.asarray(Z, dtype=float)
        r = np.sqrt((R - self.R0)**2 + Z**2)
        th = np.mod(np.arctan2(Z, R - self.R0), 2.0 * np.pi)

        # radial index into the half-integer grid; clamp to [0, Nr-2]
        fi = np.clip(r / self.dr - 0.5, 0.0, self.Nr - 1.0 - 1e-12)
        i0 = np.minimum(fi.astype(int), self.Nr - 2)
        wi = fi - i0
        fj = th / self.dth
        j0 = fj.astype(int) % self.Nt
        wj = fj - fj.astype(int)
        j1 = (j0 + 1) % self.Nt

        f = psi
        return ((1 - wi) * (1 - wj) * f[i0, j0] + (1 - wi) * wj * f[i0, j1]
                + wi * (1 - wj) * f[i0 + 1, j0] + wi * wj * f[i0 + 1, j1])

    def grad(self, psi):
        """(dpsi/dR, dpsi/dZ) on the polar grid (second-order differences)."""
        dpsi_dr = np.empty_like(psi)
        dpsi_dr[1:-1, :] = (psi[2:, :] - psi[:-2, :]) / (2.0 * self.dr)
        # across-centre at i=0, one-sided ghost at the boundary
        psi_across = np.roll(psi[0, :], self.Nt // 2)
        dpsi_dr[0, :] = (psi[1, :] - psi_across) / (2.0 * self.dr)
        psi_ghost = 2.0 * self.psi_b - psi[-1, :]
        dpsi_dr[-1, :] = (psi_ghost - psi[-2, :]) / (2.0 * self.dr)

        dpsi_dt = (np.roll(psi, -1, axis=1) - np.roll(psi, 1, axis=1)) / (2.0 * self.dth)

        ct, st = np.cos(self.th), np.sin(self.th)
        rinv = 1.0 / self.r[:, None]
        dpsi_dR = ct[None, :] * dpsi_dr - st[None, :] * rinv * dpsi_dt
        dpsi_dZ = st[None, :] * dpsi_dr + ct[None, :] * rinv * dpsi_dt
        return dpsi_dR, dpsi_dZ

    def find_extremum(self, field, kind):
        """Locate the interior extremum of a polar-grid field.

        kind = 'min' or 'max'. Grid search plus local paraboloid refinement
        (this mirrors what the in-code axis finder must do per class).
        Returns (R_ax, Z_ax, value).
        """
        f = field if kind == 'max' else -field
        i, j = np.unravel_index(np.argmax(f), f.shape)
        R_c, Z_c = self.RR[i, j], self.ZZ[i, j]

        # paraboloid fit on neighbouring points in (R, Z)
        ii = np.clip(np.arange(i - 2, i + 3), 0, self.Nr - 1)
        jj = np.arange(j - 2, j + 3) % self.Nt
        I, J = np.meshgrid(ii, jj, indexing='ij')
        x = self.RR[I, J].ravel() - R_c
        y = self.ZZ[I, J].ravel() - Z_c
        z = field[I, J].ravel()
        M = np.column_stack([np.ones_like(x), x, y, x * x, x * y, y * y])
        c, *_ = np.linalg.lstsq(M, z, rcond=None)
        H = np.array([[2 * c[3], c[4]], [c[4], 2 * c[5]]])
        try:
            dx, dy = np.linalg.solve(H, [-c[1], -c[2]])
        except np.linalg.LinAlgError:
            dx = dy = 0.0
        lim = 2.0 * self.dr
        dx, dy = np.clip(dx, -lim, lim), np.clip(dy, -lim, lim)
        val = (c[0] + c[1] * dx + c[2] * dy + c[3] * dx * dx
               + c[4] * dx * dy + c[5] * dy * dy)
        return R_c + dx, Z_c + dy, val


# ===========================================================================
# The common profile function Nprof
# ===========================================================================

class Nprof:
    """Common normalized profile function Nprof(l), l in [0,1], shared by all
    classes: tabulated on a uniform l grid, linear interpolation, clipped."""

    def __init__(self, values, l_grid=None):
        self.N = np.asarray(values, dtype=float).copy()
        self.l = (np.linspace(0.0, 1.0, len(self.N))
                  if l_grid is None else np.asarray(l_grid, dtype=float))

    def __call__(self, l):
        return np.interp(np.clip(l, 0.0, 1.0), self.l, self.N)

    def scale(self, factor):
        self.N *= factor


# ===========================================================================
# The equilibrium solver: inner Picard + outer q-profile transplant
# ===========================================================================

class REEquilibrium:
    """Multi-class RE drift-surface equilibrium with q-profile matching."""

    def __init__(self, classes, solver, F0, q_target=None,
                 match_mode='full_q', I_RE=None,
                 alpha_in=0.5, tol_in=1e-10, max_it_in=200,
                 alpha_out=0.3, tol_q=1e-3, max_it_out=50,
                 transplant='cumulative', edge_taper=0.2,
                 l_beam=1.0, l_beam_width=0.1,
                 n_l=101, n_theta_q=256, verbose=True):
        self.cl = classes
        self.gs = solver
        self.F0 = F0
        self.qt = q_target
        self.match_mode = match_mode
        self.I_RE = I_RE
        if match_mode == 'q_shape' and I_RE is None:
            raise ValueError("match_mode='q_shape' requires a prescribed I_RE")
        self.alpha_in, self.tol_in, self.max_it_in = alpha_in, tol_in, max_it_in
        self.alpha_out, self.tol_q, self.max_it_out = alpha_out, tol_q, max_it_out
        if transplant not in ('cumulative', 'pointwise'):
            raise ValueError(f"unknown transplant variant '{transplant}'")
        self.transplant = transplant
        # Truncation policy for drift surfaces leaving the domain (Ahat > 1):
        # their current is REMOVED (RE orbits crossing the wall are lost),
        # with a linear taper of width edge_taper in the label for numerical
        # smoothness. The clamp alternative (edge_taper=None: keep Nprof(1)
        # outside) creates an uncontrollable halo current that destabilizes
        # the outer iteration once the lost fraction is more than a few
        # percent (observed as stall at ~3e-3 followed by slow divergence).
        self.edge_taper = edge_taper
        # Beam-edge label: current confined to Ahat < l_beam with a
        # smoothstep roll-off of width l_beam_width, leaving a current-free
        # (vacuum) annulus to the wall. For < 1 all current-carrying drift
        # orbits are closed (no scrape-off). In the zero-drift-orbit limit
        # l_beam equals the normalized poloidal flux of the beam edge.
        self.l_beam = l_beam
        self.l_beam_width = l_beam_width
        self.n_theta_q = n_theta_q
        self.verbose = verbose

        self.nprof = None
        self.psi = np.zeros((solver.Nr, solver.Nt))
        # per-class label state, refreshed every Picard iteration
        self.A_ax = np.zeros(self.cl.n_s)     # A_s/e at the class drift axis [Wb]
        self.A_edge = np.zeros(self.cl.n_s)
        self.ax_RZ = np.zeros((self.cl.n_s, 2))
        self.lost_fraction = np.zeros(self.cl.n_s)
        self.log = []                          # convergence log records

    # --- per-class invariant labels ------------------------------------------
    def _A_field(self, s, psi):
        """A_s/e = alpha_s R - psi on the grid [Wb]."""
        return self.cl.alpha[s] * self.gs.RR - psi

    def _extremum_kind(self, psi):
        """The extremum type of A_s at its drift axis is opposite to psi's
        (A ~ -psi near the axis, the alpha R term only shifts the extremum)."""
        psi_ax = psi[0, :].mean()
        return 'min' if psi_ax > self.gs.psi_b else 'max'

    def update_labels(self, psi):
        """Per-class drift axes and edge values -> normalization of Ahat_s.
        Must be called every time psi changes (the axes move!)."""
        kind = self._extremum_kind(psi)
        R_edge = self.gs.R0 + self.gs.a       # outboard midplane boundary
        for s in range(self.cl.n_s):
            A = self._A_field(s, psi)
            R_ax, Z_ax, A_ax = self.gs.find_extremum(A, kind)
            self.ax_RZ[s] = (R_ax, Z_ax)
            self.A_ax[s] = A_ax
            self.A_edge[s] = self.cl.alpha[s] * R_edge - self.gs.psi_b
            if self.A_edge[s] == self.A_ax[s]:
                raise RuntimeError(f"class {s}: degenerate label normalization")

    def Ahat(self, s, psi, R=None, Z=None, clip=True):
        """Normalized invariant label of class s (on the grid, or at points)."""
        if R is None:
            A = self._A_field(s, psi)
        else:
            A = self.cl.alpha[s] * np.asarray(R) - self.gs.interp(psi, R, Z)
        l = (A - self.A_ax[s]) / (self.A_edge[s] - self.A_ax[s])
        return np.clip(l, 0.0, 1.0) if clip else l

    # --- GS source ------------------------------------------------------------
    def _taper(self, l_raw):
        """Edge truncation factor: 1 on closed drift surfaces (l <= 1),
        C1 smoothstep decay to 0 over edge_taper beyond, 0 outside (a
        linear taper has slope kinks that the discretization and the
        transplant iteration ring against)."""
        if self.edge_taper is None:
            return np.ones_like(l_raw)
        t = np.clip((np.maximum(l_raw, 1.0) - 1.0) / self.edge_taper, 0.0, 1.0)
        return 1.0 - t*t*(3.0 - 2.0*t)

    def source(self, psi):
        """RHS of Delta* psi = mu0 e sum_s v_par,s w_s Nprof(Ahat_s).
        Note: no explicit R factor (the 1/R of n_s cancels the R of the GS RHS).
        Also accumulates the per-class lost-current fraction (drift surfaces
        with Ahat > 1, whose current is removed by the taper policy)."""
        S = np.zeros((self.gs.Nr, self.gs.Nt))
        for s in range(self.cl.n_s):
            l_raw = self.Ahat(s, psi, clip=False)
            N = self.nprof(l_raw) * self._taper(l_raw)
            S += MU_ZERO * EL_CHG * self.cl.v_par[s] * self.cl.w[s] * N
            with np.errstate(invalid='ignore'):
                N_unt = self.nprof(l_raw)
                tot = np.abs(N_unt).sum()
                lost = np.abs(N_unt[l_raw > 1.0]).sum()
            self.lost_fraction[s] = lost / tot if tot > 0 else 0.0
        return S

    def total_current(self, psi):
        """I_RE = integral of j_phi over the poloidal cross section [A]."""
        S = self.source(psi)                     # = -mu0 R j_phi ... /R absorbed
        # j_phi = -S/(mu0 R); dA = r dr dth
        j = -S / (MU_ZERO * self.gs.RR)
        dA = self.gs.r[:, None] * self.gs.dr * self.gs.dth
        return np.sum(j * dA)

    # --- inner loop: Picard on psi ---------------------------------------------
    def picard(self):
        """Solve GS with Ahat frozen from the previous psi, refresh labels,
        under-relax; iterate to tol_in."""
        psi = self.psi
        for it in range(1, self.max_it_in + 1):
            if self.match_mode == 'q_shape':
                self._rescale_to_I_RE(psi)
            S = self.source(psi)
            psi_new = self.gs.solve(S)
            dpsi = psi_new - psi
            psi = psi + self.alpha_in * dpsi
            self.update_labels(psi)
            norm = np.abs(psi - self.gs.psi_b).max()
            res = np.abs(dpsi).max() / norm if norm > 0 else np.inf
            if res < self.tol_in:
                break
        self.psi = psi
        return it, res

    def _rescale_to_I_RE(self, psi):
        I_now = self.total_current(psi)
        if I_now != 0.0:
            self.nprof.scale(self.I_RE / I_now)

    # --- q profile from the solved psi -----------------------------------------
    def psi_axis(self):
        kind = 'max' if self.psi[0, :].mean() > self.gs.psi_b else 'min'
        return self.gs.find_extremum(self.psi, kind)

    def q_profile(self, psihat_levels=None):
        """q(psihat_n) by flux-surface contour integration:
        q = F0/(2 pi) * closed-int dl_p / (R |grad psi|).
        Surfaces are traced by 1D root finding along rays from the psi axis
        (valid for the nested surfaces of this fixed-boundary prototype)."""
        if psihat_levels is None:
            psihat_levels = np.linspace(0.02, 0.985, 80)
        R_ax, Z_ax, psi_ax = self.psi_axis()
        dpsi = self.gs.psi_b - psi_ax
        gR, gZ = self.gs.grad(self.psi)

        th = np.linspace(0.0, 2.0 * np.pi, self.n_theta_q, endpoint=False)
        ct, st = np.cos(th), np.sin(th)

        # maximum ray length to the circular boundary
        b = (R_ax - self.gs.R0) * ct + Z_ax * st
        cc = (R_ax - self.gs.R0)**2 + Z_ax**2 - self.gs.a**2
        t_max = -b + np.sqrt(b * b - cc)

        qs = np.empty(len(psihat_levels))
        for n, lev in enumerate(psihat_levels):
            target = psi_ax + lev * dpsi
            t = self._ray_roots(R_ax, Z_ax, ct, st, t_max, target)
            Rp, Zp = R_ax + t * ct, Z_ax + t * st
            gRl = self.gs.interp(gR, Rp, Zp)
            gZl = self.gs.interp(gZ, Rp, Zp)
            gradpsi = np.sqrt(gRl**2 + gZl**2)
            dRp = np.roll(Rp, -1) - Rp
            dZp = np.roll(Zp, -1) - Zp
            dl = np.sqrt(dRp**2 + dZp**2)
            # midpoint values of the integrand
            f = 1.0 / (Rp * gradpsi)
            fm = 0.5 * (f + np.roll(f, -1))
            qs[n] = abs(self.F0) / (2.0 * np.pi) * np.sum(fm * dl)
        return np.asarray(psihat_levels), qs

    def _ray_roots(self, R_ax, Z_ax, ct, st, t_max, target, n_bisect=40):
        """Bisection for psi = target along rays from the axis (vectorized)."""
        lo = np.zeros_like(ct)
        hi = t_max.copy()
        psi_lo = self.gs.interp(self.psi, R_ax + lo * ct, Z_ax + lo * st)
        sign = np.sign(self.gs.interp(self.psi, R_ax + hi * ct, Z_ax + hi * st)
                       - psi_lo)
        for _ in range(n_bisect):
            mid = 0.5 * (lo + hi)
            pm = self.gs.interp(self.psi, R_ax + mid * ct, Z_ax + mid * st)
            below = sign * (pm - target) < 0.0
            lo = np.where(below, mid, lo)
            hi = np.where(below, hi, mid)
        return 0.5 * (lo + hi)

    # --- label map Ahat <-> psihat (swappable modelling choice) ----------------
    def label_to_psihat(self, l_values, alpha=None):
        """Midplane-average label map:
            psihat_m(l) = 0.5 * [psihat(R_out(l)) + psihat(R_in(l))]
        where R_out/R_in are the outboard/inboard midplane radii of the
        drift surface with label l for a class with the given alpha
        (default: the current-density-weighted effective class). Handles the
        near-axis degeneracy where both radii sit on the same side of the
        magnetic axis."""
        alpha_e = self.cl.alpha_eff() if alpha is None else alpha
        R_ax, Z_ax, psi_ax = self.psi_axis()
        dpsi = self.gs.psi_b - psi_ax

        # effective-class labels normalized like the per-class ones
        A = lambda R: alpha_e * R - self.gs.interp(self.psi, R, np.zeros_like(R))
        kind = self._extremum_kind(self.psi)
        # effective drift axis on the midplane: extremum of A along Z=0
        Rg = np.linspace(self.gs.R0 - self.gs.a + 1e-6,
                         self.gs.R0 + self.gs.a - 1e-6, 4001)
        Ag = A(Rg)
        i_ax = np.argmin(Ag) if kind == 'min' else np.argmax(Ag)
        R_dax = Rg[i_ax]
        A_ax = Ag[i_ax]
        A_edge = alpha_e * (self.gs.R0 + self.gs.a) - self.gs.psi_b
        lhat = (Ag - A_ax) / (A_edge - A_ax)

        # small negative values can occur when the fitted psi_ax differs
        # from the interpolated midplane maximum at the grid-error level
        psihat_mid = np.clip((self.gs.interp(self.psi, Rg, np.zeros_like(Rg))
                              - psi_ax) / dpsi, 0.0, 1.0)

        out = np.empty(len(l_values))
        for k, l in enumerate(np.clip(l_values, 1e-9, 1.0)):
            # outboard branch: R >= R_dax
            mo = Rg >= R_dax
            R_out = np.interp(l, lhat[mo], Rg[mo])
            # inboard branch: R <= R_dax (labels increase towards the wall)
            mi = Rg <= R_dax
            li, Ri = lhat[mi][::-1], Rg[mi][::-1]     # make label increasing
            if l <= li[-1]:
                R_in = np.interp(l, li, Ri)
            else:
                # drift surface exits the inboard midplane inside the domain
                # only on the outboard side (near-axis degeneracy / strong
                # shift): fall back to the outboard value
                R_in = R_out
            ph_out = np.interp(R_out, Rg, psihat_mid)
            ph_in = np.interp(R_in, Rg, psihat_mid)
            out[k] = 0.5 * (ph_out + ph_in)
        return out

    # --- outer loop: q-profile matching by profile transplant -------------------
    # Per-iteration clamp of the transplant ratio: without it, a poor initial
    # profile can drive q/q_t over several decades and the nested iteration
    # limit-cycles (the "under-relax both loops" pitfall). A factor-2 cap per
    # outer iteration keeps the update monotone while barely slowing
    # convergence from a good initial guess.
    RATIO_CLAMP = 2.0

    def match_q(self):
        """Outer q-matching loop: multiplicative profile transplant with
        under-relaxation, in one of two variants (self.transplant):

        'pointwise' (the literal update of the design document):
            Nprof_new(l) = Nprof_old(l) * ratio(l)^alpha_out,
            ratio(l) = q_now(psihat_m(l)) / q_t(psihat_m(l)).
          M0 finding: converges, but with a very slow tail whenever the
          residual sits near the edge, because the local Nprof there carries
          little current and q responds mainly to the *enclosed* current
          (the response is lower-triangular, not diagonal). ~150 iterations
          for tol 1e-3 on ITER-like cases.

        'cumulative' (default): transplant the cumulative profile
            C(l) = int_0^l Nprof dl'   (enclosed-current proxy),
            C_new(l) = C_old(l) * ratio(l)^alpha_out,  Nprof_new = dC_new/dl
          clipped at zero. In cylindrical geometry q ~ r^2/I_enc, so this is
          the exact Newton direction of the same fixed point; every label has
          direct leverage and convergence is geometric (<20 iterations).

        The ratio is clamped to [1/RATIO_CLAMP, RATIO_CLAMP] per iteration.
        Both variants share the fixed point q(psihat_m(l)) = q_t(psihat_m(l))."""
        if self.qt is None:
            raise ValueError("q matching requires a target q profile")
        if self.nprof is None:
            self.init_nprof_from_qt()

        from scipy.interpolate import PchipInterpolator
        best_err, best_N, n_stall = np.inf, None, 0

        # current-density weights of the classes for the per-class map
        cw = np.abs(self.cl.w * self.cl.v_par)
        cw = cw / cw.sum()

        for outer in range(1, self.max_it_out + 1):
            n_in, res_in = self.picard()
            ph, q_now = self.q_profile()
            # monotone cubic interpolation: piecewise-linear q(psihat) has
            # kinks that give the transplant a locally wrong response and
            # stall the iteration at the interpolation-error level
            q_i = PchipInterpolator(ph, q_now)

            # --- per-class label maps, combined as a current-weighted
            # geometric mean. A single effective-class map mis-models which
            # psihat the update at label l actually controls when the class
            # drift shifts differ strongly, and the outer iteration then
            # oscillates/diverges near the edge (M0 finding, 3-class case).
            # q_now and q_t are evaluated at the SAME clamped argument (an
            # unclamped target evaluation introduces a bias floor).
            log_ratio = np.zeros(len(self.nprof.l))
            q_at = np.zeros(len(self.nprof.l))
            qt_at = np.zeros(len(self.nprof.l))
            ph_ctl = ph[0]
            for s in range(self.cl.n_s):
                phm_s = self.label_to_psihat(self.nprof.l,
                                             alpha=self.cl.alpha[s])
                # q outside the beam edge is the vacuum annulus, and q
                # INSIDE the envelope roll-off is equally uncontrollable
                # (the envelope forces j -> 0 there): clamp evaluations to
                # the full-current beam interior
                l_ctl = 1.0 if self.l_beam >= 1.0 else \
                        max(self.l_beam - self.l_beam_width, 0.0)
                ph_beam = np.interp(l_ctl, self.nprof.l, phm_s)
                pe = np.clip(phm_s, ph[0], min(ph[-1], ph_beam))
                log_ratio += cw[s] * np.log(q_i(pe) / self.qt(pe))
                q_at += cw[s] * q_i(pe)
                qt_at += cw[s] * self.qt(pe)
                ph_ctl = max(ph_ctl, min(ph_beam, ph[-1]))

            if self.match_mode == 'q_shape':
                # compare shapes only; report the achieved amplitude
                c = np.sum(q_at * qt_at) / np.sum(qt_at**2)
                log_ratio = log_ratio - np.log(c)
            else:
                c = 1.0
            ratio = np.exp(log_ratio)

            # Convergence metric: compare q and q_t directly on the psihat
            # levels of the q diagnostic (this is the actual goal), restricted
            # to the range controllable through the label maps. The mapped
            # ratio above is only used to PLACE the update on the l grid.
            mctl = ph <= ph_ctl
            err = np.abs(q_now[mctl] / (c * self.qt(ph[mctl])) - 1.0).max()

            # --- best-iterate tracking and stagnation detection: with one
            # common Nprof and strongly different class maps, exactly
            # matching q_t can be outside the range of the ansatz; keep the
            # best profile and stop when no longer improving.
            if err < 0.98 * best_err:
                best_err, best_N, n_stall = err, self.nprof.N.copy(), 0
            else:
                n_stall += 1
            I_now = self.total_current(self.psi)
            self.log.append(dict(outer=outer, inner_iters=n_in,
                                 inner_res=res_in, q_err=err, I_RE=I_now,
                                 q_amplitude=c,
                                 lost=self.lost_fraction.max()))
            if self.verbose:
                print(f"  outer {outer:3d}: inner {n_in:3d} (res {res_in:.2e})"
                      f"  max|q/q_t-1| = {err:.3e}  I_RE = {I_now/1e6:8.4f} MA"
                      + (f"  q-ampl = {c:.4f}" if self.match_mode == 'q_shape' else "")
                      + (f"  lost = {self.lost_fraction.max():.2e}"
                         if self.lost_fraction.max() > 0 else ""))
            if err < self.tol_q or n_stall >= 15:
                # finishing pass: best profile + edge null-space polish of
                # Nprof (short-wavelength structure near l=1 is nearly
                # invisible to q but imprints oscillations on the edge
                # current density), then re-converge psi and re-evaluate
                if best_err < err:
                    self.nprof.N = best_N
                self._smooth_nprof()
                self.picard()
                ph, q_now = self.q_profile()
                mctl = ph <= ph_ctl
                err = np.abs(q_now[mctl] / (c * self.qt(ph[mctl])) - 1.0).max()
                if self.verbose:
                    print(f"  finishing after {outer} outer iterations: "
                          f"polished Nprof, final max|q/q_t-1| = {err:.3e}")
                self.log.append(dict(outer=outer, inner_iters=0, inner_res=0.0,
                                     q_err=err, I_RE=self.total_current(self.psi),
                                     q_amplitude=c, lost=self.lost_fraction.max()))
                return err < self.tol_q
            ratio = np.clip(ratio, 1.0 / self.RATIO_CLAMP, self.RATIO_CLAMP)
            # smooth the log-ratio (Nprof is smooth; single-point features in
            # the measured ratio are q-evaluation artifacts, and feeding them
            # to the transplant makes the iteration chase noise)
            lr = np.log(ratio)
            lr[1:-1] = 0.25 * lr[:-2] + 0.5 * lr[1:-1] + 0.25 * lr[2:]
            factor = np.exp(lr)**self.alpha_out
            if self.transplant == 'pointwise':
                self.nprof.N *= factor
            else:
                l = self.nprof.l
                dl = np.diff(l)
                C = np.concatenate([[0.0], np.cumsum(
                    0.5 * (self.nprof.N[1:] + self.nprof.N[:-1]) * dl)])
                C *= factor
                N_new = np.gradient(C, l, edge_order=2)
                self.nprof.N = np.clip(N_new, 0.0, None)
            self._apply_beam_envelope()
        if best_N is not None and best_err < np.inf:
            self.nprof.N = best_N
            self._smooth_nprof()
            self.picard()
        return False

    def _smooth_nprof(self, n_pass=6):
        """Remove the null-space ripple of Nprof: blended [1/4,1/2,1/4]
        smoothing, full strength towards l = 1, off below l = 0.5."""
        N = self.nprof.N
        w = np.clip((self.nprof.l - 0.5) / 0.3, 0.0, 1.0)
        for _ in range(n_pass):
            Ns = N.copy()
            Ns[1:-1] = 0.25*N[:-2] + 0.5*N[1:-1] + 0.25*N[2:]
            N = (1.0 - w)*N + w*Ns
        self.nprof.N = N
        self._apply_beam_envelope()

    def _apply_beam_envelope(self):
        """Confine the current to labels below l_beam (smoothstep roll-off
        over l_beam_width), leaving a vacuum annulus to the wall. Applied
        to the stored table so all consumers inherit it."""
        if self.l_beam >= 1.0:
            return
        t = np.clip((self.nprof.l - (self.l_beam - self.l_beam_width))
                    / self.l_beam_width, 0.0, 1.0)
        self.nprof.N = self.nprof.N * (1.0 - t*t*(3.0 - 2.0*t))

    def solve_fixed_nprof(self, nprof):
        """Inner solve only, for a prescribed Nprof (no q matching)."""
        self.nprof = nprof
        if not np.any(self.psi):
            self._init_psi_estimate()
        return self.picard()

    # --- initial guesses --------------------------------------------------------
    def _I_from_nprof_cyl(self):
        """Current implied by the present Nprof when the labels are
        approximated by the cylindrical psihat ~ (r/a)^2 (used only to build
        the initial psi estimate, before any labels exist)."""
        l_cyl = (self.gs.r[:, None] / self.gs.a)**2 * np.ones((1, self.gs.Nt))
        vbar_signed = np.sum(self.cl.w * self.cl.v_par)
        j = -EL_CHG * vbar_signed * self.nprof(l_cyl) / self.gs.RR
        dA = self.gs.r[:, None] * self.gs.dr * self.gs.dth
        return np.sum(j * dA)

    def _init_psi_estimate(self, I0=None):
        """Rough cylindrical psi to start the Picard iteration: parabolic
        current channel carrying ~the expected current. The sign follows the
        GS convention: j_phi > 0  =>  Delta* psi < 0  =>  psi max at axis."""
        if I0 is None:
            I0 = self.I_RE if self.I_RE is not None else self._I_from_nprof_cyl()
        if I0 == 0.0:
            I0 = 1.0e6
        a, R0 = self.gs.a, self.gs.R0
        rh = self.gs.r[:, None] / a
        # psi(r) for j ~ (1 - (r/a)^2), normalized to psi(a) = 0
        psi_prof = np.sign(I0) * MU_ZERO * abs(I0) * R0 / (4.0 * np.pi) \
            * (0.5 - rh**2 + 0.5 * rh**4)
        self.psi = np.broadcast_to(psi_prof, (self.gs.Nr, self.gs.Nt)).copy() \
            + self.gs.psi_b
        self.update_labels(self.psi)

    def init_nprof_from_qt(self, n_l=101):
        """Initial Nprof from the target q in cylindrical approximation:
        q = r B0 / (R0 Btheta)  =>  I(r) = 2 pi r^2 B0 / (mu0 R0 q(r)),
        j(r) = I'(r) / (2 pi r), Nprof^0(l(r)) ~ R0 j(r) / (e vbar).
        The label l is approximated by the cylindrical psihat(r). A crude
        guess only costs a few extra outer iterations."""
        B0 = abs(self.F0) / self.gs.R0
        a, R0 = self.gs.a, self.gs.R0
        r = np.linspace(a / 512, a, 512)
        # first pass: use psihat ~ (r/a)^2 as the argument of q_t
        for _ in range(2):
            ph = getattr(self, '_ph_cyl', (r / a)**2)
            q = self.qt(ph)
            Ienc = 2.0 * np.pi * r**2 * B0 / (MU_ZERO * R0 * q)
            Bth = MU_ZERO * Ienc / (2.0 * np.pi * r)
            psi_pol = np.concatenate([[0.0], np.cumsum(
                0.5 * (Bth[1:] + Bth[:-1]) * np.diff(r))]) * R0
            self._ph_cyl = psi_pol / psi_pol[-1]
        # analytic derivative of I = 2 pi r^2 B0/(mu0 R0 q):
        #   j = I'/(2 pi r) = (B0/(mu0 R0)) * (2/q - (r/q^2) dq/dr)
        # (dq/dr from the smooth q(r) array; avoids the r->0 blow-up of a
        # finite-difference I'/r evaluation)
        dq_dr = np.gradient(q, r, edge_order=2)
        j = B0 / (MU_ZERO * R0) * (2.0 / q - r * dq_dr / q**2)
        j = np.clip(j, 0.0, None)
        l_grid = np.linspace(0.0, 1.0, n_l)
        N0 = np.interp(l_grid, self._ph_cyl, j)
        vbar = np.abs(np.sum(self.cl.w * self.cl.v_par))
        N0 = N0 * R0 / (EL_CHG * vbar)
        self.nprof = Nprof(N0, l_grid)
        self._apply_beam_envelope()
        if not np.any(self.psi):
            # cylindrical q_t implies the current magnitude; the sign follows
            # from j_phi = -e sum_s v_par,s w_s Nprof/R with Nprof >= 0
            sign_j = -np.sign(np.sum(self.cl.w * self.cl.v_par))
            self._init_psi_estimate(I0=sign_j * Ienc[-1])

    # --- diagnostics -------------------------------------------------------------
    def drift_shifts(self):
        """Per-class diagnostics of the drift-surface shifts:
        - axis shift: distance between class drift axis and psi axis
        - boundary shift: inboard midplane gap between the drift surface
          through the outboard midplane edge and the plasma boundary."""
        R_ax, Z_ax, _ = self.psi_axis()
        out = []
        for s in range(self.cl.n_s):
            ax_shift = float(np.hypot(self.ax_RZ[s, 0] - R_ax,
                                      self.ax_RZ[s, 1] - Z_ax))
            # label of the outboard edge is 1 by construction; find where the
            # same drift surface crosses the inboard midplane
            Rg = np.linspace(self.gs.R0 - self.gs.a, self.ax_RZ[s, 0], 2001)
            l = self.Ahat(s, self.psi, Rg, np.zeros_like(Rg), clip=False)
            # innermost crossing of l = 1 on the inboard branch
            idx = np.where(l <= 1.0)[0]
            R_in = Rg[idx[0]] if len(idx) else np.nan
            bnd_shift = R_in - (self.gs.R0 - self.gs.a)
            out.append(dict(E_MeV=self.cl.E_kin[s] / 1e6,
                            axis_shift=ax_shift, boundary_shift=float(bnd_shift),
                            d_s=float(self.cl.d_shift(abs(self.F0) / self.gs.R0,
                                                      self.gs.a)[s]),
                            lost_fraction=float(self.lost_fraction[s])))
        return out

    def write_output(self, basename='re_eq_prototype'):
        """Write the converged state: psi map, Nprof table, per-class data,
        and the convergence log (plain text; the in-code version writes the
        JOREK restart plus an HDF5 per-class file)."""
        np.savez(basename + '_psi.npz', psi=self.psi, R=self.gs.RR, Z=self.gs.ZZ,
                 psi_b=self.gs.psi_b, F0=self.F0)
        with open(basename + '_classes.dat', 'w') as f:
            f.write("# s  E_kin[eV]  xi  weight  gamma  v_par[m/s]  "
                    "A_axis[Wb]  A_edge[Wb]  R_ax[m]  Z_ax[m]  lost_fraction\n")
            for s in range(self.cl.n_s):
                f.write(f"{s + 1:3d} {self.cl.E_kin[s]:14.6e} {self.cl.xi[s]:9.5f} "
                        f"{self.cl.w[s]:12.6e} {self.cl.gamma[s]:10.4f} "
                        f"{self.cl.v_par[s]:14.6e} {self.A_ax[s]:14.6e} "
                        f"{self.A_edge[s]:14.6e} {self.ax_RZ[s, 0]:10.5f} "
                        f"{self.ax_RZ[s, 1]:10.5f} {self.lost_fraction[s]:10.3e}\n")
        with open(basename + '_nprof.dat', 'w') as f:
            f.write("# l  Nprof(l) [m^-2]\n")
            for l, N in zip(self.nprof.l, self.nprof.N):
                f.write(f"{l:10.6f} {N:14.6e}\n")
        with open(basename + '_convergence.log', 'w') as f:
            f.write("# outer  inner_iters  inner_res  max|q/qt-1|  I_RE[A]  "
                    "q_amplitude  max_lost_fraction\n")
            for rec in self.log:
                f.write(f"{rec['outer']:5d} {rec['inner_iters']:6d} "
                        f"{rec['inner_res']:12.4e} {rec['q_err']:12.4e} "
                        f"{rec['I_RE']:14.6e} {rec['q_amplitude']:10.5f} "
                        f"{rec['lost']:12.4e}\n")


# ===========================================================================
# Command-line interface
# ===========================================================================

def main():
    import argparse
    p = argparse.ArgumentParser(
        description="Kinetic RE drift-surface equilibrium prototype (M0)")
    p.add_argument('--dist', required=True,
                   help="RE distribution table (columns: E_kin[eV] xi weight)")
    p.add_argument('--dist-format', default='ekin_xi_w',
                   choices=sorted(_DIST_READERS))
    p.add_argument('--q-target', help="target q profile (columns: psihat_n q)")
    p.add_argument('--R0', type=float, default=6.2)
    p.add_argument('--a', type=float, default=2.0)
    p.add_argument('--B0', type=float, default=5.3)
    p.add_argument('--I-RE', type=float, help="prescribed RE current [A] "
                   "(required for --match-mode q_shape)")
    p.add_argument('--match-mode', default='full_q',
                   choices=['full_q', 'q_shape'])
    p.add_argument('--Nr', type=int, default=128)
    p.add_argument('--Nt', type=int, default=256)
    p.add_argument('--xi-min', type=float, default=0.9)
    p.add_argument('--alpha-in', type=float, default=0.5)
    p.add_argument('--alpha-out', type=float, default=0.3)
    p.add_argument('--tol-q', type=float, default=1e-3)
    p.add_argument('--max-it-out', type=int, default=50)
    p.add_argument('--output', default='re_eq_prototype')
    args = p.parse_args()

    classes = REClasses.from_file(args.dist, fmt=args.dist_format,
                                  xi_min=args.xi_min)
    gs = PolarGSSolver(args.R0, args.a, Nr=args.Nr, Nt=args.Nt)
    qt = read_q_target(args.q_target) if args.q_target else None
    eq = REEquilibrium(classes, gs, F0=args.B0 * args.R0, q_target=qt,
                       match_mode=args.match_mode, I_RE=args.I_RE,
                       alpha_in=args.alpha_in, alpha_out=args.alpha_out,
                       tol_q=args.tol_q, max_it_out=args.max_it_out)

    B0 = args.B0
    print("Per-class drift parameters d_s = gamma m v_par / (e B0 a):")
    for s, d in enumerate(classes.d_shift(B0, args.a)):
        print(f"  class {s + 1}: E = {classes.E_kin[s] / 1e6:8.3f} MeV, "
              f"xi = {classes.xi[s]:+.3f}, w = {classes.w[s]:.4f}, "
              f"d_s = {d:+.4e}")

    if qt is None:
        raise SystemExit("A target q profile is required (--q-target); "
                         "for fixed-Nprof runs use the API directly.")
    converged = eq.match_q()
    eq.write_output(args.output)
    if not converged:
        raise SystemExit(f"ERROR: q matching did not converge to "
                         f"{args.tol_q} in {args.max_it_out} outer iterations "
                         f"(residual written to {args.output}_convergence.log)")
    print(f"Converged. Output written to {args.output}_*.")


if __name__ == '__main__':
    main()
