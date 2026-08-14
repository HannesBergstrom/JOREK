#!/usr/bin/env python3
"""Per-circuit vacuum field response at the magnetic axis.

Answers: for unit current in each CIRCUIT, what is B_z at the axis (which moves
the plasma radially) and what is dB_z/dR there (the field decay index, which
sets plasma size at fixed position and current)?

That split is the whole point. A uniform scaling of all coil currents is mostly
a B_z knob, so it displaces the plasma and the radial position feedback then
cancels it -- measured on the 10 MeV JET case, a 40% swing of every coil bought
2.6% of q amplitude (gain ~0.053, saturating) while c_glob = 1 needed a scale of
~2.0, at which the free-boundary iteration no longer converges at all. What is
wanted instead is a POSITION-NEUTRAL combination: change the gradient while
leaving B_z at the axis alone, so the plasma compresses instead of moving and
the position loop is never excited.

Circuits, not coils: on JET 20 coils are wired into 10 circuits, so coil
currents are not independently drivable and any actuator must live in circuit
space or it asks the machine for something it cannot do.

Physics is the same as JOREK's B_coil_unit (mod_plasma_response.f90): sum the
circular-loop Green's functions over filaments, weighted by
weight(i_f) = n_turns_fila / n_turns_coil, times mu_zero.

  ./coil_circuit_response.py JET_coils.txt --axis 2.59 -0.0946
  ./coil_circuit_response.py coils.txt --axis 2.59 -0.09 --circuits my.csv
"""

import argparse
import csv
import sys

import numpy as np
from scipy.special import ellipk, ellipe

MU0 = 4.0e-7 * np.pi

# JET: alpha(circuit, coil) transcribed from find_Icoils_JET
# (models/mod_plasma_response.f90). Entries are turns, sign = winding sense.
JET_ALPHA = {
    (1, 1): 710.0, (2, 2): 426.0,
    (3, 3): -8.0, (3, 4): -20.0, (3, 5): -8.0, (3, 6): -20.0,
    (4, 7): -8.0, (4, 8): 8.0,
    (3, 9): 30.0, (3, 10): 30.0,
    (4, 11): -20.0, (4, 12): 20.0,
    (1, 13): 2.0, (1, 14): 2.0,
    (5, 15): 61.0, (5, 16): 61.0,
    (6, 15): -61.0, (6, 16): 61.0,
    (7, 17): 15.99, (8, 18): 15.0, (9, 19): 15.0, (10, 20): 21.0,
}
JET_LABELS = {
    1: "ohmic (P1/ME + P3/MU + P3/ML)",
    2: "P1/MC",
    3: "shaping, up/down SYM (-P2 +P3)",
    4: "shaping, up/down ANTISYM",
    5: "P4 symmetric  -> vertical field",
    6: "P4 antisymmetric -> radial field",
    7: "D1", 8: "D2", 9: "D3", 10: "D4",
}


def read_coils_geo(path):
    """Parse JOREK coils_geo.txt exactly as read_coils does."""
    with open(path) as fh:
        lines = [l for l in fh.read().split("\n")]
    k = 0

    def nxt():
        nonlocal k
        v = lines[k]
        k += 1
        return v

    nxt()
    ncoils = int(nxt().split()[0])
    nxt()
    hdr = []
    for _ in range(ncoils):
        parts = nxt().split('"')
        hdr.append((int(parts[0].split()[0]), parts[1].strip(),
                    float(parts[2].split()[0])))
    coils = []
    for n_fila, name, n_turns in hdr:
        nxt(); nxt(); nxt()
        R = np.empty(n_fila); Z = np.empty(n_fila); W = np.empty(n_fila)
        for i in range(n_fila):
            a, b, c = nxt().split()[:3]
            R[i], Z[i], W[i] = float(a), float(b), float(c) / n_turns
        coils.append(dict(name=name, R=R, Z=Z, w=W, n_turns=n_turns))
    return coils


def b_coil_unit(coil, R0, Z0):
    """(B_R, B_Z) at (R0,Z0) for unit current, mirroring B_coil_unit."""
    Rp, Zp, w = coil["R"], coil["Z"], coil["w"]
    rho2 = (Rp + R0) ** 2 + (Zp - Z0) ** 2
    kk = np.sqrt(4.0 * Rp * R0 / rho2)
    K, E = ellipk(kk ** 2), ellipe(kk ** 2)      # scipy takes m = k^2
    d2 = (Rp - R0) ** 2 + (Zp - Z0) ** 2
    g_br = (0.5 / np.pi) / np.sqrt(rho2) * (Z0 - Zp) / R0 * \
           ((Rp ** 2 + R0 ** 2 + (Zp - Z0) ** 2) / d2 * E - K)
    g_bz = (0.5 / np.pi) / np.sqrt(rho2) * \
           ((Rp ** 2 - R0 ** 2 - (Zp - Z0) ** 2) / d2 * E + K)
    # signs and mu_zero as in B_coil_unit
    return -float(np.sum(g_br * w)) * MU0, -float(np.sum(g_bz * w)) * MU0


def load_alpha(path, n_coils):
    """circuit,coil,turns CSV -> {(circuit,coil): turns}."""
    alpha = {}
    with open(path) as fh:
        for row in csv.reader(fh):
            if not row or row[0].lstrip().startswith("#"):
                continue
            c, k, t = int(row[0]), int(row[1]), float(row[2])
            if not 1 <= k <= n_coils:
                sys.exit("ERROR: %s references coil %d, file has %d" % (path, k, n_coils))
            alpha[(c, k)] = t
    return alpha


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("coils", help="JOREK coils_geo.txt")
    p.add_argument("--axis", nargs=2, type=float, required=True,
                   metavar=("R", "Z"), help="magnetic axis position")
    p.add_argument("--dR", type=float, default=0.01,
                   help="finite-difference step for dB_z/dR [m] (default 0.01)")
    p.add_argument("--circuits", default=None,
                   help="circuit matrix CSV 'circuit,coil,turns'; "
                        "default: the JET wiring from find_Icoils_JET")
    p.add_argument("--per-coil", action="store_true",
                   help="also list the individual coils")
    p.add_argument("--symmetric", action="store_true",
                   help="restrict the actuator to up/down SYMMETRIC circuits, so "
                        "it cannot distort the plasma vertically. Weaker than the "
                        "unrestricted optimum but still far stronger than a "
                        "uniform scale, and safe for an up/down symmetric target.")
    p.add_argument("--exclude", type=int, nargs="+", default=[], metavar="C",
                   help="circuits unavailable as actuators (e.g. the ohmic "
                        "transformer circuits, which drive the plasma current "
                        "and are not free for shaping)")
    args = p.parse_args()

    R0, Z0 = args.axis
    coils = read_coils_geo(args.coils)
    n = len(coils)
    print("coil file : %s  (%d coils, %d filaments)"
          % (args.coils, n, sum(len(c["R"]) for c in coils)))
    print("axis      : R = %.4f  Z = %.4f\n" % (R0, Z0))

    # per-coil response and its radial derivative
    bz = np.empty(n); dbz = np.empty(n); br = np.empty(n)
    for i, c in enumerate(coils):
        br[i], bz[i] = b_coil_unit(c, R0, Z0)
        _, bz_p = b_coil_unit(c, R0 + args.dR, Z0)
        _, bz_m = b_coil_unit(c, R0 - args.dR, Z0)
        dbz[i] = (bz_p - bz_m) / (2.0 * args.dR)

    if args.per_coil:
        print("per COIL, unit current:")
        print("  #  name            B_z [T/A]      dB_z/dR [T/A/m]")
        for i, c in enumerate(coils):
            print("  %2d %-14s %+.6e   %+.6e" % (i + 1, c["name"][:14], bz[i], dbz[i]))
        print()

    if args.circuits:
        alpha = load_alpha(args.circuits, n)
        labels = {}
    else:
        if n != 20:
            sys.exit("ERROR: the built-in circuit matrix is JET's (20 coils); "
                     "this file has %d. Supply --circuits." % n)
        alpha, labels = JET_ALPHA, JET_LABELS
        print("circuit matrix: built-in JET wiring (find_Icoils_JET)\n")

    n_circ = max(c for c, _ in alpha)
    A = np.zeros((n_circ, n))
    for (c, k), t in alpha.items():
        A[c - 1, k - 1] = t

    # unit CIRCUIT current -> coil currents are alpha(c,:) [turns]
    Bz_c = A @ bz
    dBz_c = A @ dbz
    Br_c = A @ br

    print("per CIRCUIT, unit circuit current:")
    print("  #   B_z [T/A]      dB_z/dR [T/A/m]   n_index    |B_r| [T/A]   description")
    for c in range(n_circ):
        # decay index n = -(R/B_z) dB_z/dR; only meaningful when B_z is not ~0
        nidx = (-R0 / Bz_c[c] * dBz_c[c]) if abs(Bz_c[c]) > 1e-12 else float("nan")
        print("  %2d  %+.6e   %+.6e   %8.3f   %.3e   %s"
              % (c + 1, Bz_c[c], dBz_c[c], nidx, abs(Br_c[c]),
                 labels.get(c + 1, "")))

    # --- position-neutral shaping directions -----------------------------
    # A combination x (in circuit space) is position-neutral at the axis when
    # it produces no net B_z (no radial force change) and no net B_r (no
    # vertical force change). Among those, we want the largest |dB_z/dR|:
    # that is the compression knob. Solve by projecting the gradient vector
    # onto the null space of the two position constraints.
    # circuits that are not available as actuators are pinned to zero by
    # adding a unit row per excluded circuit to the constraint matrix -- the
    # null space then contains only directions that leave them untouched
    # B_r = 0 pins vertical force; B_z = 0 pins radial force. Under --symmetric
    # the B_r constraint is dropped: an up/down symmetric circuit produces no
    # vertical force by construction, and the small residual B_r seen here is
    # only because the axis sits slightly off the midplane. Imposing it exactly
    # would leave 2 free circuits with 2 constraints and hence no actuator at
    # all -- the constraint would be spending the entire subspace on an effect
    # the symmetry restriction has already removed.
    C = np.vstack([Bz_c]) if args.symmetric else np.vstack([Bz_c, Br_c])
    if args.symmetric:
        # Pin every up/down ANTISYMMETRIC circuit to zero. Symmetry is decided
        # from the wiring itself: for each coil in the circuit, find its mirror
        # partner across Z = 0 and compare winding signs. A circuit whose coils
        # pair up with OPPOSITE signs drives a radial field / vertical motion
        # and would tilt the plasma; one with matching signs cannot.
        Zc = np.array([c["Z"].mean() for c in coils])
        Rc = np.array([c["R"].mean() for c in coils])
        for c in range(n_circ):
            memb = [k for k in range(n) if abs(A[c, k]) > 1e-12]
            if not memb:
                continue
            tot = 0.0
            for k in memb:
                d = [abs(Zc[j] + Zc[k]) + abs(Rc[j] - Rc[k]) for j in memb]
                j = memb[int(np.argmin(d))]
                if min(d) < 0.35:
                    tot += np.sign(A[c, k]) * np.sign(A[c, j])
            if tot <= 0.0:                       # antisymmetric or unpaired
                row = np.zeros(n_circ); row[c] = 1.0
                C = np.vstack([C, row])
        print("\nrestricted to up/down symmetric circuits")
    for c in args.exclude:
        if not 1 <= c <= n_circ:
            sys.exit("ERROR: --exclude %d outside 1..%d" % (c, n_circ))
        row = np.zeros(n_circ); row[c - 1] = 1.0
        C = np.vstack([C, row])
    if args.exclude:
        print("\nexcluded from the actuator: circuits %s"
              % ", ".join(str(c) for c in sorted(args.exclude)))
    u, s, vt = np.linalg.svd(C)
    rank = int((s > 1e-12 * max(s.max(), 1e-30)).sum())
    N = vt[rank:].T                                  # null-space basis
    print("\nposition-neutral subspace: %d of %d circuit directions"
          % (N.shape[1], n_circ))
    if N.shape[1] == 0:
        print("  (none -- every circuit combination moves the plasma)")
        return

    g = N.T @ dBz_c                                  # gradient within the null space
    if np.linalg.norm(g) < 1e-30:
        print("  no direction in it changes dB_z/dR")
        return
    x = N @ (g / np.linalg.norm(g))
    x = x / np.max(np.abs(x))                        # normalise to unit peak circuit current

    print("  best compression direction (normalised to peak circuit current 1):")
    for c in range(n_circ):
        if abs(x[c]) > 1e-6:
            print("    circuit %2d : %+8.4f   %s" % (c + 1, x[c], labels.get(c + 1, "")))
    print("  yields  dB_z/dR = %+.6e T/A/m   with  B_z = %+.2e  B_r = %+.2e (both ~0)"
          % (x @ dBz_c, x @ Bz_c, x @ Br_c))

    # how it compares to the uniform-scale knob, which is what we measured
    ones = np.ones(n_circ)
    print("\n  for reference, all circuits at unit current:")
    print("    dB_z/dR = %+.6e   B_z = %+.6e  (this is the knob that failed)"
          % (ones @ dBz_c, ones @ Bz_c))

    # --- per-COIL actuator vector -------------------------------------------
    # Coils wired into one circuit all carry that circuit's current; the turns
    # count only sets how much field it makes, which is already in alpha. So
    # alpha(c,k) = +/- n_turns(k), and converting a circuit direction to coil
    # currents in A/turn is just the winding sign:
    #     amp(k) = sum_c x(c) * sign(alpha(c,k))
    # Verified against the JET table: every non-zero alpha equals +/- the
    # coil's n_turns, and circuit 5 reproduces the rad_FB_amp pattern the user
    # already has on coils 15/16.
    amp = np.zeros(n)
    for (c, k), t in alpha.items():
        amp[k - 1] += x[c - 1] * np.sign(t)
    peak = np.max(np.abs(amp))
    if peak > 0:
        amp = amp / peak          # unit peak coil current, control carries the scale

    print("\n  per-COIL actuator vector [A/turn per unit control], "
          "normalised to peak 1:")
    line = "  re_eq_coil_amp = " + " ".join("%.6f" % v for v in amp)
    for i in range(0, len(line), 100):
        print(line[i:i + 100])
    nz = [(i + 1, coils[i]["name"], amp[i]) for i in range(n) if abs(amp[i]) > 1e-9]
    print("  non-zero on %d of %d coils:" % (len(nz), n))
    for i, nm, v in nz:
        print("    %2d %-14s %+8.4f" % (i, nm[:14], v))

    with open("re_coil_amp.txt", "w") as fh:
        fh.write("# per-coil actuator direction, A/turn per unit control\n")
        fh.write("# position-neutral at R=%.4f Z=%.4f; excluded circuits: %s\n"
                 % (R0, Z0, sorted(args.exclude) if args.exclude else "none"))
        fh.write("# dB_z/dR = %+.6e T/A/m at B_z = %+.2e, B_r = %+.2e\n"
                 % (x @ dBz_c, x @ Bz_c, x @ Br_c))
        fh.write(line.strip() + "\n")
    print("  written to re_coil_amp.txt")


if __name__ == "__main__":
    main()
