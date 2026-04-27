"""
KdV Bound States Solver
=======================
Finds bound states of the Schrödinger equation with a double delta-function potential:

    v(x,0) = -alpha * [delta(x - a) + delta(x + a)]

Each bound state energy E_i < 0 produces a KdV soliton with amplitude U_i = (2/3)|E_i|.

Usage:
    python kdv_bound_states.py
    python kdv_bound_states.py --alpha 3.0 --a 1.5
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.optimize import brentq
import argparse


# ── Transcendental equations ───────────────────────────────────────────────────

def even_eq(kappa, alpha, a):
    """Even-parity condition: kappa - alpha*(1 + exp(-2*kappa*a)) = 0"""
    return kappa - alpha * (1 + np.exp(-2 * kappa * a))


def odd_eq(kappa, alpha, a):
    """Odd-parity condition: kappa - alpha*(1 - exp(-2*kappa*a)) = 0"""
    return kappa - alpha * (1 - np.exp(-2 * kappa * a))


def find_bound_states(alpha, a, kappa_max=None, n_scan=10000):
    """
    Numerically find all bound states by scanning for sign changes
    and refining with bisection (brentq).

    Returns list of dicts: {kappa, E, parity}
    """
    if kappa_max is None:
        kappa_max = alpha * 2.5

    kappas = np.linspace(1e-6, kappa_max, n_scan)
    states = []

    for label, eq in [("even", even_eq), ("odd", odd_eq)]:
        vals = eq(kappas, alpha, a)
        for i in range(len(vals) - 1):
            if vals[i] * vals[i + 1] < 0:
                try:
                    kappa_sol = brentq(eq, kappas[i], kappas[i + 1], args=(alpha, a))
                    if kappa_sol > 1e-6:
                        states.append({
                            "kappa": kappa_sol,
                            "E": -kappa_sol**2,
                            "parity": label,
                        })
                except ValueError:
                    pass

    states.sort(key=lambda s: s["E"])
    return states


# ── Wavefunction ───────────────────────────────────────────────────────────────

def wavefunction(x, kappa, parity, a):
    """
    Analytic wavefunction for double-delta bound state.
    Normalised so that psi(a) = exp(-kappa*a) (continuity at boundary).
    """
    psi = np.zeros_like(x, dtype=float)
    boundary_val = np.exp(-kappa * a)

    if parity == "even":
        denom = np.cosh(kappa * a)
        A = boundary_val / denom if denom > 1e-12 else 0.0
        inside = np.abs(x) <= a
        psi[inside] = A * np.cosh(kappa * x[inside])
    else:
        denom = np.sinh(kappa * a)
        A = boundary_val / denom if abs(denom) > 1e-12 else 0.0
        inside = np.abs(x) <= a
        psi[inside] = A * np.sinh(kappa * x[inside])

    psi[x > a]  = np.exp(-kappa * x[x > a])
    psi[x < -a] = np.exp( kappa * x[x < -a])
    return psi


# ── Graphical equation plot ────────────────────────────────────────────────────

def plot_graphical_solution(ax, alpha, a, states):
    """
    Plot lhs = kappa and rhs curves for even/odd equations so the
    intersections (bound states) are visible.
    """
    kappa_max = alpha * 2.2
    kappas = np.linspace(1e-3, kappa_max, 500)

    rhs_even = alpha * (1 + np.exp(-2 * kappas * a))
    rhs_odd  = alpha * (1 - np.exp(-2 * kappas * a))

    ax.plot(kappas, kappas,    color="#555", lw=1.5, label=r"$\kappa$ (lhs)", ls="--")
    ax.plot(kappas, rhs_even,  color="#3266ad", lw=2,   label=r"$\alpha(1+e^{-2\kappa a})$ even")
    ax.plot(kappas, rhs_odd,   color="#c85820", lw=2,   label=r"$\alpha(1-e^{-2\kappa a})$ odd")

    colors_map = {"even": "#3266ad", "odd": "#c85820"}
    for s in states:
        ax.axvline(s["kappa"], color=colors_map[s["parity"]],
                   lw=1, ls=":", alpha=0.7)
        ax.plot(s["kappa"], s["kappa"], "o",
                color=colors_map[s["parity"]], ms=8, zorder=5)
        ax.annotate(
            f"  $\\kappa={s['kappa']:.3f}$\n  $E={s['E']:.3f}$",
            xy=(s["kappa"], s["kappa"]),
            fontsize=8, color=colors_map[s["parity"]],
        )

    ax.set_xlim(0, kappa_max)
    ax.set_ylim(0, kappa_max * 1.1)
    ax.set_xlabel(r"$\kappa$", fontsize=12)
    ax.set_ylabel("value", fontsize=12)
    ax.set_title("Graphical solution of transcendental equations", fontsize=11)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)


# ── Wavefunction plot ──────────────────────────────────────────────────────────

def plot_wavefunctions(ax, alpha, a, states):
    """Plot the wavefunctions for all bound states."""
    x = np.linspace(-4, 4, 800)

    colors = ["#3266ad", "#c85820", "#228b22", "#8b008b", "#8b4513"]
    for i, s in enumerate(states):
        psi = wavefunction(x, s["kappa"], s["parity"], a)
        label = (f"$\\psi_{i+1}$: {s['parity']}, "
                 f"$E={s['E']:.3f}$, "
                 f"$U={2/3*abs(s['E']):.3f}$")
        ax.plot(x, psi, lw=2, color=colors[i % len(colors)], label=label)

    # Mark delta-function positions
    for sign, xpos in [(+1, a), (-1, -a)]:
        ax.axvline(xpos, color="gray", lw=1, ls="--", alpha=0.5)
        ax.annotate(f"$\\delta(x{'+' if sign<0 else '-'}{a})$",
                    xy=(xpos, 0.85), fontsize=8, ha="center",
                    color="gray")

    ax.axhline(0, color="gray", lw=0.5)
    ax.set_xlim(-4, 4)
    ax.set_ylim(-1.3, 1.3)
    ax.set_xlabel("$x$", fontsize=12)
    ax.set_ylabel(r"$\psi(x)$", fontsize=12)
    ax.set_title("Bound-state wavefunctions", fontsize=11)
    ax.legend(fontsize=9, loc="upper right")
    ax.grid(True, alpha=0.3)


# ── Parameter scan ─────────────────────────────────────────────────────────────

def plot_energy_vs_alpha(ax, a):
    """Show how bound-state energies vary with alpha (for fixed a)."""
    alphas = np.linspace(0.1, 6, 300)
    even_energies, odd_energies = [], []

    for alp in alphas:
        sts = find_bound_states(alp, a)
        even = [s["E"] for s in sts if s["parity"] == "even"]
        odd  = [s["E"] for s in sts if s["parity"] == "odd"]
        even_energies.append(even[0] if even else np.nan)
        odd_energies.append(odd[0]  if odd  else np.nan)

    ax.plot(alphas, even_energies, color="#3266ad", lw=2, label="Even state $E_1$")
    ax.plot(alphas, odd_energies,  color="#c85820", lw=2, label="Odd state $E_2$")
    ax.axhline(0, color="gray", lw=0.8, ls="--")
    ax.set_xlabel(r"$\alpha$", fontsize=12)
    ax.set_ylabel("$E$", fontsize=12)
    ax.set_title(f"Bound-state energies vs. $\\alpha$ (fixed $a={a}$)", fontsize=11)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)


# ── Soliton profiles ───────────────────────────────────────────────────────────

def soliton_profile(x, t, kappa):
    """
    KdV 1-soliton: v(x,t) = -2*kappa^2 * sech^2(kappa*(x - 4*kappa^2*t))
    Speed c = 4*kappa^2, amplitude U = 2*kappa^2 = (2/3)*|E|  [with E = -kappa^2]
    """
    xi = kappa * (x - 4 * kappa**2 * t)
    return -2 * kappa**2 / np.cosh(xi)**2


def plot_solitons(ax, states):
    """Show the solitons that emerge from the bound states at a few times."""
    x = np.linspace(-10, 30, 1000)
    times = [0.0, 0.5, 1.5]
    cmaps = plt.cm.Blues, plt.cm.Oranges, plt.cm.Greens

    for i, s in enumerate(states[:3]):   # plot at most 3 solitons
        cmap = cmaps[i]
        for j, t in enumerate(times):
            prof = soliton_profile(x, t, s["kappa"])
            alpha_val = 0.4 + 0.6 * j / max(len(times) - 1, 1)
            label = (f"$U={2/3*s['kappa']**2:.3f}$, t={t}"
                     if i == 0 else (f"t={t}" if i == 0 else None))
            ax.plot(x, prof, color=cmap(0.5 + 0.3*j),
                    lw=1.5, alpha=alpha_val,
                    label=f"soliton {i+1}, t={t}" if j == 0 else None)

    ax.axhline(0, color="gray", lw=0.5)
    ax.set_xlabel("$x$", fontsize=12)
    ax.set_ylabel("$v(x,t)$", fontsize=12)
    ax.set_title("Emerging KdV solitons (t = 0, 0.5, 1.5)", fontsize=11)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)


# ── Main ───────────────────────────────────────────────────────────────────────

def main(alpha=2.0, a=1.0):
    states = find_bound_states(alpha, a)

    print("=" * 55)
    print(f"Double delta potential: alpha={alpha}, a={a}")
    print("=" * 55)
    if not states:
        print("No bound states found.")
    for i, s in enumerate(states):
        print(f"  State {i+1}: parity={s['parity']:4s}  "
              f"kappa={s['kappa']:.6f}  "
              f"E={s['E']:.6f}  "
              f"Soliton amplitude U={2/3*abs(s['E']):.6f}")
    print("=" * 55)

    fig = plt.figure(figsize=(14, 10))
    fig.suptitle(
        rf"KdV Inverse Scattering — Double $\delta$ potential  "
        rf"($\alpha={alpha}$, $a={a}$)",
        fontsize=13, fontweight="bold",
    )
    gs = gridspec.GridSpec(2, 2, hspace=0.45, wspace=0.35)

    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[1, 0])
    ax4 = fig.add_subplot(gs[1, 1])

    plot_graphical_solution(ax1, alpha, a, states)
    plot_wavefunctions(ax2, alpha, a, states)
    plot_energy_vs_alpha(ax3, a)
    if states:
        plot_solitons(ax4, states)
    else:
        ax4.text(0.5, 0.5, "No solitons\n(no bound states)",
                 ha="center", va="center", transform=ax4.transAxes, fontsize=12)
        ax4.set_title("KdV solitons", fontsize=11)

    plt.savefig("kdv_bound_states.png", dpi=150, bbox_inches="tight")
    print("Saved: kdv_bound_states.png")
    plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="KdV bound states solver")
    parser.add_argument("--alpha", type=float, default=2.0,
                        help="Delta-function strength (default: 2.0)")
    parser.add_argument("--a", type=float, default=1.0,
                        help="Half-separation of delta functions (default: 1.0)")
    args = parser.parse_args()
    main(alpha=args.alpha, a=args.a)
