"""
OFC Earthquake Model - Analysis & Visualisation
================================================
Expects three output files from the Fortran simulation:
  ofc_avalanches.dat   - one integer per line (avalanche size per quake)
  ofc_stress.dat       - L x L whitespace-separated stress values
  ofc_footprint.dat    - L x L binary (0/1) mask of the largest avalanche

Usage:
  python plot_ofc.py                        # looks in current directory
  python plot_ofc.py --dir /path/to/data    # specify data directory
  python plot_ofc.py --save                 # save PNGs instead of showing

Author: Siim Erik Pugal, June 2026
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colors import LogNorm
from scipy import stats

# --------------------------------------------------------------------------- #
# CLI                                                                          #
# --------------------------------------------------------------------------- #

parser = argparse.ArgumentParser(description="OFC model diagnostics")
parser.add_argument("--dir",  default=".", help="Directory containing .dat files")
parser.add_argument("--save", action="store_true", help="Save figures to PNG instead of showing")
parser.add_argument("--burnin", type=float, default=0.2,
                    help="Fraction of quakes to discard as transient warm-up (default 0.2)")
args = parser.parse_args()

data_dir = Path(args.dir)

# --------------------------------------------------------------------------- #
# Load data                                                                    #
# --------------------------------------------------------------------------- #

def load(filename, **kwargs):
    path = data_dir / filename
    if not path.exists():
        print(f"ERROR: '{path}' not found. Run the Fortran code first.")
        sys.exit(1)
    return np.loadtxt(path, **kwargs)

print("Loading data...")
sizes_raw  = load("ofc_avalanches.dat", dtype=int)
stress     = load("ofc_stress.dat")
footprint  = load("ofc_footprint.dat", dtype=int)

L = stress.shape[0]
n_quakes = len(sizes_raw)

# Discard burn-in transient before steady state is reached
burnin_idx = int(args.burnin * n_quakes)
sizes = sizes_raw[burnin_idx:]
print(f"Grid:      {L} x {L}")
print(f"Quakes:    {n_quakes} total, {len(sizes)} after {int(args.burnin*100)}% burn-in")
print(f"Largest avalanche: {sizes.max()} sites ({100*sizes.max()/L**2:.1f}% of grid)")

# --------------------------------------------------------------------------- #
# Matplotlib style                                                             #
# --------------------------------------------------------------------------- #

plt.rcParams.update({
    "figure.facecolor":  "#0f1117",
    "axes.facecolor":    "#0f1117",
    "axes.edgecolor":    "#3a3d4a",
    "axes.labelcolor":   "#c8ccd8",
    "axes.titlecolor":   "#e8eaf0",
    "axes.grid":         True,
    "grid.color":        "#2a2d3a",
    "grid.linewidth":    0.5,
    "xtick.color":       "#8890a4",
    "ytick.color":       "#8890a4",
    "text.color":        "#c8ccd8",
    "font.family":       "sans-serif",
    "font.size":         11,
    "axes.titlesize":    13,
    "axes.titleweight":  "500",
    "axes.labelsize":    11,
    "legend.framealpha": 0.0,
    "legend.fontsize":   10,
    "lines.linewidth":   1.6,
})

TEAL   = "#5DCAA5"
AMBER  = "#EF9F27"
CORAL  = "#D85A30"
BLUE   = "#378ADD"
MUTED  = "#8890a4"

# --------------------------------------------------------------------------- #
# Figure layout                                                                #
# --------------------------------------------------------------------------- #

fig = plt.figure(figsize=(16, 5.5))
fig.patch.set_facecolor("#0f1117")
gs = gridspec.GridSpec(1, 3, figure=fig, wspace=0.38, left=0.07, right=0.97,
                       top=0.88, bottom=0.14)

ax1 = fig.add_subplot(gs[0])   # power-law distribution
ax2 = fig.add_subplot(gs[1])   # stress histogram
ax3 = fig.add_subplot(gs[2])   # avalanche footprint

# =========================================================================== #
# PLOT 1 - Avalanche size distribution (log-log)                               #
# =========================================================================== #

# Logarithmic binning: equal-width bins on a log scale avoid the artefact
# where large-s bins are nearly empty simply because they are narrow.
s_min = max(1, sizes.min())
s_max = sizes.max()
n_bins = 40
log_edges = np.logspace(np.log10(s_min), np.log10(s_max), n_bins + 1)
counts, edges = np.histogram(sizes, bins=log_edges)
bin_centres = np.sqrt(edges[:-1] * edges[1:])   # geometric centre

# Normalise to probability density P(s)
bin_widths = np.diff(edges)
prob = counts / (counts.sum() * bin_widths)

# Keep only non-zero bins for fitting and plotting
mask = counts > 0
s_plot = bin_centres[mask]
p_plot = prob[mask]

ax1.scatter(s_plot, p_plot, color=TEAL, s=28, zorder=3, label="log-binned P(s)")

# Power-law fit on the upper 60% of the size range (tail, away from transients)
fit_mask = s_plot > np.percentile(s_plot, 40)
if fit_mask.sum() >= 3:
    slope, intercept, r, _, _ = stats.linregress(
        np.log10(s_plot[fit_mask]), np.log10(p_plot[fit_mask])
    )
    s_fit = np.array([s_plot[fit_mask].min(), s_plot[fit_mask].max()])
    p_fit = 10**intercept * s_fit**slope
    ax1.plot(s_fit, p_fit, color=AMBER, lw=1.8, ls="--", zorder=4,
             label=f"power-law fit  τ = {-slope:.2f}")
    ax1.text(0.97, 0.97, f"slope = {slope:+.2f}\n$R^2$ = {r**2:.3f}",
             transform=ax1.transAxes, ha="right", va="top",
             fontsize=9.5, color=AMBER,
             bbox=dict(boxstyle="round,pad=0.3", fc="#1a1d26", ec="#3a3d4a", lw=0.5))

ax1.set_xscale("log")
ax1.set_yscale("log")
ax1.set_xlabel("Avalanche size  $s$  (sites fired)")
ax1.set_ylabel("Probability density  $P(s)$")
ax1.set_title("Avalanche size distribution")
ax1.legend(loc="lower left")

# Annotate Gutenberg-Richter reference slope −1.8
ax1.axvline(np.percentile(sizes, 95), color=MUTED, lw=0.6, ls=":")
ax1.text(np.percentile(sizes, 95)*1.05, p_plot.max()*0.4,
         "95th\npercentile", fontsize=8, color=MUTED, va="top")

# =========================================================================== #
# PLOT 2 - Stress histogram                                                    #
# =========================================================================== #

stress_flat = stress.ravel()
n_hist_bins = 60
counts_h, edges_h = np.histogram(stress_flat, bins=n_hist_bins, range=(0.0, 1.0))
centres_h = 0.5 * (edges_h[:-1] + edges_h[1:])
freq_h = counts_h / counts_h.sum()   # relative frequency

ax2.bar(centres_h, freq_h,
        width=edges_h[1] - edges_h[0],
        color=BLUE, alpha=0.75, zorder=3, label="stress values")

# Expected flat level if distribution were perfectly uniform
flat_level = 1.0 / n_hist_bins
ax2.axhline(flat_level, color=AMBER, lw=1.4, ls="--",
            label=f"uniform baseline  ({flat_level:.4f})")

# Mark mean and std
mu  = stress_flat.mean()
sig = stress_flat.std()
ax2.axvline(mu,      color=TEAL,  lw=1.2, ls="-",  label=f"mean = {mu:.3f}")
ax2.axvline(mu - sig, color=CORAL, lw=0.9, ls=":",  label=f"±1σ = {sig:.3f}")
ax2.axvline(mu + sig, color=CORAL, lw=0.9, ls=":")

ax2.set_xlabel("Local stress value")
ax2.set_ylabel("Relative frequency")
ax2.set_title("Stress field distribution")
ax2.set_xlim(0.0, 1.0)
ax2.legend(loc="upper left", fontsize=9)

# Summary statistics inset
stats_text = (
    f"sites: {L} x {L} = {L**2:,}\n"
    f"mean:  {mu:.4f}\n"
    f"std:   {sig:.4f}\n"
    f"min:   {stress_flat.min():.4f}\n"
    f"max:   {stress_flat.max():.4f}"
)
ax2.text(0.97, 0.97, stats_text,
         transform=ax2.transAxes, ha="right", va="top",
         fontsize=8.5, color=MUTED, family="monospace",
         bbox=dict(boxstyle="round,pad=0.4", fc="#1a1d26", ec="#3a3d4a", lw=0.5))

# =========================================================================== #
# PLOT 3 - Largest avalanche footprint                                         #
# =========================================================================== #

# Use a perceptually uniform sequential colormap on the footprint.
# footprint is 0/1; imshow will render 0=background, 1=fired.
fired_count = footprint.sum()
fired_pct   = 100.0 * fired_count / L**2

# Custom two-tone colormap: dark background + teal for fired sites
from matplotlib.colors import ListedColormap
fp_cmap = ListedColormap(["#0f1117", TEAL])

im = ax3.imshow(footprint, cmap=fp_cmap, vmin=0, vmax=1,
                interpolation="none", origin="upper",
                extent=[0, L, L, 0])

ax3.set_xlabel("Fault X coordinate")
ax3.set_ylabel("Fault Y coordinate")
ax3.set_title(f"Largest avalanche footprint\n"
              f"{fired_count:,} sites fired  ({fired_pct:.1f}% of grid)")

# Minimal colourbar
cbar = fig.colorbar(im, ax=ax3, shrink=0.6, pad=0.02, ticks=[0, 1])
cbar.ax.set_yticklabels(["resting", "fired"], fontsize=8.5, color=MUTED)
cbar.outline.set_edgecolor("#3a3d4a")

# =========================================================================== #
# Supertitle and export                                                         #
# =========================================================================== #

fig.suptitle(
    f"OFC Model diagnostics  —  {L} x {L} grid,  α = 0.2,  "
    f"{n_quakes:,} quakes  ({int(args.burnin*100)}% burn-in discarded)",
    fontsize=12, color="#e8eaf0", y=0.98
)

if args.save:
    out = data_dir / "ofc_diagnostics.png"
    fig.savefig(out, dpi=300, bbox_inches="tight", facecolor=fig.get_facecolor())
    print(f"Saved → {out}")
else:
    plt.show()