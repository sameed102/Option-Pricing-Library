"""Plots Monte Carlo error vs number of paths (log-log) for plain Monte Carlo,
antithetic variates and the control variate.

Reads the CSV written by ./build/convergence and saves a PNG:
    python3 scripts/plot_convergence.py convergence.csv convergence.png
"""
import csv
import sys

import matplotlib

matplotlib.use("Agg")  # draw straight to a file, no window needed
import matplotlib.pyplot as plt

csv_path = sys.argv[1] if len(sys.argv) > 1 else "convergence.csv"
png_path = sys.argv[2] if len(sys.argv) > 2 else "convergence.png"

with open(csv_path, newline="") as f:
    rows = list(csv.DictReader(f))

n = [int(r["n_paths"]) for r in rows]


def column(name):
    return [float(r[name]) for r in rows]


# (label, SE column, |error| column, colour, marker). Colour and marker both identify
# the method, so the chart still reads without colour.
methods = [
    ("Plain MC", "std_error", "abs_error", "#2a78d6", "o"),
    ("Antithetic", "av_std_error", "av_abs_error", "#eb6834", "s"),
    ("Control variate", "cv_std_error", "cv_abs_error", "#1baf7a", "^"),
]

fig, ax = plt.subplots(figsize=(7, 4.8))
for label, se_col, err_col, colour, marker in methods:
    # Standard error is the expected size of the error, so it is drawn as a line.
    # On log-log axes SE = s / sqrt(n) is a straight line with slope -1/2.
    ax.plot(n, column(se_col), color=colour, linewidth=2, label=f"{label}, standard error")
    # The absolute error |MC - BS| is one random draw per row, so it is drawn as points.
    ax.plot(n, column(err_col), marker, color=colour, markersize=8, markerfacecolor="none",
            markeredgewidth=2, label=f"{label}, |MC − BS|")

ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel("Number of paths")
ax.set_ylabel("Error in option price")
ax.set_title("Monte Carlo convergence: ATM European call, seed 42")
ax.grid(True, which="major", color="#e4e3df", linewidth=0.8)
ax.set_axisbelow(True)
for side in ("top", "right"):
    ax.spines[side].set_visible(False)
# Below the axes so it never covers data; 3 columns = one column per method.
ax.legend(frameon=False, fontsize=8, ncol=3, loc="upper center", bbox_to_anchor=(0.5, -0.14))

fig.tight_layout()
fig.savefig(png_path, dpi=150, bbox_inches="tight")  # keep the legend inside the image
print(f"saved {png_path}")
