"""Plots Monte Carlo error vs number of paths (log-log), with and without antithetic variates.

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


PLAIN, ANTI = "#2a78d6", "#eb6834"  # blue = plain MC, orange = antithetic

fig, ax = plt.subplots(figsize=(7, 4.5))

# Standard error is the expected size of the error, so it is drawn as a line.
# On log-log axes SE = s / sqrt(n) is a straight line with slope -1/2.
ax.plot(n, column("std_error"), color=PLAIN, linewidth=2, label="Plain MC, standard error")
ax.plot(n, column("av_std_error"), color=ANTI, linewidth=2, label="Antithetic, standard error")

# The absolute error |MC - BS| is one random draw per row, so it is drawn as points.
ax.plot(n, column("abs_error"), "o", color=PLAIN, markersize=8,
        markerfacecolor="none", markeredgewidth=2, label="Plain MC, |MC − BS|")
ax.plot(n, column("av_abs_error"), "s", color=ANTI, markersize=8,
        markerfacecolor="none", markeredgewidth=2, label="Antithetic, |MC − BS|")

ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel("Number of paths")
ax.set_ylabel("Error in option price")
ax.set_title("Monte Carlo convergence: ATM European call, seed 42")
ax.grid(True, which="major", color="#e4e3df", linewidth=0.8)
ax.set_axisbelow(True)
for side in ("top", "right"):
    ax.spines[side].set_visible(False)
ax.legend(frameon=False)

fig.tight_layout()
fig.savefig(png_path, dpi=150)
print(f"saved {png_path}")
