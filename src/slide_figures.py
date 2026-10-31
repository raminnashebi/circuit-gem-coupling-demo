"""Slide figures: one message per figure, drawn from results/results.json (run src/coupling_demo.py first).

Writes results/slides/:
  step1_limits.png    the host caps the circuit, and the cap falls with age
  step2_window.png    a safe operating window that closes with age
  step4_in_time.png   the host fails before the circuit does
  step6_measure.png   the data fix the product of expression and kcat, not each one
Uses Calibri when it is installed (as on a Mac with Microsoft Office), otherwise matplotlib's default font.
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import font_manager
from matplotlib.ticker import FixedLocator, NullLocator

ROOT = Path(__file__).resolve().parents[1]
R = json.loads((ROOT / "results" / "results.json").read_text())
OUT = ROOT / "results" / "slides"
OUT.mkdir(parents=True, exist_ok=True)

OFFICE_FONTS = Path("/Applications/Microsoft PowerPoint.app/Contents/Resources/DFonts")
if (OFFICE_FONTS / "Calibri.ttf").exists():
    for f in ("Calibri.ttf", "Calibrib.ttf"):
        font_manager.fontManager.addfont(str(OFFICE_FONTS / f))
    plt.rcParams["font.family"] = "Calibri"
plt.rcParams.update({"font.size": 15, "axes.titlesize": 15, "axes.labelsize": 15})

INK, MUTED, GRID, AXIS = "#1E2430", "#5B6472", "#E1E5EA", "#C9D0D9"
CRIMSON, TEAL = "#B5243C", "#008C99"  # validated pair (light surface)
TEAL_T, CRIMSON_T, GREY_T = "#CFE5E7", "#F4D9DE", "#ECEEF1"
SIZE = (7.75, 4.75)

NEED = R["level1_window"]["need_flux"]


def tidy(ax):
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(AXIS)
    ax.tick_params(colors=MUTED, labelsize=13, length=0, pad=5)
    ax.grid(axis="y", color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)


def crossing(age, y, level):
    """First age at which y falls below level (linear interpolation)."""
    below = np.nonzero(y < level)[0]
    if not len(below):
        return None
    i = below[0]
    return float(np.interp(level, [y[i], y[i - 1]], [age[i], age[i - 1]]))


def save(fig, name):
    fig.savefig(OUT / name, dpi=300, facecolor="white")
    plt.close(fig)
    print("wrote", OUT / name)


# Step 1: the host caps the circuit, and the cap falls with age
s = R["level1_capacity_series"]
age = np.array(s["age_d"])
fig, ax = plt.subplots(figsize=SIZE)
tidy(ax)
ax.axhline(NEED, color=MUTED, ls="--", lw=1.4)
ax.text(2, NEED + 0.02, "flux the circuit needs", ha="left", va="bottom", fontsize=13.5, color=MUTED)
for key, color, label, xy in (("ATP", TEAL, "ATP effector: hexokinase sets the cap", (40, 0.90)),
                              ("NADPH", CRIMSON, "NADPH effector: G6PD sets the cap", (4, 0.24))):
    y = np.array(s[key])
    ax.plot(age, y, color=color, lw=3)
    ax.text(*xy, label, color=color, fontsize=15, fontweight="bold")
    d = crossing(age, y, NEED)
    if d is not None:
        ax.plot([d], [NEED], "o", ms=9, color=color, mec="white", mew=1.5, zorder=5)
ax.set_xlim(0, 120)
ax.set_ylim(0, 1.25)
ax.set_xlabel("Cell age (days)", color=MUTED)
ax.set_ylabel("Most flux the cell\ncan give the circuit", color=INK)
fig.tight_layout()
save(fig, "step1_limits.png")

# Step 2: a safe operating window that closes with age
w = R["level1_window"]
ser = w["series"]
age, ceil, upper = np.array(ser["age_d"]), np.array(ser["ceiling"]), np.array(ser["upper_safe"])
close = w["closing_age_days"]
top = 0.75
end = crossing(age, ceil, NEED)  # after this age no enzyme level reaches the need
fig, ax = plt.subplots(figsize=SIZE)
tidy(ax)
ax.grid(False)
a = np.linspace(0, end, 600)
u = np.clip(np.interp(a, age, upper), NEED, top)
ax.fill_between(a, 0, NEED, color=GREY_T, lw=0)
ax.fill_between(a, NEED, u, color=TEAL_T, lw=0)
ax.fill_between(a, u, top, color=CRIMSON_T, lw=0)
ax.axvspan(end, 120, color=GREY_T, lw=0)
ax.plot(age, ceil, color=INK, lw=2)
ax.axvline(close, color=CRIMSON, ls=":", lw=2)
ax.text(close + 1.5, 0.56, f"window closes\nat about day {close:.0f}", color=CRIMSON, fontsize=14, va="top",
        fontweight="bold")
ax.text(4, 0.255, "safe operating\nwindow", color=TEAL, fontsize=18, fontweight="bold", va="center")
ax.text(3, 0.70, "too much: the reserve runs low", color=CRIMSON, fontsize=13.5)
ax.text(4, 0.05, "too little: the toxin is not cleared", color=INK, fontsize=13.5)
ax.annotate("most the network can supply", (90, float(np.interp(90, age, ceil))), xytext=(78, 0.42), fontsize=13,
            color=INK, arrowprops=dict(arrowstyle="-", color=INK, lw=0.9))
ax.set_xlim(0, 120)
ax.set_ylim(0, top)
ax.set_xlabel("Cell age (days)", color=MUTED)
ax.set_ylabel("Flux the enzyme could carry\n(kcat × enzyme amount)", color=INK)
fig.tight_layout()
save(fig, "step2_window.png")

# Step 4: the host fails before the circuit does
s2 = R["level2_series"]
t, X, reserve = np.array(s2["age_d"]), np.array(s2["X"]), np.array(s2["reserve"])
truth = R["level2_truth_design"]
lost, low = truth["toxin_loses_control_at_age_d"], truth["reserve_below_threshold_at_age_d"]
fig, (a1, a2) = plt.subplots(2, 1, figsize=SIZE, sharex=True)
for ax in (a1, a2):
    tidy(ax)
    ax.axvspan(-6, 0, color=CRIMSON_T, zorder=0)
a1.text(-5.5, 11.6, "Shaded: expression stage (circuit made)", fontsize=12.5, color=CRIMSON, clip_on=False)
a1.axhline(truth["x_without_circuit"], color=MUTED, ls=":", lw=1.3)
a1.text(119, truth["x_without_circuit"] - 0.5, "without the circuit", ha="right", va="top", fontsize=12.5,
        color=MUTED)
a1.axhline(4.0, color=MUTED, ls="--", lw=1.2)
a1.text(25, 4.4, "target", fontsize=12.5, color=MUTED)
a1.plot(t, X, color=INK, lw=2.6)
a1.plot([lost], [np.interp(lost, t, X)], "o", ms=9, color=INK, mec="white", mew=1.5, zorder=5)
a1.text(lost - 3, 7.0, f"toxin control lost: day {lost:.0f}", ha="right", fontsize=14.5, fontweight="bold",
        color=INK)
a1.set_ylim(0, 11)
a1.set_yticks([0, 5, 10])
a1.set_ylabel("Toxin", color=INK)
a2.axhline(R["level1_window"]["reserve_min_fold"], color=MUTED, ls="--", lw=1.2)
a2.text(119, R["level1_window"]["reserve_min_fold"] + 0.25, "minimum reserve", ha="right", va="bottom",
        fontsize=12.5, color=MUTED)
a2.plot(t, reserve, color=CRIMSON, lw=2.6)
a2.plot([low], [np.interp(low, t, reserve)], "o", ms=9, color=CRIMSON, mec="white", mew=1.5, zorder=5)
a2.text(low + 3, 4.4, f"host reserve gone: day {low:.0f}", fontsize=14.5, fontweight="bold", color=CRIMSON)
a2.set_ylim(0, 7.2)
a2.set_yticks([0, 2.5, 5])
a2.set_ylabel("Redox reserve", color=INK)
a2.set_xlim(-6, 120)
a2.set_xticks([0, 20, 40, 60, 80, 100, 120])
a2.set_xlabel("Cell age (days after the expression stage)", color=MUTED)
fig.tight_layout(h_pad=0.8)
save(fig, "step4_in_time.png")

# Step 6: the data fix the product of expression and kcat, not each one
ident = R["identifiability"]
rng = ident["log10_ranges"]
acc, acc2, true = ident["accepted_log10"], ident["accepted_with_enzyme_amount_log10"], ident["true_log10"]
fig, ax = plt.subplots(figsize=SIZE)
tidy(ax)
ax.grid(False)
ax.add_patch(plt.Rectangle((rng["a_E"][0], rng["kcat"][0]), rng["a_E"][1] - rng["a_E"][0],
                           rng["kcat"][1] - rng["kcat"][0], color=GREY_T, lw=0, zorder=0))
ax.scatter(acc["a_E"], acc["kcat"], s=34, color=TEAL, edgecolor="white", linewidth=0.6, zorder=3)
ax.scatter(acc2["a_E"], acc2["kcat"], s=46, color=CRIMSON, edgecolor="white", linewidth=0.6, zorder=4)
ax.plot([true["a_E"]], [true["kcat"]], marker="*", ms=20, color=INK, mec="white", mew=1.2, zorder=5)
ax.text(rng["a_E"][0] + 0.05, rng["kcat"][0] + 0.05, "all designs tried", va="bottom", fontsize=13, color=MUTED)
ax.text(-2.33, -0.30, f"fit the toxin data:\n{ident['n_accepted']} designs on one line", fontsize=14.5,
        fontweight="bold", color=TEAL, va="top")
ax.text(-1.45, 0.42, f"add the enzyme amount:\n{ident['n_accepted_with_enzyme_amount']} designs near the truth (star)",
        fontsize=14.5, fontweight="bold", color=CRIMSON, va="bottom")
ticks_x = [0.005, 0.01, 0.02, 0.05, 0.1, 0.2]
ticks_y = [0.2, 0.5, 1, 2, 5]
ax.xaxis.set_major_locator(FixedLocator(np.log10(ticks_x)))
ax.set_xticklabels([f"{v:g}" for v in ticks_x])
ax.yaxis.set_major_locator(FixedLocator(np.log10(ticks_y)))
ax.set_yticklabels([f"{v:g}" for v in ticks_y])
ax.xaxis.set_minor_locator(NullLocator())
ax.set_xlim(rng["a_E"][0] - 0.02, rng["a_E"][1] + 0.02)
ax.set_ylim(rng["kcat"][0] - 0.02, rng["kcat"][1] + 0.02)
ax.set_xlabel("Expression strength (log scale)", color=MUTED)
ax.set_ylabel("kcat (log scale)", color=INK)
fig.tight_layout()
save(fig, "step6_measure.png")
