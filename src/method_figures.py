"""Method diagrams: how a gene circuit model is coupled to a metabolic network model, in three pictures.

Writes results/method/ (PNG and SVG):
  d1_network.png   steps 1 to 3: the circuit enters the network as reactions and bounds
  d2_coupling.png  step 4: the two models exchange a few variables at every time step
  d3_loop.png      steps 5 and 6: design, then measure, in a loop with the lab
These are schematics of the method, not results. Uses Calibri when it is installed, otherwise the default font.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.patches import Circle, FancyArrowPatch, FancyBboxPatch

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "method"
OUT.mkdir(parents=True, exist_ok=True)

OFFICE_FONTS = Path("/Applications/Microsoft PowerPoint.app/Contents/Resources/DFonts")
if (OFFICE_FONTS / "Calibri.ttf").exists():
    for f in ("Calibri.ttf", "Calibrib.ttf"):
        font_manager.fontManager.addfont(str(OFFICE_FONTS / f))
    plt.rcParams["font.family"] = "Calibri"
plt.rcParams.update({"font.size": 15, "svg.fonttype": "none"})

INK, MUTED, LINE = "#1E2430", "#5B6472", "#C9D0D9"
CRIMSON, CRIMSON_T = "#B5243C", "#F7E8EB"
TEAL, TEAL_T = "#008C99", "#DCEFF1"
CELL = "#F3F5F7"
W, H = 12.15, 4.55


def canvas():
    fig = plt.figure(figsize=(W, H))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, W)
    ax.set_ylim(0, H)
    ax.axis("off")
    return fig, ax


def box(ax, x, y, w, h, text, fc="white", ec=LINE, color=INK, size=15, bold=False, lw=1.4, r=0.12):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle=f"round,pad=0,rounding_size={r}", fc=fc, ec=ec, lw=lw,
                                zorder=2))
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", color=color, fontsize=size,
            fontweight="bold" if bold else "normal", zorder=3, linespacing=1.15)


def node(ax, x, y, text, fc="white", ec=MUTED, color=INK, r=0.36, size=14, bold=False):
    ax.add_patch(Circle((x, y), r, fc=fc, ec=ec, lw=1.6, zorder=3))
    ax.text(x, y, text, ha="center", va="center", fontsize=size, color=color, zorder=4,
            fontweight="bold" if bold else "normal")


def arrow(ax, p, q, color=MUTED, lw=2.2, style="-|>", rad=0.0, ms=16, ls="-", z=1):
    ax.add_patch(FancyArrowPatch(p, q, arrowstyle=style, mutation_scale=ms, color=color, lw=lw, ls=ls,
                                 connectionstyle=f"arc3,rad={rad}", shrinkA=0, shrinkB=0, zorder=z))


def badge(ax, x, y, n, r=0.24):
    ax.add_patch(Circle((x, y), r, fc=CRIMSON, ec="white", lw=1.5, zorder=6))
    ax.text(x, y - 0.01, str(n), ha="center", va="center", color="white", fontsize=16, fontweight="bold",
            zorder=7)


def note(ax, x, y, text, color=INK, size=14, ha="left", va="center", bold=False):
    ax.text(x, y, text, ha=ha, va=va, color=color, fontsize=size, fontweight="bold" if bold else "normal",
            zorder=8, linespacing=1.15)


def save(fig, name):
    for ext in ("png", "svg"):
        fig.savefig(OUT / f"{name}.{ext}", dpi=300, facecolor="white")
    plt.close(fig)
    print("wrote", OUT / f"{name}.png")


# ---------------------------------------------------------------------------------------------
# Diagram 1, steps 1 to 3: the circuit enters the network as reactions and bounds
# ---------------------------------------------------------------------------------------------
fig, ax = canvas()
ax.add_patch(FancyBboxPatch((1.35, 0.25), 6.95, 4.05, boxstyle="round,pad=0,rounding_size=0.35", fc=CELL, ec=LINE,
                            lw=1.6, zorder=0))
note(ax, 1.55, 0.40, "Host network: one stage-specific genome-scale model", color=MUTED, size=13)

GLC, G6P, LAC, LACX = (0.55, 3.25), (2.55, 3.25), (2.55, 1.05), (0.55, 1.05)
NADPH, ATP = (5.05, 3.25), (5.05, 1.45)
for xy, t in ((GLC, "glucose"), (LACX, "lactate")):
    node(ax, *xy, t, fc="white", r=0.48, size=13)
node(ax, *G6P, "G6P")
node(ax, *LAC, "lactate", r=0.44, size=13)
node(ax, *NADPH, "NADPH", fc=TEAL_T, ec=TEAL, bold=True, r=0.46)
node(ax, *ATP, "ATP", fc=TEAL_T, ec=TEAL, bold=True, r=0.46)

arrow(ax, (GLC[0] + 0.5, GLC[1]), (G6P[0] - 0.38, G6P[1]))
note(ax, 1.62, 3.48, "HK", color=TEAL, size=13, ha="center", bold=True)
arrow(ax, (G6P[0], G6P[1] - 0.38), (LAC[0], LAC[1] + 0.46))
note(ax, 2.42, 2.10, "glycolysis", color=MUTED, size=13, ha="right")
arrow(ax, (LAC[0] - 0.46, LAC[1]), (LACX[0] + 0.5, LACX[1]))
arrow(ax, (2.62, 1.75), (ATP[0] - 0.47, ATP[1] + 0.05), color=TEAL, lw=1.8, rad=0.0)
arrow(ax, (G6P[0] + 0.38, G6P[1]), (NADPH[0] - 0.48, NADPH[1]), color=TEAL, lw=4.2)
note(ax, 3.80, 3.52, "G6PD (PPP)", color=TEAL, size=13, ha="center", bold=True)
arrow(ax, (NADPH[0] + 0.3, NADPH[1] + 0.36), (6.15, 3.95), color=MUTED, lw=1.8)
note(ax, 6.22, 3.95, "glutathione", color=MUTED, size=13)
arrow(ax, (ATP[0] + 0.33, ATP[1] - 0.33), (6.15, 0.72), color=MUTED, lw=1.8)
note(ax, 6.22, 0.70, "ion pumps", color=MUTED, size=13)

# the effector reaction, added by the circuit
box(ax, 6.05, 2.05, 2.05, 0.95, "effector\ntoxin → product", fc=CRIMSON_T, ec=CRIMSON, color=CRIMSON, size=14,
    bold=True, lw=2.4)
arrow(ax, (NADPH[0] + 0.40, NADPH[1] - 0.25), (6.05, 2.70), color=CRIMSON, lw=3.2, rad=-0.15)

# the circuit, outside the network, sets the effector's bound through a gene rule
box(ax, 9.05, 1.55, 2.95, 1.95, "", fc="white", ec=CRIMSON, lw=2.0)
note(ax, 10.52, 3.20, "Circuit", color=CRIMSON, size=15, ha="center", bold=True)
note(ax, 10.52, 2.65, "transgene → enzyme E", size=14, ha="center")
note(ax, 10.52, 2.05, "v  ≤  kcat · E", color=CRIMSON, size=17, ha="center", bold=True)
arrow(ax, (9.05, 2.52), (8.12, 2.52), color=CRIMSON, lw=2.6)
note(ax, 8.58, 2.80, "gene rule", color=CRIMSON, size=12.5, ha="center")

# step badges and their one-line labels
badge(ax, 9.30, 3.85, 1)
note(ax, 9.62, 3.85, "add the circuit: reactions,\ngene rule, flux bound", size=13.5)
badge(ax, 3.08, 2.42, 2)
note(ax, 3.38, 2.42, "linear programming: the ceiling\non v, and the enzyme that limits it", size=13)
badge(ax, 0.55, 2.40, 3)
note(ax, 0.55, 1.97, "data become\nbounds", size=13.5, ha="center")
for xy in ((1.62, 2.98), (1.62, 1.30)):
    note(ax, *xy, "[min, max]", color=TEAL, size=12, ha="center", bold=True)
save(fig, "d1_network")

# ---------------------------------------------------------------------------------------------
# Diagram 2, step 4: the two models exchange a few variables at every time step
# ---------------------------------------------------------------------------------------------
fig, ax = canvas()
# circuit model
box(ax, 0.25, 1.75, 4.35, 2.55, "", fc=CRIMSON_T, ec=CRIMSON, lw=2.0, r=0.25)
note(ax, 2.42, 3.95, "Circuit model (ODE)", color=CRIMSON, size=16, ha="center", bold=True)
for i, t in enumerate(("sensor", "regulator", "enzyme E")):
    box(ax, 0.45 + i * 1.38, 2.70, 1.20, 0.70, t, fc="white", ec=CRIMSON, size=14)
    if i:
        arrow(ax, (0.45 + i * 1.38 - 0.18, 3.05), (0.45 + i * 1.38, 3.05), color=CRIMSON, lw=1.8, ms=13, z=4)
note(ax, 2.42, 2.20, "toxin level feeds the sensor", color=INK, size=13.5, ha="center")

# network model
box(ax, 7.55, 1.75, 4.35, 2.55, "", fc=CELL, ec=LINE, lw=2.0, r=0.25)
note(ax, 9.72, 3.95, "Network model (linear program)", color=INK, size=16, ha="center", bold=True)
dx = 0.15
for (x, y) in ((8.05, 3.05), (9.20, 3.35), (9.20, 2.55), (10.40, 3.05), (11.25, 2.55)):
    node(ax, x + dx, y, "", r=0.17)
for p, q in (((8.22, 3.10), (9.03, 3.33)), ((8.22, 3.00), (9.03, 2.60)), ((9.37, 3.30), (10.23, 3.10)),
             ((9.37, 2.60), (10.23, 3.00)), ((10.55, 2.95), (11.08, 2.62))):
    arrow(ax, (p[0] + dx, p[1]), (q[0] + dx, q[1]), color=MUTED, lw=1.5, ms=11, z=4)
arrow(ax, (9.37 + dx, 2.55), (11.08 + dx, 2.55), color=CRIMSON, lw=2.6, ms=12, z=4)
note(ax, 9.72, 2.12, "solved again at every time step", color=MUTED, size=13.5, ha="center")

# the exchange
arrow(ax, (4.65, 3.55), (7.50, 3.55), color=CRIMSON, lw=3.2, ms=20)
note(ax, 6.08, 3.88, "enzyme level E(t)", color=CRIMSON, size=14.5, ha="center", bold=True)
note(ax, 6.08, 3.25, "sets the bound  v ≤ kcat·E", color=INK, size=13.5, ha="center")
arrow(ax, (7.50, 2.30), (4.65, 2.30), color=TEAL, lw=3.2, ms=20)
note(ax, 6.08, 2.62, "flux delivered, ATP, NADPH", color=TEAL, size=14.5, ha="center", bold=True)
note(ax, 6.08, 1.98, "toxin removed, expression cost", color=INK, size=13.5, ha="center")
badge(ax, 6.08, 4.30, 4)

# the two stages of the red cell
box(ax, 0.25, 0.25, 4.80, 1.05, "Erythroblast: the circuit is made\nexpression costs ATP", fc="white", ec=CRIMSON,
    size=14)
arrow(ax, (5.15, 0.78), (6.95, 0.78), color=INK, lw=2.4, ms=18)
note(ax, 6.05, 1.05, "maturation", color=MUTED, size=13, ha="center")
box(ax, 7.05, 0.25, 4.85, 1.05, "Mature red cell: no ribosomes\nonly a decaying enzyme pool acts", fc="white", ec=INK,
    size=14)
save(fig, "d2_coupling")

# ---------------------------------------------------------------------------------------------
# Diagram 3, steps 5 and 6: design, then measure, in a loop with the lab
# ---------------------------------------------------------------------------------------------
fig, ax = canvas()
bw, bh = 3.15, 1.05
top_y, bot_y = 3.20, 0.30
xs = (0.35, 4.50, 8.65)
model_kw = dict(fc=CRIMSON_T, ec=CRIMSON, size=14.5)
lab_kw = dict(fc=TEAL_T, ec=TEAL, size=14.5)
box(ax, xs[0], top_y, bw, bh, "Sample circuit designs\nexpression, kcat, threshold", **model_kw)
box(ax, xs[1], top_y, bw, bh, "Simulate each design\nwith the coupled model", **model_kw)
box(ax, xs[2], top_y, bw, bh, "Rank by robustness\nacross ages and scenarios", **model_kw)
box(ax, xs[2], bot_y, bw, bh, "Identifiability: what\nthe data cannot fix", **model_kw)
box(ax, xs[1], bot_y, bw, bh, "Lab measures the most\ninformative quantity", **lab_kw)
box(ax, xs[0], bot_y, bw, bh, "New data become\nbounds (step 3)", **lab_kw)
for a, b in ((0, 1), (1, 2)):
    arrow(ax, (xs[a] + bw + 0.05, top_y + bh / 2), (xs[b] - 0.05, top_y + bh / 2), color=INK, lw=2.4, ms=18)
for a, b in ((2, 1), (1, 0)):
    arrow(ax, (xs[a] - 0.05, bot_y + bh / 2), (xs[b] + bw + 0.05, bot_y + bh / 2), color=INK, lw=2.4, ms=18)
arrow(ax, (xs[2] + bw / 2, top_y - 0.05), (xs[2] + bw / 2, bot_y + bh + 0.05), color=INK, lw=2.4, ms=18)
arrow(ax, (xs[0] + bw / 2, bot_y + bh + 0.05), (xs[0] + bw / 2, top_y - 0.05), color=INK, lw=2.4, ms=18)
note(ax, W / 2, 2.27, "design  →  build  →  test", color=INK, size=20, ha="center", bold=True)
note(ax, W / 2, 1.78, "the best designs go to the lab with the measurement that separates them",
     color=MUTED, size=14, ha="center")
badge(ax, xs[0] + 0.05, top_y + bh - 0.02, 5)
badge(ax, xs[2] + 0.05, bot_y + bh - 0.02, 6)
save(fig, "d3_loop")
