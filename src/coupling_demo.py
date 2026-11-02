#!/usr/bin/env python3
"""Toy demonstration: coupling a gene circuit model to a metabolic network model.

Everything here is invented for teaching. The network is a small stylised red cell
(glycolysis, a pentose phosphate route, pumps and a redox demand) and every number,
including the "measurements", is made up or simulated. Nothing below is a result
about real red cells.

Pipeline:

  Level 1, stoichiometric coupling
    1. Add the circuit's effector reaction to the network, tied to a transgene by a
       gene rule, and bound its flux by the enzyme level (v <= kcat * E).
    2. Use linear programming to ask whether the effector is feasible, which
       resource limits it (marginal values), and where the safe operating window is.
    3. Constrain the network with (synthetic) exchange rates, tracer fractions and
       enzyme activities; use flux sampling instead of a single objective.
  Level 2, dynamic coupling
    4. Couple an ODE circuit model to the network in time steps: enzyme level sets
       the flux bound, the network returns the flux the cell can really deliver,
       and the toxin level feeds back to the sensor.
    5. Search circuit designs as an ensemble with a fast lookup of exact LP results.
    6. Use synthetic time course data to see what is and is not identifiable, and
       which measurement would resolve it.

Run from the repository root:  python src/coupling_demo.py   (or: make results)
Writes results/results.json and results/figures/*.png.
"""
from __future__ import annotations

import json
import time
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import ListedColormap
from scipy.interpolate import RegularGridInterpolator
from scipy.linalg import null_space
from scipy.optimize import linprog
from scipy.stats import qmc
from scipy.stats import t as student_t

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results"
FIG = OUT / "figures"
FIG.mkdir(parents=True, exist_ok=True)
SEED = 2026
T0 = time.time()

TEAL, RED, GREY, INK = "#1D6A73", "#9E1B32", "#8A93A0", "#1E2430"
TEAL_L, RED_L, GREY_L = "#C5DFE1", "#F4C7CF", "#E3E7EC"
plt.rcParams.update({
    "font.size": 10, "axes.spines.top": False, "axes.spines.right": False,
    "axes.titlesize": 11, "axes.titleweight": "bold", "axes.labelsize": 10,
    "figure.dpi": 150, "savefig.dpi": 150, "savefig.bbox": "tight", "legend.frameon": False,
})
RESULTS: dict = {}


def log(msg: str) -> None:
    print(f"[{time.time() - T0:6.1f}s] {msg}", flush=True)


# ====================================================================================
# 1. Toy metabolic network (stoichiometric model)
# ====================================================================================
METS = ["glc", "g6p", "pyr", "lac", "atp", "adp", "nadph", "nadp", "x", "y"]


def build_reactions(effector: str = "NADPH"):
    """Reaction list: (name, stoichiometry, gene rule). Gene rules are illustrative."""
    detox = ({"x": -1, "nadph": -1, "y": 1, "nadp": 1} if effector == "NADPH"
             else {"x": -1, "atp": -1, "y": 1, "adp": 1})
    return [
        ("EX_glc", {"glc": 1}, "(boundary)"),
        ("HK", {"glc": -1, "atp": -1, "g6p": 1, "adp": 1}, "HK1"),
        ("GLYC", {"g6p": -1, "adp": -3, "pyr": 2, "atp": 3}, "PFK and lower glycolysis (lumped)"),
        ("PPP", {"g6p": -1, "nadp": -2, "adp": -8 / 3, "nadph": 2, "pyr": 5 / 3, "atp": 8 / 3},
         "G6PD (oxidative branch, lumped with recycling)"),
        ("LDH", {"pyr": -1, "lac": 1}, "LDH"),
        ("EX_lac", {"lac": -1}, "(boundary)"),
        ("ATPASE", {"atp": -1, "adp": 1}, "Na/K pump and other maintenance"),
        ("GSH_REDOX", {"nadph": -1, "nadp": 1}, "glutathione reductase (basal oxidative load)"),
        ("EX_x", {"x": 1}, "(boundary: toxin supply)"),
        ("DETOX", detox, "detoxE (transgene from the circuit)"),
        ("EX_y", {"y": -1}, "(boundary: product)"),
        ("EXPR_COST", {"atp": -1, "adp": 1}, "(ATP cost of making circuit proteins)"),
    ]


class ToyGEM:
    def __init__(self, effector: str = "NADPH"):
        rx = build_reactions(effector)
        self.effector = effector
        self.names = [r[0] for r in rx]
        self.gpr = {r[0]: r[2] for r in rx}
        self.stoich = {r[0]: r[1] for r in rx}
        self.idx = {n: i for i, n in enumerate(self.names)}
        self.n = len(rx)
        self.S = np.zeros((len(METS), self.n))
        for j, (_, st, _) in enumerate(rx):
            for met, c in st.items():
                self.S[METS.index(met), j] = c

    def bounds(self, **kw):
        p = {**BASE, **kw}
        i = self.idx
        lb, ub = np.zeros(self.n), np.full(self.n, 1000.0)
        ub[i["EX_glc"]] = p["glc_max"]
        if "glc_lb" in p:
            lb[i["EX_glc"]] = p["glc_lb"]
        if "glc_ub" in p:
            ub[i["EX_glc"]] = p["glc_ub"]
        ub[i["HK"]], ub[i["PPP"]] = p["hk_cap"], p["ppp_cap"]
        lb[i["ATPASE"]], lb[i["GSH_REDOX"]] = p["m_atp"], p["m_nadph"]
        ub[i["EX_x"]], ub[i["DETOX"]] = p["x_supply"], p["detox_ub"]
        lb[i["EXPR_COST"]] = ub[i["EXPR_COST"]] = p["expr_cost"]
        if "lac_lb" in p:
            lb[i["EX_lac"]] = p["lac_lb"]
        if "lac_ub" in p:
            ub[i["EX_lac"]] = p["lac_ub"]
        return lb, ub


# Invented parameters, in arbitrary flux units (young cell glucose uptake is about 0.4)
BASE = dict(glc_max=1.5, hk_cap=1.0, ppp_cap=0.40, m_atp=0.8, m_nadph=0.12,
            detox_ub=0.0, expr_cost=0.0, x_supply=1000.0)
HALF_LIFE = dict(G6PD=62.0, HK=100.0)  # days, enzyme activity declines with cell age


def caps_at(age_d, t_g: float = HALF_LIFE["G6PD"], t_hk: float = HALF_LIFE["HK"]):
    return dict(ppp_cap=BASE["ppp_cap"] * 2.0 ** (-age_d / t_g),
                hk_cap=BASE["hk_cap"] * 2.0 ** (-age_d / t_hk))


# ---- linear programming helpers -------------------------------------------------------
def _lp(m, c, lb, ub, A_ub=None, b_ub=None):
    return linprog(c, A_ub=A_ub, b_ub=b_ub, A_eq=m.S, b_eq=np.zeros(m.S.shape[0]),
                   bounds=np.column_stack([lb, ub]), method="highs")


def max_flux(m, lb, ub, name, A_ub=None, b_ub=None):
    c = np.zeros(m.n)
    c[m.idx[name]] = -1.0
    r = _lp(m, c, lb, ub, A_ub, b_ub)
    return (-r.fun, r.x) if r.status == 0 else (np.nan, None)


def min_total_flux(m, lb, ub):
    """Parsimonious solution: smallest total flux that meets all bounds."""
    r = _lp(m, np.ones(m.n), lb, ub)
    return r.x if r.status == 0 else None


def fva(m, lb, ub, names, A_ub=None, b_ub=None):
    out = {}
    for nme in names:
        c = np.zeros(m.n)
        c[m.idx[nme]] = 1.0
        lo = _lp(m, c, lb, ub, A_ub, b_ub)
        hi = _lp(m, -c, lb, ub, A_ub, b_ub)
        out[nme] = (lo.fun if lo.status == 0 else np.nan, -hi.fun if hi.status == 0 else np.nan)
    return out


def metabolic_max(m, **p):
    """Largest DETOX flux the cell can deliver while meeting maintenance demands."""
    lb, ub = m.bounds(**{**p, "detox_ub": 1000.0})
    return max_flux(m, lb, ub, "DETOX")[0]


def oxidative_reserve(m, v_detox, **p):
    """Fold increase of the basal glutathione-reductase flux still possible at this DETOX flux."""
    lb, ub = m.bounds(**{**p, "detox_ub": v_detox})
    lb[m.idx["DETOX"]] = v_detox
    g, _ = max_flux(m, lb, ub, "GSH_REDOX")
    return g / BASE["m_nadph"]


def atp_spare(m, v_detox, **p):
    lb, ub = m.bounds(**{**p, "detox_ub": v_detox})
    lb[m.idx["DETOX"]] = v_detox
    a, _ = max_flux(m, lb, ub, "ATPASE")
    return a - BASE["m_atp"]


def marginal_values(m, p, eps=1e-3):
    """Gain in maximal DETOX flux per unit relaxation of each constraint (a shadow price)."""
    base = metabolic_max(m, **p)
    tests = {"G6PD capacity (PPP)": ("ppp_cap", +1), "hexokinase capacity": ("hk_cap", +1),
             "glucose uptake limit": ("glc_max", +1), "ATP supply (less pump demand)": ("m_atp", -1),
             "NADPH supply (less basal demand)": ("m_nadph", -1)}
    out = {}
    for label, (key, sgn) in tests.items():
        q = dict(p)
        q[key] = {**BASE, **p}[key] + sgn * eps
        out[label] = (metabolic_max(m, **q) - base) / eps
    return base, out


# ---- objective free flux sampling (hit and run in the null space) ---------------------
def sample_fluxes(m, lb, ub, A_ub=None, b_ub=None, n=3000, thin=20, burn=500, seed=0):
    rng = np.random.default_rng(seed)
    free = (ub - lb) > 1e-9
    Sf, b = m.S[:, free], -m.S[:, ~free] @ lb[~free]
    nf = int(free.sum())
    if A_ub is not None:
        Af, bf = A_ub[:, free], b_ub - A_ub[:, ~free] @ lb[~free]
    else:
        Af, bf = np.zeros((0, nf)), np.zeros(0)
    lbf, ubf = lb[free], ub[free]
    pts = []
    for _ in range(6 * nf):
        r = linprog(rng.normal(size=nf), A_ub=Af if len(bf) else None, b_ub=bf if len(bf) else None,
                    A_eq=Sf, b_eq=b, bounds=np.column_stack([lbf, ubf]), method="highs")
        if r.status == 0:
            pts.append(r.x)
    if not pts:
        raise RuntimeError("empty feasible region")
    x = np.mean(pts, axis=0)
    N = null_space(Sf)
    out = []
    for it in range(burn + n * thin):
        d = N @ rng.normal(size=N.shape[1])
        d /= np.linalg.norm(d)
        with np.errstate(divide="ignore", invalid="ignore"):
            t1, t2 = (lbf - x) / d, (ubf - x) / d
        mk = np.abs(d) > 1e-12
        tlo, thi = np.max(np.minimum(t1, t2)[mk]), np.min(np.maximum(t1, t2)[mk])
        if len(bf):
            Ad, slack = Af @ d, bf - Af @ x
            pos, neg = Ad > 1e-12, Ad < -1e-12
            if pos.any():
                thi = min(thi, np.min(slack[pos] / Ad[pos]))
            if neg.any():
                tlo = max(tlo, np.max(slack[neg] / Ad[neg]))
        if thi - tlo < 1e-12:
            continue
        x = x + rng.uniform(tlo, thi) * d
        if it >= burn and (it - burn) % thin == 0:
            out.append(x.copy())
    X = np.zeros((len(out), m.n))
    X[:, free] = np.array(out)
    X[:, ~free] = lb[~free]
    return X


# ====================================================================================
# 2. Circuit model (ODE) and its coupling to the network
# ====================================================================================
# sensor (toxin level X) -> regulator R -> effector enzyme E; E removes X using NADPH
CIRC = dict(a_R=0.2, d_R=0.1, K_R=1.0, Km=1.0, L=0.20, kc=0.02, cbR=1.0, cbE=10.0,
            dE_A=np.log(2) / 24.0, T_exp=144.0, T_mat=2880.0, X_target=4.0)
TRUE_THETA = dict(a_E=0.02, kcat=1.0, K_X=2.0, hl=150.0)  # hl: effector half-life (days) in mature cell
NEED = CIRC["L"] - CIRC["kc"] * CIRC["X_target"]  # DETOX flux that keeps toxin at its target
RESERVE_MIN = 2.5                                  # required fold of spare redox capacity
DT = 0.5                                           # hours
PER_DAY = int(round(24 / DT))
REC_A = int(round(6 / DT))  # record every 6 h during the expression stage
N_STEPS = int((CIRC["T_exp"] + CIRC["T_mat"]) / DT)
DAYS_A = int(CIRC["T_exp"] / 24)                   # days of expression (erythroblast stage)
N_DAYS = DAYS_A + int(CIRC["T_mat"] / 24)          # recorded days: 0 .. N_DAYS


def hill(u, K, n=2):
    return u ** n / (K ** n + u ** n)


def run_coupled_lp(m, theta, t_g=HALF_LIFE["G6PD"]):
    """Level 2: operator splitting. At every step the circuit proposes a flux bound, the LP returns
    what the network can deliver, and the delivered flux goes back into the ODEs."""
    c = CIRC
    X, R, E = c["L"] / c["kc"], 0.0, 0.0
    rec = {k: [] for k in ["t_d", "X", "R", "E", "kin", "v", "vmet", "reserve", "atp_spare", "cost", "crisis"]}
    crisis_any = False
    for k in range(N_STEPS):
        t = k * DT
        on = t < c["T_exp"]
        hX, hR = hill(X, theta["K_X"]), hill(R, c["K_R"])
        cost = (c["cbR"] * c["a_R"] * hX + c["cbE"] * theta["a_E"] * hR) if on else 0.0
        age = 0.0 if on else (t - c["T_exp"]) / 24
        kin = theta["kcat"] * E * X / (c["Km"] + X)
        p = dict(caps_at(age, t_g=t_g), expr_cost=cost)
        lb, ub = m.bounds(**{**p, "detox_ub": kin})
        v, _ = max_flux(m, lb, ub, "DETOX")
        crisis = bool(np.isnan(v))
        crisis_any |= crisis
        v = 0.0 if crisis else v
        if (on and k % REC_A == 0) or (not on and k % PER_DAY == 0):
            vmet = np.nan if crisis else metabolic_max(m, **p)
            rec["t_d"].append(t / 24)
            for key, val in zip(["X", "R", "E", "kin", "v", "vmet", "cost", "crisis"],
                                [X, R, E, kin, v, vmet, cost, crisis]):
                rec[key].append(val)
            rec["reserve"].append(np.nan if crisis else oxidative_reserve(m, v, **p))
            rec["atp_spare"].append(np.nan if crisis else atp_spare(m, v, **p))
        v = min(v, max(X / DT + c["L"] - c["kc"] * X, 0.0))
        dE = c["dE_A"] if on else np.log(2) / (theta["hl"] * 24)
        X = max(X + DT * (c["L"] - v - c["kc"] * X), 1e-6)
        R += DT * ((c["a_R"] * hX if on else 0.0) - c["d_R"] * R)
        E += DT * ((theta["a_E"] * hR if on else 0.0) - dE * E)
    out = {k: np.array(v) for k, v in rec.items()}
    out["crisis_any"] = crisis_any
    return out


def simulate_table(theta, tabA, tabB, dt=DT, n_steps=None):
    """Same dynamics for many designs at once, using tabulated exact LP results:
    delivered flux = min(enzyme bound, network capacity), which is exact because the set of
    feasible DETOX fluxes is an interval [0, network capacity]."""
    c = CIRC
    a_E, kcat, K_X, hl = (np.atleast_1d(theta[k]).astype(float) for k in ("a_E", "kcat", "K_X", "hl"))
    N = max(len(a_E), len(kcat), len(K_X), len(hl))
    a_E, kcat, K_X, hl = (np.broadcast_to(v, (N,)).copy() for v in (a_E, kcat, K_X, hl))
    cost_grid, VA, cost_max = tabA
    age_grid, VB = tabB
    steps = n_steps or int((c["T_exp"] + c["T_mat"]) / dt)
    per_day = int(round(24 / dt))
    X = np.full(N, c["L"] / c["kc"])
    R = np.zeros(N)
    E = np.zeros(N)
    crisis = np.zeros(N, bool)
    cost_seen = np.zeros(N)
    days = steps // per_day + 1
    Xd, Ed, Vd, Md = (np.zeros((N, days)) for _ in range(4))
    dE_mat = np.log(2) / (hl * 24)
    for k in range(steps):
        t = k * dt
        on = t < c["T_exp"]
        hX, hR = hill(X, K_X), hill(R, c["K_R"])
        if on:
            cost = c["cbR"] * c["a_R"] * hX + c["cbE"] * a_E * hR
            crisis |= cost > cost_max
            cost_seen = np.maximum(cost_seen, cost)
            vmet = np.interp(cost, cost_grid, VA)
        else:
            cost = 0.0
            vmet = np.full(N, np.interp((t - c["T_exp"]) / 24, age_grid, VB))
        kin = kcat * E * X / (c["Km"] + X)
        v = np.minimum(kin, vmet)
        if k % per_day == 0:
            d = k // per_day
            Xd[:, d], Ed[:, d], Vd[:, d], Md[:, d] = X, E, v, vmet
        v = np.minimum(v, np.maximum(X / dt + c["L"] - c["kc"] * X, 0.0))
        X = np.maximum(X + dt * (c["L"] - v - c["kc"] * X), 1e-6)
        R = R + dt * ((c["a_R"] * hX if on else 0.0) - c["d_R"] * R)
        E = E + dt * ((a_E * hR if on else 0.0) - (c["dE_A"] if on else dE_mat) * E)
    return dict(X=Xd, E=Ed, v=Vd, vmet=Md, crisis=crisis, cost_max=cost_seen)


def protected_fraction(sim):
    """Fraction of mature life (age 0..119 d) with toxin at target AND redox reserve intact."""
    X, v, vmet = (sim[k][:, DAYS_A:DAYS_A + 120] for k in ("X", "v", "vmet"))
    reserve = 1 + (vmet - v) / BASE["m_nadph"]  # exact for this network: NADPH is the only shared pool
    ok = (X <= CIRC["X_target"]) & (reserve >= RESERVE_MIN)
    frac = ok.mean(axis=1)
    frac[sim["crisis"]] = 0.0
    return frac, (X <= CIRC["X_target"]).mean(axis=1), (reserve >= RESERVE_MIN).mean(axis=1)


# ====================================================================================
# 3. Build everything
# ====================================================================================
m = ToyGEM("NADPH")
mA = ToyGEM("ATP")
rng = np.random.default_rng(SEED)

# ---- self checks of the network ------------------------------------------------------
log("checking the toy network")
lb, ub = m.bounds(**caps_at(0.0))
v0 = min_total_flux(m, lb, ub)
RESULTS["checks"] = {
    "adenylate_pool_conserved": bool(np.allclose(m.S[METS.index("atp")] + m.S[METS.index("adp")], 0)),
    "nadp_pool_conserved": bool(np.allclose(m.S[METS.index("nadph")] + m.S[METS.index("nadp")], 0)),
    "baseline_steady_state_residual": float(np.abs(m.S @ v0).max()),
    "rank_S": int(np.linalg.matrix_rank(m.S)),
    "n_reactions": m.n, "n_metabolites": len(METS),
}
RESULTS["network"] = {"metabolites": METS,
                      "reactions": [{"name": n, "stoichiometry": m.stoich[n], "gene_rule": m.gpr[n]} for n in m.names]}

# ====================================================================================
# 4. LEVEL 1: stoichiometric coupling
# ====================================================================================
log("level 1: capacity, marginal values, operating window")
ages_tab = [0, 60, 120]
L1 = {"NADPH": {}, "ATP": {}}
for label, mod in (("NADPH", m), ("ATP", mA)):
    for a in ages_tab:
        p = caps_at(a)
        vmax, mv = marginal_values(mod, p)
        L1[label][a] = {"network_capacity": vmax, "marginal": mv,
                        "limiting": [k for k, val in sorted(mv.items(), key=lambda kv: -kv[1]) if val > 1e-3]}
RESULTS["level1_capacity"] = {lab: {str(a): d for a, d in dd.items()} for lab, dd in L1.items()}

# Figure 1: realised DETOX flux against the enzyme bound, for two kinds of effector
kgrid = np.linspace(0, 1.4, 141)
fig, axs = plt.subplots(1, 2, figsize=(9.6, 3.6), sharey=True)
cols = {0: TEAL, 60: "#E39A2D", 120: RED}
for ax, (label, mod, title) in zip(axs, (("NADPH", m, "Effector that uses NADPH"), ("ATP", mA, "Effector that uses ATP"))):
    for a in ages_tab:
        V = L1[label][a]["network_capacity"]
        ax.plot(kgrid, np.minimum(kgrid, V), color=cols[a], lw=2.2, label=f"cell age {a} d")
        ax.plot(kgrid[kgrid >= V][:1], [V], "o", color=cols[a], ms=4.5)
    top = L1[label][0]["limiting"][0].split(" (")[0]
    ax.set_title(f"{title}\nplateau set by: {top}", fontsize=10)
    ax.axhline(NEED, color=GREY, ls="--", lw=1)
    ax.text(1.38, NEED + 0.015, "flux needed to hold the toxin at its target", color=GREY, fontsize=8, ha="right")
    ax.set_xlabel("flux the enzyme could carry (kcat × E)")
    ax.set_xlim(0, 1.4)
axs[0].set_ylabel("flux the cell can actually deliver")
axs[0].legend(loc="upper left", fontsize=8)
axs[0].set_ylim(0, 1.35)
fig.tight_layout()
fig.savefig(FIG / "fig1_capacity.png")
plt.close(fig)

# Operating window for the NADPH effector: enzyme level x cell age
ages_w = np.linspace(0, 120, 61)
k_w = np.linspace(0, 0.78, 79)
cat = np.zeros((len(k_w), len(ages_w)), int)
res_fold = np.zeros_like(cat, float)
Vw = np.zeros(len(ages_w))
for j, a in enumerate(ages_w):
    p = caps_at(a)
    Vw[j] = metabolic_max(m, **p)
    for i, kk in enumerate(k_w):
        v = min(kk, Vw[j] - 1e-9)
        fold = oxidative_reserve(m, v, **p)
        res_fold[i, j] = fold
        cat[i, j] = 0 if v < NEED - 1e-9 else (1 if fold >= RESERVE_MIN else 2)


def window_closing_age():
    lo, hi = 0.0, 150.0
    for _ in range(30):
        mid = 0.5 * (lo + hi)
        ok = metabolic_max(m, **caps_at(mid)) >= NEED + (RESERVE_MIN - 1) * BASE["m_nadph"]
        lo, hi = (mid, hi) if ok else (lo, mid)
    return 0.5 * (lo + hi)


close_age = window_closing_age()
upper_safe = Vw - (RESERVE_MIN - 1) * BASE["m_nadph"]  # largest detox flux that keeps the reserve
grid_edge = np.array([k_w[cat[:, j] == 1].max() if (cat[:, j] == 1).any() else np.nan for j in range(len(ages_w))])
ok = np.isfinite(grid_edge)
RESULTS["checks"]["window_edge_grid_vs_formula_max_abs_diff"] = float(np.abs(grid_edge[ok] - upper_safe[ok]).max())
RESULTS["level1_window"] = {"need_flux": NEED, "reserve_min_fold": RESERVE_MIN, "closing_age_days": close_age,
                            "network_capacity_at_age": dict(zip([int(a) for a in ages_w[::6]], [float(x) for x in Vw[::6]])),
                            "series": {"age_d": [float(a) for a in ages_w], "network_capacity": [float(x) for x in Vw],
                                       "upper_safe": [float(x) for x in upper_safe]}}
RESULTS["level1_capacity_series"] = {"age_d": [float(a) for a in ages_w],
                                     "NADPH": [float(marginal_values(m, caps_at(a))[0]) for a in ages_w],
                                     "ATP": [float(marginal_values(mA, caps_at(a))[0]) for a in ages_w]}
fig, ax = plt.subplots(figsize=(6.6, 4.0))
cmap = ListedColormap([GREY_L, TEAL_L, RED_L])
ax.pcolormesh(ages_w, k_w, cat, cmap=cmap, shading="nearest", vmin=0, vmax=2)
ax.plot(ages_w, Vw, color=INK, lw=1.6)
ax.axvline(close_age, color=RED, ls=":", lw=1.5)
ax.text(close_age + 1.5, 0.77, f"window closes\nat about {close_age:.0f} days", color=RED, fontsize=8.5, va="top")
ax.text(3, 0.06, "too little detox", fontsize=9, color=INK)
ax.text(3, 0.30, "operating window", fontsize=9, color=TEAL, fontweight="bold")
ax.text(3, 0.735, "too much: little redox reserve left", fontsize=9, color=RED)
ax.annotate("network capacity", xy=(88, float(np.interp(88, ages_w, Vw))), xytext=(86, 0.47), fontsize=8.5, color=INK,
            arrowprops=dict(arrowstyle="-", color=INK, lw=0.8), ha="left")
ax.set_xlabel("cell age (days)")
ax.set_ylabel("flux the enzyme could carry (kcat × E)")
ax.set_title("Safe operating window of the effector")
fig.tight_layout()
fig.savefig(FIG / "fig2_window.png")
plt.close(fig)

# ====================================================================================
# 5. Synthetic data and constraining the network (level 1, steps 3)
# ====================================================================================
log("synthetic data: exchange rates, tracer fraction, enzyme activities")


def rep_lognormal(true, cv, n):
    sd = np.sqrt(np.log(1 + cv ** 2))
    return true * np.exp(rng.normal(-0.5 * sd ** 2, sd, size=n))


def ci95(x):
    h = student_t.ppf(0.975, len(x) - 1) * x.std(ddof=1) / np.sqrt(len(x))
    return x.mean() - h, x.mean() + h


# D1 exchange rates and D2 tracer PPP fraction: baseline cells (no circuit) at several ages
ages_d = np.array([0, 20, 40, 60, 80, 100, 120])
D1 = {"age": ages_d.tolist(), "glc": [], "lac": [], "ppp_frac": [], "truth": []}
for a in ages_d:
    lb, ub = m.bounds(**caps_at(a))
    v = min_total_flux(m, lb, ub)
    tr = dict(glc=v[m.idx["EX_glc"]], lac=v[m.idx["EX_lac"]], ppp_frac=v[m.idx["PPP"]] / v[m.idx["HK"]])
    D1["truth"].append(tr)
    D1["glc"].append(rep_lognormal(tr["glc"], 0.06, 4).tolist())
    D1["lac"].append(rep_lognormal(tr["lac"], 0.06, 4).tolist())
    D1["ppp_frac"].append(rep_lognormal(tr["ppp_frac"], 0.15, 4).tolist())

# D3 enzyme activity assays (relative to young cells), ages as above, 3 replicates, CV 10%
D3 = {}
for enz, t_half in (("G6PD", HALF_LIFE["G6PD"]), ("HK", HALF_LIFE["HK"])):
    D3[enz] = [rep_lognormal(2.0 ** (-a / t_half), 0.10, 3).tolist() for a in ages_d]


def fit_log2_decay(ages, acts):
    xs = np.repeat(ages, [len(a) for a in acts]).astype(float)
    ys = np.log2(np.concatenate(acts))
    Xm = np.column_stack([np.ones_like(xs), xs])
    beta = np.linalg.lstsq(Xm, ys, rcond=None)[0]
    resid = ys - Xm @ beta
    cov = resid @ resid / (len(xs) - 2) * np.linalg.inv(Xm.T @ Xm)
    return beta, cov


fitG, covG = fit_log2_decay(ages_d, D3["G6PD"])
fitH, covH = fit_log2_decay(ages_d, D3["HK"])
RESULTS["data_fits"] = {
    "G6PD_half_life_d": float(-1 / fitG[1]), "G6PD_half_life_true": HALF_LIFE["G6PD"],
    "G6PD_half_life_95": sorted(float(-1 / (fitG[1] + s * 1.96 * np.sqrt(covG[1, 1]))) for s in (+1, -1)),
    "HK_half_life_d": float(-1 / fitH[1]), "HK_half_life_true": HALF_LIFE["HK"],
}

# Flux sampling at age 60 d: exchange data only, then exchange data plus tracer fraction
i60 = list(ages_d).index(60)
glc_ci, lac_ci = ci95(np.array(D1["glc"][i60])), ci95(np.array(D1["lac"][i60]))
f_ci = ci95(np.array(D1["ppp_frac"][i60]))
open_caps = dict(ppp_cap=1.0, hk_cap=1.5, glc_max=1.5)
lbS, ubS = m.bounds(**open_caps, glc_lb=glc_ci[0], glc_ub=glc_ci[1], lac_lb=lac_ci[0], lac_ub=lac_ci[1])
A_tr = np.zeros((2, m.n))
A_tr[0, m.idx["PPP"]], A_tr[0, m.idx["HK"]] = 1.0, -f_ci[1]
A_tr[1, m.idx["PPP"]], A_tr[1, m.idx["HK"]] = -1.0, f_ci[0]
b_tr = np.zeros(2)
log("flux sampling")
smp_ex = sample_fluxes(m, lbS, ubS, seed=11)
smp_tr = sample_fluxes(m, lbS, ubS, A_tr, b_tr, seed=12)
smp_ex2 = sample_fluxes(m, lbS, ubS, seed=13)  # second chain to compare
ppp_frac = lambda X: X[:, m.idx["PPP"]] / X[:, m.idx["HK"]]
true_f60 = D1["truth"][i60]["ppp_frac"]
q = lambda x: [float(np.percentile(x, p)) for p in (2.5, 50, 97.5)]
RESULTS["sampling"] = {
    "age_d": 60, "n_samples": int(len(smp_ex)), "true_ppp_fraction": true_f60,
    "exchange_only": q(ppp_frac(smp_ex)), "with_tracer": q(ppp_frac(smp_tr)),
    "fva_ppp_flux_exchange_only": fva(m, lbS, ubS, ["PPP"])["PPP"],
    "fva_ppp_flux_with_tracer": fva(m, lbS, ubS, ["PPP"], A_tr, b_tr)["PPP"],
    "max_steady_state_violation": float(np.abs(smp_tr @ m.S.T).max()),
    "between_chain_mean_diff_over_sd": float(abs(ppp_frac(smp_ex).mean() - ppp_frac(smp_ex2).mean()) / ppp_frac(smp_ex).std()),
}

# Predicted network capacity for DETOX by cell age: knowledge from exchange and tracer data only, versus plus enzyme activity
log("tabulating exact LP network capacity over (G6PD capacity, cell age)")
ppp_axis = np.linspace(0.03, 1.2, 48)
age_axis = np.linspace(0, 120, 41)
Vgrid = np.array([[metabolic_max(m, ppp_cap=pc, hk_cap=caps_at(a)["hk_cap"]) for a in age_axis] for pc in ppp_axis])
V2D = RegularGridInterpolator((ppp_axis, age_axis), Vgrid, bounds_error=False, fill_value=np.nan)
age_fine = np.linspace(0, 125, 501)


def capacity_curve(t_g=None, ppp0=None):
    ppp = ppp0 * 2.0 ** (-age_fine / t_g)
    return V2D(np.column_stack([np.clip(ppp, ppp_axis[0], ppp_axis[-1]), np.clip(age_fine, 0, 120)]))


truth_curve = capacity_curve(HALF_LIFE["G6PD"], BASE["ppp_cap"])
# prior (exchange and tracer data only): capacity must at least cover the basal flux at every age
basal = BASE["m_nadph"] / 2
pri = []
while len(pri) < 400:
    p0, tg = 10 ** rng.uniform(np.log10(0.15), np.log10(1.0)), 10 ** rng.uniform(np.log10(20), np.log10(200))
    if p0 * 2.0 ** (-120 / tg) >= basal:
        pri.append(capacity_curve(tg, p0))
pri = np.array(pri)
post = []
for b in rng.multivariate_normal(fitG, covG, size=400):
    p0, tg = 2.0 ** b[0] * BASE["ppp_cap"], -1 / b[1]
    post.append(capacity_curve(tg, p0))
post = np.array(post)
band = lambda A: (np.nanpercentile(A, 2.5, axis=0), np.nanpercentile(A, 97.5, axis=0))
idx60 = np.argmin(abs(age_fine - 60))
RESULTS["prediction_bands"] = {
    "age_60_truth": float(truth_curve[idx60]),
    "age_60_without_enzyme_assay_95": [float(band(pri)[0][idx60]), float(band(pri)[1][idx60])],
    "age_60_with_enzyme_assay_95": [float(band(post)[0][idx60]), float(band(post)[1][idx60])],
    "truth_inside_post_band_fraction": float(np.mean(
        ((truth_curve >= band(post)[0]) & (truth_curve <= band(post)[1]))[age_fine <= 120])),
}

fig, axs = plt.subplots(1, 3, figsize=(12.2, 3.6))
ax = axs[0]
tt = np.array([d["glc"] for d in D1["truth"]])
for key, col, lab in (("glc", TEAL, "glucose uptake"), ("lac", RED, "lactate release")):
    data = np.array(D1[key])
    ax.plot(ages_d, [d[key] for d in D1["truth"]], color=col, lw=1, alpha=0.7)
    ax.errorbar(ages_d, data.mean(1), yerr=[data.mean(1) - np.array([ci95(r)[0] for r in data]), np.array([ci95(r)[1] for r in data]) - data.mean(1)],
                fmt="o", color=col, ms=4, capsize=2, label=lab)
ax.set_xlabel("cell age (days)")
ax.set_ylabel("exchange rate (flux units)")
ax.set_title("Synthetic exchange data (4 replicates)")
ax.legend(fontsize=8, loc="lower right")
ax.set_ylim(0, 1.0)
ax = axs[1]
bins = np.linspace(0, 1.0, 41)
ax.hist(ppp_frac(smp_ex), bins=bins, color=GREY_L, label="exchange rates only", density=True)
ax.hist(ppp_frac(smp_tr), bins=bins, color=TEAL, alpha=0.85, label="plus tracer fraction", density=True)
ax.axvline(true_f60, color=RED, lw=1.5, label="true value")
ax.set_xlabel("fraction of glucose through the PPP")
ax.set_ylabel("density of feasible flux states")
ax.set_title("Flux sampling at cell age 60 d")
ax.legend(fontsize=8)
ax = axs[2]
lo, hi = band(pri)
ax.fill_between(age_fine, lo, hi, color=GREY_L, ec="#B9C0CA", lw=0.6, label="exchange and tracer data only")
lo, hi = band(post)
ax.fill_between(age_fine, lo, hi, color=TEAL, alpha=0.6, label="plus enzyme activity assay")
ax.plot(age_fine, truth_curve, color=RED, lw=1.5, label="true capacity")
ax.axhline(NEED, color=GREY, ls="--", lw=1)
ax.set_xlim(0, 120)
ax.set_ylim(0, 0.8)
ax.set_xlabel("cell age (days)")
ax.set_ylabel("predicted network capacity")
ax.set_title("What each dataset tells us")
ax.legend(fontsize=8, loc="upper right")
fig.tight_layout()
fig.savefig(FIG / "fig3_data.png")
plt.close(fig)

# ====================================================================================
# 6. LEVEL 2: dynamic coupling
# ====================================================================================
log("level 2: coupled run with an LP at every step")
t_lp = time.time()
run = run_coupled_lp(m, TRUE_THETA)
RESULTS["level2_timing_s"] = time.time() - t_lp
age_rec = run["t_d"] - DAYS_A
d0 = int(np.argmax(age_rec >= -1e-9))  # first record at the start of the mature stage


def first_cross(mask, ages):
    idx = np.where(mask)[0]
    return float(ages[idx[0]]) if len(idx) else None


mature = slice(d0, d0 + 120)
tox_fail = first_cross(run["X"][mature] > CIRC["X_target"], age_rec[mature])
res_fail = first_cross(run["reserve"][mature] < RESERVE_MIN, age_rec[mature])
met_limited = run["v"] < run["kin"] - 1e-7  # the network, not the enzyme, sets the delivered flux
RESULTS["level2_truth_design"] = {
    "theta": TRUE_THETA, "x_at_end_of_expression": float(run["X"][d0]), "x_without_circuit": CIRC["L"] / CIRC["kc"],
    "toxin_loses_control_at_age_d": tox_fail, "reserve_below_threshold_at_age_d": res_fail,
    "metabolism_becomes_limiting_at_age_d": first_cross(met_limited[mature] & (run["kin"][mature] > 1e-6), age_rec[mature]),
    "max_expression_cost": float(run["cost"][:d0 + 1].max()),
    "min_atp_spare_during_expression": float(np.nanmin(run["atp_spare"][:d0 + 1])),
    "crisis_during_run": bool(run["crisis_any"]),
}

tdays = run["t_d"] - DAYS_A
mA, mB = tdays <= 1e-9, tdays >= -1e-9
fig = plt.figure(figsize=(9.2, 9.0))
gs = fig.add_gridspec(4, 2, width_ratios=[1.0, 5.2], wspace=0.05, hspace=0.2)
axl = [fig.add_subplot(gs[r, 0]) for r in range(4)]
axr = [fig.add_subplot(gs[r, 1], sharey=axl[r]) for r in range(4)]
for r in range(4):
    axl[r].set_xlim(-DAYS_A, 0)
    axr[r].set_xlim(0, 120)
    axl[r].patch.set_facecolor("#FBEEF0")
    axr[r].spines["left"].set_visible(False)
    axr[r].tick_params(left=False, labelleft=False)
    axl[r].set_xticks([-6, -3, 0])
    if r < 3:
        axl[r].tick_params(labelbottom=False)
        axr[r].tick_params(labelbottom=False)


def both(r, y, **kw):
    axl[r].plot(tdays[mA], y[mA], **kw)
    axr[r].plot(tdays[mB], y[mB], **kw)


axl[0].set_title("expression\nstage", fontsize=9, color=RED)
axr[0].set_title("mature cell: no new protein", fontsize=9)
both(0, run["R"], color="#7E5BA8", lw=1.8)
axl[0].set_ylabel("regulator R")
both(1, TRUE_THETA["kcat"] * run["E"], color=TEAL, ls="--", lw=1.6, label="enzyme bound (kcat × E)")
both(1, run["vmet"], color=INK, lw=1.8, label="network capacity")
both(1, run["v"], color=TEAL, lw=2.4, label="flux delivered")
for ax in (axl[1], axr[1]):
    ax.axhline(NEED, color=GREY, ls=":", lw=1.2)
axr[1].text(2, NEED - 0.025, "flux needed", color=GREY, fontsize=8, va="top")
axr[1].fill_between(tdays[mB], 0, run["v"][mB], where=(met_limited & (run["kin"] > 1e-6))[mB], color=RED_L, lw=0,
                    label="network is the limit")
axl[1].set_ylabel("DETOX flux")
axr[1].legend(fontsize=8, loc="upper right")
both(2, run["X"], color=INK, lw=2, label="toxin level, with circuit")
for ax in (axl[2], axr[2]):
    ax.axhline(CIRC["L"] / CIRC["kc"], color=GREY, ls=":", lw=1.4)
    ax.axhline(CIRC["X_target"], color=RED, ls="--", lw=1)
axl[2].set_yscale("log")
axr[2].set_yscale("log")
axl[2].set_ylim(1, 14)
axl[2].set_yticks([1, 2, 4, 10])
axl[2].set_yticklabels(["1", "2", "4", "10"])
axl[2].minorticks_off()
axr[2].text(2, CIRC["L"] / CIRC["kc"] * 0.88, "without circuit", color=GREY, fontsize=8, va="top")
axr[2].text(2, CIRC["X_target"] * 1.05, "target", color=RED, fontsize=8, va="bottom")
axl[2].set_ylabel("toxin level (log scale)")
both(3, run["reserve"], color=RED, lw=2)
for ax in (axl[3], axr[3]):
    ax.axhline(RESERVE_MIN, color=RED, ls="--", lw=1)
axr[3].text(118, RESERVE_MIN + 0.12, "minimum reserve", color=RED, fontsize=8, ha="right")
axl[3].set_ylabel("spare redox capacity\n(fold over basal)")
axl[3].set_xlabel("days")
axr[3].set_xlabel("cell age (days since the expression stage)")
fig.savefig(FIG / "fig4_dynamics.png")
plt.close(fig)

# ====================================================================================
# 7. Tables of exact LP results, used for fast search over many designs
# ====================================================================================
log("tabulating exact LP results for the design search")
cost_grid = np.arange(0, 2.2, 0.02)
VA = np.array([metabolic_max(m, **{**caps_at(0), "expr_cost": cst}) for cst in cost_grid])
ok = ~np.isnan(VA)
cost_max = float(cost_grid[ok][-1])
tabA = (cost_grid[ok], VA[ok], cost_max)
SCEN = {"fast G6PD decline (half-life 40 d)": 40.0, "baseline (62 d)": 62.0, "slow decline (90 d)": 90.0}
age_day = np.linspace(0, 125, 501)


def tabB_for(t_g):
    return age_day, V2D(np.column_stack([np.clip(BASE["ppp_cap"] * 2.0 ** (-age_day / t_g), ppp_axis[0], ppp_axis[-1]), np.clip(age_day, 0, 120)]))


# consistency with the exact LP at awkward points
chk_ages = np.array([0, 17, 33, 61, 95, 119.0])
direct = np.array([metabolic_max(m, **caps_at(a)) for a in chk_ages])
interp = np.interp(chk_ages, *tabB_for(HALF_LIFE["G6PD"]))
RESULTS["checks"]["table_vs_exact_lp_max_abs_diff"] = float(np.abs(direct - interp).max())
costs_chk = np.array([0.0, 0.3, 0.7, 1.1])
RESULTS["checks"]["cost_table_vs_exact_lp_max_abs_diff"] = float(np.abs(
    [metabolic_max(m, **{**caps_at(0), "expr_cost": c_}) for c_ in costs_chk] - np.interp(costs_chk, tabA[0], tabA[1])).max())
RESULTS["checks"]["max_expression_cost_feasible"] = cost_max

# the same coupled run with the table, at two time steps, to check LP coupling and time step
sim_tab = simulate_table(TRUE_THETA, tabA, tabB_for(HALF_LIFE["G6PD"]))
sim_fine = simulate_table(TRUE_THETA, tabA, tabB_for(HALF_LIFE["G6PD"]), dt=0.125)
_daily = np.isclose(run["t_d"] % 1.0, 0.0)
_xa, _xb = run["X"][_daily], sim_tab["X"][0, np.round(run["t_d"][_daily]).astype(int)]
RESULTS["checks"]["lp_run_vs_table_run_max_rel_diff_X"] = float(np.max(np.abs(_xa - _xb) / _xb))
RESULTS["checks"]["dt_0.5_vs_0.125_max_rel_diff_X"] = float(np.max(np.abs(sim_tab["X"][0, :-1] - sim_fine["X"][0, :-1]) / sim_fine["X"][0, :-1]))
# redox reserve from the closed form 1 + (capacity - v)/m against the LP value
reserve_cf = 1 + (run["vmet"] - run["v"]) / BASE["m_nadph"]
RESULTS["checks"]["reserve_closed_form_vs_lp_max_abs_diff"] = float(np.nanmax(np.abs(reserve_cf - run["reserve"])))
log("checks: " + json.dumps({k: (round(v, 6) if isinstance(v, float) else v) for k, v in RESULTS["checks"].items()}))

# ====================================================================================
# 8. Design ensemble and identifiability with synthetic time course data
# ====================================================================================
log("ensemble of circuit designs")
LOG_RANGES = dict(a_E=(np.log10(0.004), np.log10(0.25)), kcat=(np.log10(0.2), np.log10(5.0)),
                  K_X=(np.log10(0.5), np.log10(8.0)), hl=(np.log10(30.0), np.log10(400.0)))
sob = qmc.Sobol(d=4, scramble=True, seed=SEED).random_base2(14)
THETA = {k: 10 ** (LOG_RANGES[k][0] + sob[:, i] * (LOG_RANGES[k][1] - LOG_RANGES[k][0]))
         for i, k in enumerate(["a_E", "kcat", "K_X", "hl"])}
sims = {name: simulate_table(THETA, tabA, tabB_for(tg)) for name, tg in SCEN.items()}
prot = {name: protected_fraction(s) for name, s in sims.items()}
robust = np.min([p[0] for p in prot.values()], axis=0)
daily_age = np.arange(120)
RESULTS["ensemble"] = {
    "n_designs": int(len(robust)),
    "designs_failing_from_expression_cost": float(np.mean(sims["baseline (62 d)"]["crisis"])),
    "best_protected_fraction_per_scenario": {n: float(p[0].max()) for n, p in prot.items()},
    "network_bound_protected_fraction": {n: float(np.mean(
        np.interp(daily_age, *tabB_for(tg)) >= NEED + (RESERVE_MIN - 1) * BASE["m_nadph"])) for n, tg in SCEN.items()},
    "robust_best": float(robust.max()), "robust_median": float(np.median(robust)),
    "fraction_of_designs_within_5pct_of_best_robust": float(np.mean(robust >= 0.95 * robust.max())),
}

# synthetic time course of the toxin: truth design, baseline scenario, 3 replicates, 12% lognormal noise
obs_days = np.array([3, 5] + [DAYS_A + a for a in (5, 15, 30, 45, 60, 75, 90, 105, 119)])
Xtrue = sim_tab["X"][0, obs_days]
D4 = np.array([rep_lognormal(x, 0.12, 3) for x in Xtrue])
logD = np.log(D4).mean(1)
dist = np.sqrt(np.mean((np.log(np.maximum(sims["baseline (62 d)"]["X"][:, obs_days], 1e-9)) - logD) ** 2, axis=1))
dist[sims["baseline (62 d)"]["crisis"]] = np.inf
acc = dist <= 0.15
if acc.sum() < 60:
    acc = dist <= np.sort(dist)[60]
log(f"accepted designs: {int(acc.sum())} of {len(acc)}")
lg = {k: np.log10(THETA[k]) for k in THETA}
prod = lg["a_E"] + lg["kcat"]
iqr = lambda x: float(np.subtract(*np.percentile(x, [75, 25])))
contraction = {k: 1 - iqr(lg[k][acc]) / iqr(lg[k]) for k in lg}
contraction["a_E × kcat"] = 1 - iqr(prod[acc]) / iqr(prod)
# add one more measurement: effector amount at the end of the expression stage (15% noise)
E_true = sim_tab["E"][0, DAYS_A]
E_obs = rep_lognormal(E_true, 0.15, 3)
dist2 = np.sqrt((dist ** 2 * len(obs_days) + (np.log(np.maximum(sims["baseline (62 d)"]["E"][:, DAYS_A], 1e-12)) - np.log(E_obs).mean()) ** 2) / (len(obs_days) + 1))
acc2 = (dist2 <= 0.15) & np.isfinite(dist)
if acc2.sum() < 25:
    acc2 = dist2 <= np.sort(dist2[np.isfinite(dist2)])[25]
contraction2 = {k: 1 - iqr(lg[k][acc2]) / iqr(lg[k]) for k in lg}
contraction2["a_E × kcat"] = 1 - iqr(prod[acc2]) / iqr(prod)
RESULTS["identifiability"] = {
    "n_accepted": int(acc.sum()), "n_accepted_with_enzyme_amount": int(acc2.sum()),
    "corr_log_aE_log_kcat_accepted": float(np.corrcoef(lg["a_E"][acc], lg["kcat"][acc])[0, 1]),
    "corr_log_aE_log_kcat_with_enzyme_amount": float(np.corrcoef(lg["a_E"][acc2], lg["kcat"][acc2])[0, 1]),
    "contraction_toxin_data": contraction, "contraction_plus_enzyme_amount": contraction2,
    "log10_ranges": {k: [float(LOG_RANGES[k][0]), float(LOG_RANGES[k][1])] for k in ("a_E", "kcat")},
    "accepted_log10": {k: [float(x) for x in lg[k][acc]] for k in ("a_E", "kcat")},
    "accepted_with_enzyme_amount_log10": {k: [float(x) for x in lg[k][acc2]] for k in ("a_E", "kcat")},
    "true_log10": {k: float(np.log10(v)) for k, v in TRUE_THETA.items()},
    "protected_fraction_of_accepted_designs_baseline": [float(x) for x in np.percentile(prot["baseline (62 d)"][0][acc], [5, 50, 95])],
    "robust_protected_fraction_of_accepted_designs": [float(x) for x in np.percentile(robust[acc], [5, 50, 95])],
    "protected_fraction_truth_design": float(protected_fraction(sim_tab)[0][0]),
    "expression_burden_truth_design": float(sim_tab["cost_max"][0]),
    "expression_burden_accepted_designs_5_50_95": [float(x) for x in np.percentile(sims["baseline (62 d)"]["cost_max"][acc], [5, 50, 95])],
    "expression_burden_accepted_with_enzyme_amount_5_50_95": [float(x) for x in np.percentile(sims["baseline (62 d)"]["cost_max"][acc2], [5, 50, 95])],
}

fig, axs = plt.subplots(1, 3, figsize=(12.6, 3.9))
ax = axs[0]
ax.scatter(lg["a_E"], lg["kcat"], s=2, color=GREY_L, label="all designs")
ax.scatter(lg["a_E"][acc], lg["kcat"][acc], s=6, color=TEAL, label="fit the toxin time course")
ax.scatter(lg["a_E"][acc2], lg["kcat"][acc2], s=9, color="#E39A2D", label="and the enzyme amount")
ax.plot(np.log10(TRUE_THETA["a_E"]), np.log10(TRUE_THETA["kcat"]), "*", color=RED, ms=13, label="true design")
ax.set_xlabel("log10 expression strength a_E")
ax.set_ylabel("log10 kcat")
ax.set_title("Expression and kcat compensate")
ax.legend(fontsize=7.5, loc="lower left")
ax = axs[1]
names = ["a_E", "kcat", "a_E × kcat", "K_X", "hl"]
labels = ["expression\na_E", "kcat", "a_E × kcat\n(product)", "sensor\nthreshold", "enzyme\nhalf-life"]
xx = np.arange(len(names))
ax.bar(xx - 0.2, [contraction[n] for n in names], 0.38, color=TEAL, label="toxin time course")
ax.bar(xx + 0.2, [contraction2[n] for n in names], 0.38, color="#E39A2D", label="plus enzyme amount")
ax.set_xticks(xx)
ax.set_xticklabels(labels, fontsize=8)
ax.set_ylabel("narrowing of the plausible range")
ax.set_ylim(0, 1.22)
ax.set_title("What the data fix")
ax.legend(fontsize=8, ncol=2, loc="upper center")
ax = axs[2]
valid = ~sims["baseline (62 d)"]["crisis"]
hb = ax.hexbin(prod[valid], lg["hl"][valid], C=robust[valid], gridsize=28, cmap="viridis", reduce_C_function=np.mean, mincnt=3)
ax.plot(np.log10(TRUE_THETA["a_E"] * TRUE_THETA["kcat"]), np.log10(TRUE_THETA["hl"]), "*", color=RED, ms=12)
cb = fig.colorbar(hb, ax=ax)
cb.set_label("protected fraction of life,\nworst of three G6PD scenarios")
ax.set_xlabel("log10 (a_E × kcat)")
ax.set_ylabel("log10 effector half-life (d)")
ax.set_title("Design score is capped by the network")
fig.tight_layout()
fig.savefig(FIG / "fig5_designs.png")
plt.close(fig)

# ====================================================================================
# 9. Save results
# ====================================================================================
def clean(o):
    if isinstance(o, dict):
        return {str(k): clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [clean(v) for v in o]
    if isinstance(o, (np.floating, float)):
        return None if not np.isfinite(o) else float(o)
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.bool_,)):
        return bool(o)
    return o


RESULTS["circuit_parameters"] = {**{k: float(v) for k, v in CIRC.items()}, "true_design": TRUE_THETA,
                                 "need_flux": NEED, "reserve_min_fold": RESERVE_MIN, "dt_hours": DT}
RESULTS["base_parameters"] = {**BASE, "half_lives_days": HALF_LIFE}
RESULTS["synthetic_data"] = {"ages_d": ages_d.tolist(), "D1_exchange": D1, "D3_activity": D3,
                             "D4_toxin_days": obs_days.tolist(), "D4_toxin_replicates": D4.tolist(),
                             "E_amount_obs_replicates": E_obs.tolist()}
RESULTS["level2_series"] = {"age_d": age_rec.tolist(), **{k: run[k].tolist() for k in ("X", "R", "E", "kin", "v", "vmet", "reserve", "atp_spare")}}
(OUT / "results.json").write_text(json.dumps(clean(RESULTS), indent=1))
log("done: results/results.json and results/figures written")
