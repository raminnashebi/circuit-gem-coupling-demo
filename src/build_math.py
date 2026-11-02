#!/usr/bin/env python3
"""Build docs/mathematics.html: the equations and derivations behind coupling_demo.py.

Needs results/results.json (run src/coupling_demo.py first) and `pip install latex2mathml`. The equations are written
in LaTeX here and converted to MathML, so the finished page is one self contained HTML file with no
JavaScript and no internet needed. Every closed form in the text is recomputed below and checked against
a linear program or against the demo's own results, and the checks are printed in section 10.
"""
import html
import json
import re
from pathlib import Path

import numpy as np
from latex2mathml.converter import convert
from scipy.optimize import brentq, linprog

ROOT = Path(__file__).resolve().parents[1]
R = json.loads((ROOT / "results" / "results.json").read_text())
BP, CP = R["base_parameters"], R["circuit_parameters"]
TRUE = CP["true_design"]

# ---------------------------------------------------------------- parameters (all read from the demo)
g0, h0, mA, mN, uglc = BP["ppp_cap"], BP["hk_cap"], BP["m_atp"], BP["m_nadph"], BP["glc_max"]
tG, tH = BP["half_lives_days"]["G6PD"], BP["half_lives_days"]["HK"]
aR, dR, KR, Km, L, kc = CP["a_R"], CP["d_R"], CP["K_R"], CP["Km"], CP["L"], CP["kc"]
cR, cE, dEA, Texp, Tmat, Xbar = CP["cbR"], CP["cbE"], CP["dE_A"], CP["T_exp"], CP["T_mat"], CP["X_target"]
rho, dt = CP["reserve_min_fold"], CP["dt_hours"]
nu = L - kc * Xbar
Theta = nu + (rho - 1) * mN
aE, kcat, KX, hE = TRUE["a_E"], TRUE["kcat"], TRUE["K_X"], TRUE["hl"]

# ---------------------------------------------------------------- independent re-implementation
mets, rxs = R["network"]["metabolites"], R["network"]["reactions"]
names = [r["name"] for r in rxs]
idx = {n: i for i, n in enumerate(names)}


def build_S(effector="NADPH"):
    S = np.zeros((len(mets), len(rxs)))
    for j, r in enumerate(rxs):
        st = dict(r["stoichiometry"])
        if r["name"] == "DETOX" and effector == "ATP":
            st = {"x": -1, "atp": -1, "y": 1, "adp": 1}
        for m, c in st.items():
            S[mets.index(m), j] = c
    return S


def V_lp(a, c, eff="NADPH"):
    S = build_S(eff)
    lb, ub = np.zeros(S.shape[1]), np.full(S.shape[1], 1000.0)
    ub[idx["EX_glc"]], ub[idx["HK"]], ub[idx["PPP"]] = uglc, h0 * 2 ** (-a / tH), g0 * 2 ** (-a / tG)
    lb[idx["ATPASE"]], lb[idx["GSH_REDOX"]] = mA, mN
    lb[idx["EXPR_COST"]] = ub[idx["EXPR_COST"]] = c
    obj = np.zeros(S.shape[1])
    obj[idx["DETOX"]] = -1
    r = linprog(obj, A_eq=S, b_eq=np.zeros(S.shape[0]), bounds=np.column_stack([lb, ub]), method="highs")
    return -r.fun if r.status == 0 else np.nan


def Hcap(a):
    return min(h0 * 2 ** (-a / tH), uglc)


def gcap(a, tg=tG):
    return g0 * 2 ** (-a / tg)


def V_cf(a, c, eff="NADPH", tg=tG):
    H, g = Hcap(a), gcap(a, tg)
    if eff == "NADPH":
        p = min(g, H, 3 * (2 * H - mA - c))
        return 2 * p - mN if p >= mN / 2 else np.nan
    v = 2 * H - mN / 6 - mA - c
    return v if v >= 0 else np.nan


def hill(u, K):
    return u * u / (K * K + u * u)


def sim_cf(a_E, k, K_X, h_E, tg=tG):
    """Coupled run written only from the equations of sections 5 and 6 (closed form network capacity, Euler)."""
    X, Rr, E = L / kc, 0.0, 0.0
    per_day, n_steps = int(round(24 / dt)), int((Texp + Tmat) / dt)
    Xd, Ed = [], []
    for n in range(n_steps):
        t = n * dt
        on = t < Texp
        hX, hR = hill(X, K_X), hill(Rr, KR)
        c = (cR * aR * hX + cE * a_E * hR) if on else 0.0
        a = 0.0 if on else (t - Texp) / 24
        V = V_cf(a, c, tg=tg)
        kappa = k * E * X / (Km + X)
        v = 0.0 if np.isnan(V) else min(kappa, V)
        if n % per_day == 0:
            Xd.append(X)
            Ed.append(E)
        v = min(v, max(X / dt + L - kc * X, 0.0))
        X = max(X + dt * (L - v - kc * X), 1e-6)
        Rr += dt * ((aR * hX if on else 0.0) - dR * Rr)
        E += dt * ((a_E * hR if on else 0.0) - (dEA if on else np.log(2) / (h_E * 24)) * E)
    return np.array(Xd), np.array(Ed)


# ---------------------------------------------------------------- numerical verification suite
CHECKS = []


def check(label, value, tol, shown, note=""):
    CHECKS.append((label, shown, value <= tol, note))


# (1) closed form network capacity against the linear program, both effectors
worst, nfeas = 0.0, 0
for eff in ("NADPH", "ATP"):
    for a in range(0, 121, 5):
        for c in (0.0, 0.3, 0.7, 1.0, 1.1, 1.17, 1.2):
            x, y = V_lp(a, c, eff), V_cf(a, c, eff)
            assert np.isnan(x) == np.isnan(y), ("feasibility mismatch", eff, a, c, x, y)
            if not np.isnan(x):
                worst, nfeas = max(worst, abs(x - y)), nfeas + 1
check("Network capacity closed forms against the linear program (both effectors, ages 0 to 120 d, costs 0 to 1.2)",
      worst, 1e-9, f"largest difference {worst:.1e} over {nfeas} feasible cases; feasibility agrees in all cases")
cmax_cf = 2 * Hcap(0) - mA - mN / 6
check("Largest feasible expression cost, closed form against the demo", abs(cmax_cf - R["checks"]["max_expression_cost_feasible"]), 1e-9,
      f"{cmax_cf:.4f} against {R['checks']['max_expression_cost_feasible']:.4f}")
a_close = tG * np.log2(2 * g0 / (Theta + mN))
a_nu = tG * np.log2(2 * g0 / (nu + mN))
check("Window closing age, closed form against the demo's bisection", abs(a_close - R["level1_window"]["closing_age_days"]), 1e-6,
      f"{a_close:.3f} d against {R['level1_window']['closing_age_days']:.3f} d")
# (2) baseline fluxes against the synthetic truth
D1 = R["synthetic_data"]["D1_exchange"]["truth"][0]
f_hk = (mA + mN / 6) / 2
p_ppp = mN / 2
lac = 2 * (f_hk - p_ppp) + 5 / 3 * p_ppp
check("Parsimonious baseline: glucose uptake, lactate release, PPP share against the demo's true values",
      max(abs(f_hk - D1["glc"]), abs(lac - D1["lac"]), abs(p_ppp / f_hk - D1["ppp_frac"])), 1e-9,
      f"{f_hk:.3f}, {lac:.3f}, {p_ppp / f_hk:.4f} against {D1['glc']:.3f}, {D1['lac']:.3f}, {D1['ppp_frac']:.4f}")
# (3) the whole coupled run re-implemented from the equations against the demo's linear program run
Xs, Es = sim_cf(aE, kcat, KX, hE)
ser = R["level2_series"]
age_s, X_s = np.array(ser["age_d"]), np.array(ser["X"])
daily = np.isclose(age_s % 1.0, 0.0)
di = np.round(age_s[daily] + 6).astype(int)
rel = np.max(np.abs(Xs[di] - X_s[daily]) / X_s[daily])
Xneg = sim_cf(aE, 1.01 * kcat, KX, hE)[0]
rel_neg = np.max(np.abs(Xneg[di] - X_s[daily]) / X_s[daily])
check("Coupled run re-implemented from the equations (closed form network capacity, Euler) against the demo's run (a linear program every step)",
      rel, 1e-4, f"largest relative difference in the toxin level {rel:.1e}; for comparison, changing kcat by 1% gives {rel_neg:.0e}")
# (4) structural non identifiability along a_E x kcat = constant
a2, a3 = 0.0137, 0.0531
Xa, Ea = sim_cf(a2, aE * kcat / a2, KX, hE)
Xb, Eb = sim_cf(a3, aE * kcat / a3, KX, hE)
dX = max(np.max(np.abs(Xa - Xs)), np.max(np.abs(Xb - Xs)))
ratio = Ea[6] / Es[6], Eb[6] / Es[6]
dX_neg = np.max(np.abs(sim_cf(a2, 1.01 * aE * kcat / a2, KX, hE)[0] - Xs))
check("Three designs with the same product a_E x kcat: toxin trajectories",
      dX, 1e-9, f"largest difference {dX:.0e}; effector amounts at day 6 differ by factors {ratio[0]:.2f} and {ratio[1]:.2f}; for comparison, a 1% change in the product gives {dX_neg:.0e}")
# (5) quasi steady state against the simulation (mature stage, ages 0 to 119)
E0 = Es[6]
ages = np.arange(120)
Eage = E0 * 2 ** (-ages / hE)
Vage = np.array([V_cf(a, 0.0) for a in ages])
B = L - kc * Km - kcat * Eage
Xk = (B + np.sqrt(B * B + 4 * kc * L * Km)) / (2 * kc)
Xq = np.maximum(Xk, (L - Vage) / kc)
Xsim = Xs[6:126]
rel_q = np.abs(Xq - Xsim) / Xsim
check("Quasi steady state formula against the simulated toxin level (ages 0 to 119 d)", rel_q.mean(), 0.05,
      f"mean relative error {rel_q.mean() * 100:.1f}%, largest {rel_q.max() * 100:.1f}% (at age {int(ages[rel_q.argmax()])} d, just after the network becomes the limit at {a_nu:.0f} d and the toxin starts to rise faster than the formula follows)")
vq = L - kc * Xq
Rq = 1 + (Vage - vq) / mN
Rsim = np.array(ser["reserve"])[np.isclose(age_s % 1, 0) & (age_s >= 0)][:120]
check("Redox reserve from the quasi steady state against the linear program value", np.abs(Rq - Rsim).max(), 0.05,
      f"largest absolute difference {np.abs(Rq - Rsim).max():.3f} (fold over basal)")
a_kappa = hE * np.log2(E0 * kcat * Xbar / (nu * (Km + Xbar)))
a_tox = min(a_kappa, a_nu)
first = lambda mask: int(ages[np.argmax(mask)]) if mask.any() else None
sim_tox, sim_res = R["level2_truth_design"]["toxin_loses_control_at_age_d"], R["level2_truth_design"]["reserve_below_threshold_at_age_d"]
a_res = brentq(lambda a: 1 + (V_cf(a, 0.0) - (L - kc * max(((L - kc * Km - kcat * E0 * 2 ** (-a / hE)) + np.sqrt((L - kc * Km - kcat * E0 * 2 ** (-a / hE)) ** 2 + 4 * kc * L * Km)) / (2 * kc), (L - V_cf(a, 0.0)) / kc))) / mN - rho, 0, 119)
check("Age at which the toxin loses control: formula against the simulation", abs(a_tox - sim_tox), 3.0,
      f"formula {a_tox:.1f} d (network limit {a_nu:.1f} d, enzyme limit {a_kappa:.1f} d); simulation {sim_tox:.0f} d on daily records")
check("Age at which the redox reserve falls below its threshold: formula against the simulation", abs(a_res - sim_res), 1.5,
      f"formula {a_res:.1f} d; simulation {sim_res:.0f} d on daily records")
# (6) upper bound on protected life
bound_an = {}
for tg in (40.0, 62.0, 90.0):
    bound_an[tg] = np.mean([V_cf(a, 0.0, tg=tg) >= Theta for a in ages])
bound_dm = list(R["ensemble"]["network_bound_protected_fraction"].values())
check("Upper bound on protected life, formula against the demo (three G6PD scenarios)", max(abs(np.array(list(bound_an.values())) - bound_dm)), 1e-9,
      "formula " + ", ".join(f"{v:.4f}" for v in bound_an.values()) + "; demo " + ", ".join(f"{v:.4f}" for v in bound_dm))
best = list(R["ensemble"]["best_protected_fraction_per_scenario"].values())
within = max(abs(np.array(best) - np.array(list(bound_an.values())))) * 120
check("Best design in 16,384 against the upper bound (in days of mature life)", within, 1.0 + 1e-9, f"largest gap {within:.1f} day")
# (7) expression stage fixed point
def Ffix(X):
    Rs = aR / dR * hill(X, KX)
    Es_ = aE / dEA * hill(Rs, KR)
    return kcat * Es_ * X / (Km + X) + kc * X - L


Xstar = brentq(Ffix, 1e-6, 50)
Rstar = aR / dR * hill(Xstar, KX)
Estar = aE / dEA * hill(Rstar, KR)
cstar = cR * aR * hill(Xstar, KX) + cE * aE * hill(Rstar, KR)
x6 = R["level2_truth_design"]["x_at_end_of_expression"]


def relax_expression(days=40):
    X, Rr, E = L / kc, 0.0, 0.0
    for n in range(int(days * 24 / dt)):
        hX, hR = hill(X, KX), hill(Rr, KR)
        V = V_cf(0.0, cR * aR * hX + cE * aE * hR)
        v = min(kcat * E * X / (Km + X), V)
        v = min(v, max(X / dt + L - kc * X, 0.0))
        X = max(X + dt * (L - v - kc * X), 1e-6)
        Rr += dt * (aR * hX - dR * Rr)
        E += dt * (aE * hR - dEA * E)
    return X


Xrelax = relax_expression()
check("Expression stage fixed point against the ODE with expression held on for 40 days", abs(Xstar - Xrelax) / Xstar, 1e-4,
      f"fixed point {Xstar:.4f}; ODE after 40 days {Xrelax:.4f}; the demo's six day stage ends at {x6:.2f}, still on the way there")

# ---------------------------------------------------------------- numbers used in the text
SD = R["identifiability"]
NUM = {
    "g0": f"{g0:g}", "h0": f"{h0:g}", "tG": f"{tG:g}", "tH": f"{tH:g}", "uglc": f"{uglc:g}", "mA": f"{mA:g}", "mN": f"{mN:g}",
    "aR": f"{aR:g}", "dR": f"{dR:g}", "KR": f"{KR:g}", "Km": f"{Km:g}", "L": f"{L:g}", "kc": f"{kc:g}", "cR": f"{cR:g}", "cE": f"{cE:g}",
    "Xbar": f"{Xbar:g}", "nu": f"{nu:.2f}", "rho": f"{rho:g}", "Theta": f"{Theta:.2f}", "Texp": f"{Texp:g}", "dt": f"{dt:g}",
    "aE": f"{aE:g}", "kcat": f"{kcat:g}", "KX": f"{KX:g}", "hE": f"{hE:g}", "cmax": f"{cmax_cf:.2f}", "cflat": f"{2 * Hcap(0) - mA - g0 / 3:.4f}",
    "Vyoung": f"{V_cf(0, 0):.2f}", "VATPyoung": f"{V_cf(0, 0, 'ATP'):.2f}", "a_close": f"{a_close:.1f}", "a_nu": f"{a_nu:.1f}",
    "gTheta": f"{(Theta + mN) / 2:.2f}", "f_hk": f"{f_hk:.2f}", "phi": f"{p_ppp / f_hk:.3f}", "lacb": f"{lac:.2f}",
    "a_kappa": f"{a_kappa:.1f}", "a_tox": f"{a_tox:.1f}", "a_res": f"{a_res:.1f}", "sim_tox": f"{sim_tox:.0f}", "sim_res": f"{sim_res:.0f}",
    "Xstar": f"{Xstar:.2f}", "fp_err": f"{abs(Xstar - Xrelax) / Xstar:.0e}", "cstar": f"{cstar:.2f}", "E0": f"{E0:.3f}",
    "bound40": f"{bound_an[40.0]:.3f}", "bound62": f"{bound_an[62.0]:.3f}", "bound90": f"{bound_an[90.0]:.3f}",
    "n_acc": f"{SD['n_accepted']}", "burden_lo": f"{SD['expression_burden_accepted_designs_5_50_95'][0]:.2f}",
    "burden_hi": f"{SD['expression_burden_accepted_designs_5_50_95'][2]:.2f}", "corr": f"{SD['corr_log_aE_log_kcat_accepted']:.2f}",
    "Nd": f"{R['ensemble']['n_designs']:,}", "g0half": f"{g0 / 2:g}",
}

# ---------------------------------------------------------------- the document (LaTeX between $ ... $ and $$ ... $$)
def matrix_tex():
    from fractions import Fraction
    S = build_S("NADPH")

    def cell(x):
        fr = Fraction(x).limit_denominator(3)
        if fr == 0:
            return "0"
        if fr.denominator == 1:
            return str(fr.numerator)
        return ("-" if fr < 0 else "") + r"\tfrac{" + str(abs(fr.numerator)) + "}{" + str(fr.denominator) + "}"
    rows = [" & ".join(cell(x) for x in row) for row in S]
    return r"\mathbf S=\begin{pmatrix}" + r" \\ ".join(rows) + r"\end{pmatrix}"


BODY = r"""
<h1>Mathematical formulation of the coupling demonstration</h1>
<p class="sub">Equations, derivations and numerical checks for the toy model of a gene circuit coupled to a metabolic network. Companion to the demonstration report, in the same notation as the code.</p>
<div class="banner"><strong>Scope.</strong> These are exact statements about the <em>toy</em> model, whose network and parameters are invented. They show the structure of the method. They are not claims about real red cells.</div>

<div class="box"><h3>What this page establishes</h3>
<ul>
<li><strong>The network side is a linear program with an exact answer.</strong> The set of feasible effector fluxes is an interval $[0,V]$ (Lemma 1), so the cell delivers exactly $\min\{\kappa,V\}$ whatever the circuit asks for. For this network $V$ has a closed form (Proposition 1).</li>
<li><strong>The circuit side is three ODEs</strong> that meet the network only through that flux map.</li>
<li><strong>In the mature stage a quasi steady state formula</strong> predicts the toxin level, and it implies an upper bound on the protected life of <em>any</em> circuit design (Propositions 2 and 3).</li>
<li><strong>Expression strength and kcat are structurally non identifiable</strong> from toxin data, because they enter only as a product (Proposition 4).</li>
<li><strong>Every formula is recomputed</strong> against a linear program or against the demo's own output in section 10.</li>
</ul></div>

<h2>1. Notation and parameters</h2>
<table>
<tr><th>Symbol</th><th>Meaning</th><th>Value or unit</th></tr>
<tr><td>$t$, $a$</td><td>time; cell age, $a(t)=\max\{0,\,t-T_{\mathrm{exp}}\}/24$</td><td>hours; days. $T_{\mathrm{exp}}={{Texp}}$ h (6 d)</td></tr>
<tr><td>$v\in\mathbf R^{12}$, $S\in\mathbf R^{10\times 12}$</td><td>flux vector; stoichiometric matrix</td><td>flux units (young glucose uptake is about 0.4)</td></tr>
<tr><td>$g(a)$, $h(a)$</td><td>G6PD capacity (upper bound of the PPP flux); hexokinase capacity</td><td>$g_0 2^{-a/\tau_G}$, $h_0 2^{-a/\tau_H}$ with $g_0={{g0}}$, $\tau_G={{tG}}$ d, $h_0={{h0}}$, $\tau_H={{tH}}$ d</td></tr>
<tr><td>$u_{\mathrm{glc}}$, $H(a)$</td><td>glucose uptake limit; effective uptake capacity $\min\{h(a),u_{\mathrm{glc}}\}$</td><td>$u_{\mathrm{glc}}={{uglc}}$</td></tr>
<tr><td>$m_A$, $m_N$</td><td>maintenance: pump ATP; basal NADPH demand (glutathione reductase)</td><td>${{mA}}$; ${{mN}}$</td></tr>
<tr><td>$c$, $b$</td><td>ATP drawn by expression; flux bound on the effector reaction set by the circuit</td><td>$c\ge 0$; $b\ge0$</td></tr>
<tr><td>$V(a,c)$</td><td>network capacity: largest feasible effector flux</td><td>section 3</td></tr>
<tr><td>$E$, $R$, $X$</td><td>effector enzyme; regulator; toxin level</td><td>arbitrary units</td></tr>
<tr><td>$\kappa(E,X)$</td><td>kinetic bound $k\,E\,X/(K_m+X)$</td><td>$K_m={{Km}}$</td></tr>
<tr><td>$L$, $k_c$, $\bar X$</td><td>toxin load; clearance by other routes; target toxin level</td><td>${{L}}$; ${{kc}}$ h$^{-1}$; ${{Xbar}}$</td></tr>
<tr><td>$\eta$, $\rho$, $\Theta$</td><td>detox flux needed, $\eta=L-k_c\bar X$; required redox reserve (fold over basal); $\Theta=\eta+(\rho-1)m_N$</td><td>${{nu}}$; ${{rho}}$; ${{Theta}}$</td></tr>
<tr><td>$a_R,d_R,K_R$</td><td>regulator synthesis, decay, threshold</td><td>${{aR}}$, ${{dR}}$ h$^{-1}$, ${{KR}}$</td></tr>
<tr><td>$c_R,c_E$</td><td>ATP cost per unit of regulator, effector synthesis</td><td>${{cR}}$, ${{cE}}$</td></tr>
<tr><td>$\theta=(a_E,k,K_X,h_E)$</td><td>design: effector synthesis rate, kcat, sensor threshold, effector half life</td><td>true: ${{aE}}$, ${{kcat}}$, ${{KX}}$, ${{hE}}$ d</td></tr>
<tr><td>$\delta_A$, $\delta_M$</td><td>effector loss rate in the expression stage, in the mature stage</td><td>$\ln 2/24$ h$^{-1}$; $\ln 2/(24h_E)$ h$^{-1}$</td></tr>
</table>

<h2>2. The network</h2>
<h3>2.1 Reactions and balances</h3>
<p>Twelve reactions. With $\mathrm{atp}$/$\mathrm{adp}$ and $\mathrm{nadph}$/$\mathrm{nadp}$ as conserved pairs, ten metabolites give eight independent balances. Writing $\delta=1$ for an ATP using effector and $\delta=0$ for an NADPH using effector, the steady state $Sv=0$ reads</p>
$$\begin{array}{rcl}
v_{\mathrm{ex,glc}}-v_{\mathrm{HK}} &=& 0 \quad(\mathrm{glc})\\
v_{\mathrm{HK}}-v_{\mathrm{GLYC}}-v_{\mathrm{PPP}} &=& 0 \quad(\mathrm{g6p})\\
2v_{\mathrm{GLYC}}+\tfrac53 v_{\mathrm{PPP}}-v_{\mathrm{LDH}} &=& 0 \quad(\mathrm{pyr})\\
v_{\mathrm{LDH}}-v_{\mathrm{ex,lac}} &=& 0 \quad(\mathrm{lac})\\
-v_{\mathrm{HK}}+3v_{\mathrm{GLYC}}+\tfrac83 v_{\mathrm{PPP}}-v_{\mathrm{ATPase}}-v_{\mathrm{expr}}-\delta\,v_{\mathrm{DET}} &=& 0 \quad(\mathrm{atp})\\
2v_{\mathrm{PPP}}-v_{\mathrm{GSH}}-(1-\delta)\,v_{\mathrm{DET}} &=& 0 \quad(\mathrm{nadph})\\
v_{\mathrm{ex,x}}-v_{\mathrm{DET}} &=& 0 \quad(\mathrm{toxin})\\
v_{\mathrm{DET}}-v_{\mathrm{ex,y}} &=& 0 \quad(\mathrm{product})
\end{array}$$ @@balance
<p>Two consequences are used throughout. Eliminating $v_{\mathrm{GLYC}}$ with the hexokinase flux $f=v_{\mathrm{HK}}$ and the PPP flux $p=v_{\mathrm{PPP}}$ gives</p>
$$v_{\mathrm{GLYC}}=f-p\ge 0,\qquad v_{\mathrm{GSH}}+(1-\delta)\,v_{\mathrm{DET}}=2p,\qquad 2f-\tfrac13 p=v_{\mathrm{ATPase}}+v_{\mathrm{expr}}+\delta\,v_{\mathrm{DET}}.$$ @@reduced
<p>That is, NADPH is made only by the PPP, and ATP is made by glycolysis at two per glucose minus a third of the PPP flux. The full matrix is in the box below.</p>
<details><summary>Show the matrix $S$ (rows: glc, g6p, pyr, lac, atp, adp, nadph, nadp, x, y; columns: ex,glc; HK; GLYC; PPP; LDH; ex,lac; ATPase; GSH; ex,x; DET; ex,y; expr)</summary>
$$ {{MATRIX}} $$
</details>

<h3>2.2 Bounds, ageing and maintenance</h3>
<table>
<tr><th>Flux</th><th>Lower bound</th><th>Upper bound</th></tr>
<tr><td>$v_{\mathrm{ex,glc}}$</td><td>$0$</td><td>$u_{\mathrm{glc}}$</td></tr>
<tr><td>$v_{\mathrm{HK}}$</td><td>$0$</td><td>$h(a)$</td></tr>
<tr><td>$v_{\mathrm{PPP}}$</td><td>$0$</td><td>$g(a)$</td></tr>
<tr><td>$v_{\mathrm{ATPase}}$ (pump)</td><td>$m_A$</td><td>$\bar u$</td></tr>
<tr><td>$v_{\mathrm{GSH}}$ (basal redox)</td><td>$m_N$</td><td>$\bar u$</td></tr>
<tr><td>$v_{\mathrm{DET}}$ (effector)</td><td>$0$</td><td>$b$</td></tr>
<tr><td>$v_{\mathrm{expr}}$ (expression)</td><td>$c$</td><td>$c$</td></tr>
<tr><td>all others</td><td>$0$</td><td>$\bar u$ (1000, never binding)</td></tr>
</table>
<p>The cell does not grow, so there is no biomass objective. Maintenance enters as lower bounds on $v_{\mathrm{ATPase}}$ and $v_{\mathrm{GSH}}$, and the feasible set for an age $a$, an expression cost $c$ and an effector bound $b$ is</p>
$$\Pi(a,c,b)=\{\,v\in\mathbf R^{12}:\ Sv=0,\ \ \ell(a,c)\le v\le u(a,c,b)\,\},\qquad g(a)=g_0 2^{-a/\tau_G},\ \ h(a)=h_0 2^{-a/\tau_H}.$$ @@polytope

<h2>3. Level 1: what the network can deliver</h2>
<h3>3.1 Network capacity and delivered flux</h3>
<p>The network capacity is the largest effector flux the cell can support while meeting maintenance (the maximum that flux variability analysis returns for that reaction), and the delivered flux is the largest one when the circuit offers at most $b$:</p>
$$V(a,c)=\max\{\,v_{\mathrm{DET}}:\ v\in\Pi(a,c,\bar u)\,\},\qquad v^\star(a,c,b)=\max\{\,v_{\mathrm{DET}}:\ v\in\Pi(a,c,b)\,\}.$$ @@capacity
<div class="prop"><strong>Lemma 1 (interval).</strong> Let $c\le c_{\max}(a)$, so that $\Pi(a,c,\bar u)\ne\emptyset$. Then the set of feasible effector fluxes is the interval $[0,V(a,c)]$, and
$$v^\star(a,c,b)=\min\{\,b,\ V(a,c)\,\}.$$ @@lemma
If $c>c_{\max}(a)$ the problem is infeasible (an energy crisis).</div>
<div class="proof"><em>Proof.</em> $\Pi(a,c,\bar u)$ is a convex polytope, and its image under the coordinate map $v\mapsto v_{\mathrm{DET}}$ is a convex subset of the line, hence an interval. It contains $0$, because the maintenance solution with $v_{\mathrm{ex,x}}=v_{\mathrm{DET}}=v_{\mathrm{ex,y}}=0$ is feasible exactly when $\Pi\ne\emptyset$. Its upper end is $V$ by definition. If $b\ge V$ the bound $b$ is not active; if $b<V$ then $b$ lies in the interval and is attained. $\square$</div>
<p>This is what lets the demo replace a linear program at every time step by a table: the network capacity depends only on $(a,c)$, and the circuit enters through $\min\{b,V\}$.</p>

<h3>3.2 Closed forms</h3>
<div class="prop"><strong>Proposition 1 (network capacity).</strong> Let $H=H(a)$, $g=g(a)$. For the NADPH using effector,
$$V(a,c)=2\min\{\,g,\ H,\ 3(2H-m_A-c)\,\}-m_N,\qquad \text{feasible}\Longleftrightarrow c\le c_{\max}(a)=2H(a)-m_A-\tfrac16 m_N .$$ @@Vclosed
For the ATP using effector,
$$V^{\mathrm{ATP}}(a,c)=2H(a)-\tfrac16 m_N-m_A-c,\qquad \text{feasible}\Longleftrightarrow c\le c_{\max}(a).$$ @@VATP</div>
<div class="proof"><em>Proof (NADPH effector).</em> By @@reduced, $v_{\mathrm{DET}}=2p-v_{\mathrm{GSH}}\le 2p-m_N$, so $V$ is largest for the largest feasible $p$. The constraints on $p$ are $p\le g$ (G6PD capacity), $p\le f\le H$ (the PPP is fed by hexokinase), and the ATP balance $2f-\tfrac13p\ge m_A+c$. The left side increases with $f$, so take $f=H$; the ATP balance becomes $p\le 3(2H-m_A-c)$. Hence $p\le\bar p:=\min\{g,H,3(2H-m_A-c)\}$ and $V=2\bar p-m_N$. Feasibility needs $\bar p\ge m_N/2$ (so that $v_{\mathrm{GSH}}\ge m_N$ is possible); the binding term is the third, which gives $c\le 2H-m_A-m_N/6$. <em>ATP effector.</em> Now $v_{\mathrm{GSH}}=2p\ge m_N$, so the smallest $p=m_N/2$ is best, and the ATP balance gives $v_{\mathrm{DET}}=2H-\tfrac13p-m_A-c$ with $p=m_N/2$. $\square$</div>
<p>In the demo's parameters the first branch of the minimum, $g(a)$, is the active one for every age from 0 to 120 d when $c=0$ (the mature stage), and at age 0 for every cost up to ${{cflat}}$, which is $2H-m_A-g_0/3$ (the expression stage). There $V=2g(a)-m_N$, which is ${{Vyoung}}$ for a young cell, and the network capacity for the ATP effector is ${{VATPyoung}}$ at age 0. Table interpolation in the code is exact up to the kinks of this piecewise linear function.</p>

<h3>3.3 Marginal values (what one more unit of a resource buys)</h3>
<p>The marginal value of a model parameter $\vartheta$ is the derivative of the network capacity, evaluated in the demo by a forward difference $[V(\vartheta+\epsilon)-V(\vartheta)]/\epsilon$ with $\epsilon=10^{-3}$. From Proposition 1 the exact values in each regime are</p>
<table>
<tr><th>Effector and active term</th><th>$\partial V/\partial g$</th><th>$\partial V/\partial H$</th><th>$\partial V/\partial m_A$</th><th>$\partial V/\partial m_N$</th><th>$\partial V/\partial c$</th></tr>
<tr><td>NADPH, $\bar p=g$ (the demo)</td><td>$2$</td><td>$0$</td><td>$0$</td><td>$-1$</td><td>$0$</td></tr>
<tr><td>NADPH, $\bar p=H$</td><td>$0$</td><td>$2$</td><td>$0$</td><td>$-1$</td><td>$0$</td></tr>
<tr><td>NADPH, $\bar p=3(2H-m_A-c)$</td><td>$0$</td><td>$12$</td><td>$-6$</td><td>$-1$</td><td>$-6$</td></tr>
<tr><td>ATP</td><td>$0$</td><td>$2$</td><td>$-1$</td><td>$-\tfrac16$</td><td>$-1$</td></tr>
</table>
<p>The demo's values (2.0, 1.0 and 0.1667 for the corresponding resources) are these numbers, with the sign of the demand parameters read as "supply".  The marginal value is the shadow price of the constraint in the dual linear program, so the active term in the minimum names the limiting resource.</p>

<h3>3.4 Redox reserve and the operating window</h3>
<p>Fix the effector at flux $v\le V$ and ask how far the basal redox flux could be raised. Since $v_{\mathrm{GSH}}=2p-v$ and $p\le\bar p$ whatever $v$ is,</p>
$$\Lambda(a,c;v)=\frac{\max\{v_{\mathrm{GSH}}:\ v\in\Pi,\ v_{\mathrm{DET}}=v\}}{m_N}=\frac{2\bar p-v}{m_N}=1+\frac{V(a,c)-v}{m_N}.$$ @@reserve
<p>With the offered bound $b$ the delivered flux is $v=\min\{b,V\}$. The effector is useful if $v\ge\eta$ and safe if $\Lambda\ge\rho$. If $b>V$ the reserve is $1<\rho$, so only $b\le V$ qualifies, and then $V-b\ge(\rho-1)m_N$. Therefore</p>
$$W(a)=\left[\,\eta,\ \ V(a)-(\rho-1)m_N\,\right],\qquad W(a)\neq\emptyset\Longleftrightarrow V(a)\ge\Theta=\eta+(\rho-1)m_N.$$ @@window
<p>With $V=2g-m_N$ the window closes where $g(a^\ast)=(\Theta+m_N)/2$, and the effector cannot even reach $\eta$ after $a_\eta$ where $g(a_\eta)=(\eta+m_N)/2$:</p>
$$a^\ast=\tau_G\log_2\frac{2g_0}{\Theta+m_N}={{a_close}}\ \text{d},\qquad a_\eta=\tau_G\log_2\frac{2g_0}{\eta+m_N}={{a_nu}}\ \text{d}.$$ @@closing

<h2>4. Data and constraining the network</h2>
<h3>4.1 How the synthetic data are generated</h3>
<p>For a true value $\theta_i$ with coefficient of variation $\mathrm{CV}$, replicate $j$ is lognormal with mean $\theta_i$:</p>
$$y_{ij}=\theta_i\exp\!\left(\sigma\xi_{ij}-\tfrac12\sigma^2\right),\quad \xi_{ij}\sim\mathrm{N}(0,1),\quad \sigma^2=\ln(1+\mathrm{CV}^2),$$ @@lognormal
<p>with $n=4$, $\mathrm{CV}=6\%$ for exchange rates (D1), $n=4$, $15\%$ for the tracer fraction (D2), $n=3$, $10\%$ for enzyme activity (D3), $n=3$, $12\%$ for the toxin level (D4) and $n=3$, $15\%$ for the enzyme amount (D5). The 95% interval used as a flux bound is the $t$ interval $\bar y\pm t_{0.975,\,n-1}\,s/\sqrt n$.</p>
<p>The true baseline fluxes (no circuit) are the parsimonious solution $\min\mathbf 1^\top v$ subject to $v\in\Pi(a,0,0)$. With no effector, NADPH is used only by the basal reaction, so the minimal solution has</p>
$$v_{\mathrm{PPP}}=\frac{m_N}{2},\qquad v_{\mathrm{HK}}=\frac{m_A+m_N/6}{2}={{f_hk}},\qquad v_{\mathrm{ex,lac}}=2(v_{\mathrm{HK}}-v_{\mathrm{PPP}})+\tfrac53v_{\mathrm{PPP}}={{lacb}},$$ @@baseline
<p>and a PPP share $\phi=v_{\mathrm{PPP}}/v_{\mathrm{HK}}={{phi}}$, independent of age while hexokinase is not limiting. This is why the exchange data in the demo are flat.</p>

<h3>4.2 Constraining and sampling the flux space</h3>
<p>The exchange intervals $[\underline g,\bar g]$ and $[\underline\ell,\bar\ell]$ for glucose and lactate become bounds on $v_{\mathrm{ex,glc}}$ and $v_{\mathrm{ex,lac}}$. The tracer interval $[\underline\phi,\bar\phi]$ for the PPP share is a pair of linear inequalities, because $\phi=v_{\mathrm{PPP}}/v_{\mathrm{HK}}$ and $v_{\mathrm{HK}}\ge0$:</p>
$$v_{\mathrm{PPP}}-\bar\phi\,v_{\mathrm{HK}}\le 0,\qquad -v_{\mathrm{PPP}}+\underline\phi\,v_{\mathrm{HK}}\le 0.$$ @@tracer
<p>Rather than optimise, the demo samples the polytope $\Omega=\{v:\ Sv=0,\ \ell\le v\le u,\ Av\le0\}$ by hit and run. Let $x$ be the free coordinates (those with $u_j>\ell_j$), $S_f$ the corresponding columns and $Z$ an orthonormal basis of $\operatorname{ker}S_f$. From a point $x\in\Omega$ draw $w\sim\mathrm{N}(0,I)$, set $d=Zw/\|Zw\|$ and move to $x+td$ with $t$ uniform on</p>
$$[t_{\mathrm{lo}},t_{\mathrm{hi}}]=\left[\max_{i:\,d_i\ne0}\min\!\left\{\tfrac{\ell_i-x_i}{d_i},\tfrac{u_i-x_i}{d_i}\right\},\ \min_{i:\,d_i\ne0}\max\!\left\{\tfrac{\ell_i-x_i}{d_i},\tfrac{u_i-x_i}{d_i}\right\}\right]\ \cap\ \{t:\ A(x+td)\le 0\}.$$ @@hitrun
<p>The chain has the uniform distribution on $\Omega$ as its stationary law, so it describes the feasible flux states without choosing an objective. The starting point is the mean of random vertices (an interior point), with 500 burn in steps and a thinning of 20. Flux variability analysis is the companion problem $\min/\max\,v_j$ over $\Omega$.</p>

<h3>4.3 Enzyme activity and the predicted network capacity</h3>
<p>Activities relative to young cells follow $\log_2A_{ij}=\beta_0+\beta_1a_i+e_{ij}$. Least squares gives</p>
$$\hat\beta=(X^\top X)^{-1}X^\top y,\qquad \hat\Sigma=s^2(X^\top X)^{-1},\qquad \hat\tau=-1/\hat\beta_1 .$$ @@regression
<p>Draws $\beta^{(s)}\sim\mathrm{N}(\hat\beta,\hat\Sigma)$ give G6PD capacity curves $g^{(s)}(a)=g_02^{\beta^{(s)}_0+\beta^{(s)}_1a}$ and, through Proposition 1, network capacity curves $V\left(g^{(s)}(a),a\right)$ whose 2.5th and 97.5th percentiles form the band in the demo. Without the enzyme assay the same band is computed from a prior on $(g_0,\tau_G)$ that is log uniform on $[0.15,1]\times[20,200]$ d and conditioned on $g(120)\ge m_N/2$ (capacity must cover the basal redox need).</p>

<h2>5. The circuit model</h2>
<p>With the stage indicator $\chi(t)=1$ for $t<T_{\mathrm{exp}}$ and $0$ afterwards, and $h(u;K)=u^2/(K^2+u^2)$:</p>
$$\begin{array}{rcl}
\dfrac{dR}{dt} &=& \chi(t)\,a_R\,h(X;K_X)-d_R\,R, \\
\dfrac{dE}{dt} &=& \chi(t)\,a_E\,h(R;K_R)-\delta(t)\,E, \\
\dfrac{dX}{dt} &=& L-v_{\mathrm{DET}}-k_c\,X .
\end{array}$$ @@ode
<p>Here $\delta(t)=\delta_A$ for $t<T_{\mathrm{exp}}$ and $\delta_M$ afterwards. Expression costs ATP, drawn only in the expression stage, and the enzyme can carry at most the kinetic bound:</p>
$$c(t)=\chi(t)\left[c_R\,a_R\,h(X;K_X)+c_E\,a_E\,h(R;K_R)\right],\qquad \kappa(t)=\frac{k\,E\,X}{K_m+X}.$$ @@costkappa
<p><strong>Expression stage.</strong> Setting the derivatives to zero with $\chi=1$ and $v_{\mathrm{DET}}=\kappa$ (the network is not limiting) gives a fixed point $X^\ast$ that solves one scalar equation:</p>
$$\kappa\left(E^\ast(X^\ast),X^\ast\right)+k_cX^\ast=L,\qquad R^\ast(X)=\frac{a_R}{d_R}h(X;K_X),\quad E^\ast(X)=\frac{a_E}{\delta_A}h\left(R^\ast(X);K_R\right).$$ @@fixedpoint
<p>This is the self regulation of the sensor: more toxin raises $R$ and then $E$, which removes toxin. For the true design $X^\ast={{Xstar}}$ at an expression cost of ${{cstar}}$. Held on for 40 days the ODE converges to it (relative difference {{fp_err}}). The six day stage in the demo ends at {{x6_sim}}, still on the way, because the approach is a damped oscillation.</p>
<p><strong>Mature stage.</strong> With $\chi=0$ the regulator decays, $R(t)=R(T_{\mathrm{exp}})e^{-d_R(t-T_{\mathrm{exp}})}$, within about two days, and the effector follows $E(a)=E_0\,2^{-a/h_E}$ with $E_0={{E0}}$ for the true design. No new protein is made, so the circuit is open loop.</p>

<h2>6. The coupling</h2>
<p>At each time $t$ the circuit hands the network a bound $b=\kappa(t)$ and an expression cost $c(t)$, and the network returns the delivered flux. By Lemma 1</p>
$$v_{\mathrm{DET}}(t)=\begin{cases}\min\{\kappa(t),\ V(a(t),c(t))\} & c(t)\le c_{\max}\\ 0\ \ (\text{crisis, flagged}) & c(t)>c_{\max}.\end{cases}$$ @@fluxmap
<p>This closes the loop. The ODEs set $\kappa$ and $c$, and the network sets $v_{\mathrm{DET}}$, which drives the toxin level that the sensor reads. The scheme is operator splitting with forward Euler, $\Delta t={{dt}}$ h, and a guard that stops the toxin level going negative:</p>
$$v_n\leftarrow\min\left\{v_n,\ \max\{X_n/\Delta t+L-k_cX_n,\,0\}\right\},\qquad X_{n+1}=\max\{X_n+\Delta t\,(L-v_n-k_cX_n),\ 10^{-6}\}.$$ @@euler
<p>$R$ and $E$ are advanced with the same step. The demo solves the linear program at every step for the true design and uses a table of $V$ for the design search; Lemma 1 says the two are the same map, and section 10 confirms it.</p>

<h2>7. Quasi steady state of the mature stage</h2>
<p>The toxin relaxes on a time scale of at most $1/k_c=50$ h, while the enzyme and the network capacity change over days to months. Setting $dX/dt=0$ at fixed age gives $L=\min\{\kappa(E,X),V\}+k_cX$. Both $F_\kappa(X)=\kappa+k_cX$ and $F_V(X)=V+k_cX$ are strictly increasing in $X$, so $\min\{F_\kappa,F_V\}=L$ has a unique root, and it is the larger of the two single roots.</p>
<div class="prop"><strong>Proposition 2 (quasi steady state).</strong>
$$X_{\mathrm{qs}}(a)=\max\{\,X_\kappa(a),\ X_V(a)\,\},\qquad X_V=\frac{L-V(a)}{k_c},\qquad X_\kappa=\frac{B+\sqrt{B^2+4k_cLK_m}}{2k_c},\quad B=L-k_cK_m-kE(a).$$ @@qss
Consequently the toxin is at or below target if and only if the enzyme and the network can each deliver the needed flux at the target level:
$$X_{\mathrm{qs}}\le\bar X\Longleftrightarrow \kappa(E,\bar X)\ge\eta\ \ \text{and}\ \ V(a)\ge\eta .$$ @@qsscontrol</div>
<div class="proof"><em>Proof.</em> $X_\kappa$ solves $(L-k_cX)(K_m+X)=kEX$, that is $k_cX^2-BX-LK_m=0$, whose positive root is the formula. The equivalence follows from $\min\{F_\kappa,F_V\}(\bar X)\ge L$ and monotonicity. $\square$</div>
<p>The first condition fails at $a_\kappa=h_E\log_2\dfrac{E_0\,k\,\bar X}{\eta\,(K_m+\bar X)}={{a_kappa}}$ d, the second at $a_\eta={{a_nu}}$ d, so the toxin is lost at $a_{\mathrm{tox}}=\min\{a_\kappa,a_\eta\}={{a_tox}}$ d. The network, not the enzyme, is the first to fail. The simulation gives ${{sim_tox}}$ d on daily records. The reserve in the quasi steady state is $\Lambda=1+\left(V-(L-k_cX_{\mathrm{qs}})\right)/m_N$, which crosses $\rho$ at ${{a_res}}$ d (simulation: ${{sim_res}}$ d).</p>

<h2>8. An upper bound on every design</h2>
<p>A design is <em>protected</em> at age $a$ if the toxin is at target and the redox reserve is intact. The protected fraction of the 120 days of mature life is</p>
$$P(\theta)=\frac1{120}\sum_{a=0}^{119}\mathbf 1\left[\,X(a)\le\bar X\ \text{ and }\ \Lambda(a)\ge\rho\,\right].$$ @@protected
<div class="prop"><strong>Proposition 3 (upper bound on design performance).</strong> In the quasi steady state, protection at age $a$ implies $V(a)\ge\Theta$. Hence for every design $\theta$
$$P(\theta)\ \le\ \bar P=\frac1{120}\left|\,\{a\in\{0,\dots,119\}:\ V(a)\ge\Theta\}\,\right|,$$ @@boundP
and $\bar P$ depends only on the network and the G6PD decline, not on the circuit.</div>
<div class="proof"><em>Proof.</em> Protection gives $X_{\mathrm{qs}}\le\bar X$, so the delivered flux is $v=L-k_cX_{\mathrm{qs}}\ge L-k_c\bar X=\eta$. The reserve condition $1+(V-v)/m_N\ge\rho$ then gives $V\ge v+(\rho-1)m_N\ge\eta+(\rho-1)m_N=\Theta$. $\square$</div>
<p>With $V=2g-m_N$ the set $\{a:V(a)\ge\Theta\}$ is $\{0,\dots,\lfloor a^\ast\rfloor\}$, so $\bar P=(\lfloor a^\ast\rfloor+1)/120$ with $a^\ast=\tau_G\log_2\!\left(2g_0/(\Theta+m_N)\right)$. For $\tau_G=40,\,62,\,90$ d this gives ${{bound40}}$, ${{bound62}}$ and ${{bound90}}$. These are exactly the values the demo computes from its table of the network capacity, and the best of {{Nd}} designs in each scenario is within one day of them (a transient can overshoot the quasi steady state by a day).</p>

<h2>9. Design ensemble, inference and identifiability</h2>
<p><strong>Ensemble.</strong> The design vector $\theta=(a_E,k,K_X,h_E)$ is drawn from a scrambled Sobol sequence, log uniform on $a_E\in[0.004,0.25]$, $k\in[0.2,5]$, $K_X\in[0.5,8]$, $h_E\in[30,400]$ d ({{Nd}} points). The robust score is the worst case over three G6PD half lives, $P_{\mathrm{rob}}(\theta)=\min_{\tau_G\in\{40,62,90\}}P_{\tau_G}(\theta)$. Designs with $c(t)>c_{\max}$ at any time score $0$.</p>
<p><strong>Inference.</strong> With log toxin observations $\bar y_j$ at $n$ days, a design is accepted if</p>
$$d(\theta)=\left[\tfrac1n\sum_{j=1}^n\left(\ln X_j(\theta)-\bar y_j\right)^2\right]^{1/2}\le\varepsilon=0.15,$$ @@abc
<p>(rejection sampling, {{n_acc}} of {{Nd}} accepted). Adding the enzyme amount at the end of expression extends the sum by the term $\left(\ln E_6(\theta)-\overline{\ln E_{\mathrm{obs}}}\right)^2$ with $n+1$ in the denominator. The narrowing of a parameter is reported as $1-\mathrm{IQR}(\log_{10}\theta_{\mathrm{acc}})/\mathrm{IQR}(\log_{10}\theta_{\mathrm{prior}})$.</p>
<div class="prop"><strong>Proposition 4 (structural non identifiability).</strong> Let two designs share $K_X$ and $h_E$ and satisfy $a_Ek=a_E'k'$, and suppose their expression costs stay below $c_{\mathrm{flat}}=2H(0)-m_A-g_0/3={{cflat}}$. Then $R$, $X$ and $v_{\mathrm{DET}}$ coincide for all $t$, while $E'=(a_E'/a_E)\,E$.</div>
<div class="proof"><em>Proof.</em> Write $E=a_E\bar E$, where $d\bar E/dt=\chi\,h(R;K_R)-\delta\bar E$ does not involve $a_E$. Then $\kappa=(a_Ek)\bar EX/(K_m+X)$, so $R$, $X$ and $\kappa$ depend on $(a_E,k)$ only through the product. The expression cost does depend on $a_E$ alone, but in the expression stage $V(0,c)=2\min\{g_0,H,3(2H-m_A-c)\}-m_N=2g_0-m_N$ for every $c\le c_{\mathrm{flat}}$, so the cost does not change the delivered flux. $\square$</div>
<p>So, from toxin data alone, $\log a_E$ and $\log k$ are confounded along a line of slope $-1$. The accepted designs have an expression cost (5th to 95th percentile) of ${{burden_lo}}$ to ${{burden_hi}}$, below $c_{\mathrm{flat}}$, which is why the accepted set lies on that line (correlation of the logs ${{corr}}$). The only observation that separates them is one that depends on $a_E$ itself, the effector amount $E(T_{\mathrm{exp}})=a_E\bar E(T_{\mathrm{exp}})$, which is why the demo adds it. The product $a_Ek$, the sensor threshold and the half life are identifiable, and so is the protected fraction, which depends on the design only through them.</p>

<h2>10. Numerical verification</h2>
<p>Each formula above was recomputed from the demo's own results and parameters, and checked against the linear program or against the demo's output. The reference run is the true design and baseline G6PD decline.</p>
{{CHECKTABLE}}

<h2>11. Map from equations to code</h2>
<table>
<tr><th>Equation</th><th>In <code>coupling_demo.py</code></th></tr>
<tr><td>Network @@balance, bounds @@polytope</td><td><code>build_reactions</code>, <code>ToyGEM.bounds</code></td></tr>
<tr><td>Network capacity and delivered flux @@capacity, @@lemma</td><td><code>max_flux</code>, <code>metabolic_max</code>; the table <code>tabB_for</code>, <code>V2D</code></td></tr>
<tr><td>Marginal values (3.3)</td><td><code>marginal_values</code></td></tr>
<tr><td>Redox reserve @@reserve</td><td><code>oxidative_reserve</code>; <code>protected_fraction</code> uses the closed form</td></tr>
<tr><td>Window @@window and closing age @@closing</td><td>operating window block; <code>window_closing_age</code></td></tr>
<tr><td>Data models @@lognormal, @@regression</td><td><code>rep_lognormal</code>, <code>ci95</code>, <code>fit_log2_decay</code></td></tr>
<tr><td>Constraints and sampling @@tracer, @@hitrun</td><td><code>fva</code>, <code>sample_fluxes</code></td></tr>
<tr><td>Circuit @@ode, @@costkappa and coupling @@fluxmap, @@euler</td><td><code>run_coupled_lp</code> (linear program each step), <code>simulate_table</code> (table, many designs)</td></tr>
<tr><td>Ensemble and inference @@protected, @@abc</td><td>design ensemble block (<code>qmc.Sobol</code>), <code>protected_fraction</code></td></tr>
</table>
<p>The closed forms, the quasi steady state, the network capacity and the bound on design performance are derived here and are not in the code; the code computes the same quantities numerically, and section 10 shows that they agree.</p>

<h2>12. What is exact and what is approximate</h2>
<ul>
<li><strong>Exact for this toy model:</strong> Lemma 1; Proposition 1 and its marginal values; the reserve formula; the window and its closing ages; Proposition 4.</li>
<li><strong>Exact under the quasi steady state:</strong> Proposition 2 and the bound of Proposition 3. The check in section 10 measures the error (a few percent near the transition, where the toxin level moves faster than the enzyme and the network capacity).</li>
<li><strong>Numerical:</strong> forward Euler with $\Delta t=0.5$ h (changes the toxin level by less than 2% against $\Delta t=0.125$ h); hit and run sampling; rejection sampling for inference, which is not a full Bayesian posterior.</li>
<li><strong>Assumptions:</strong> the network is in steady state at every instant (so no cofactor ratios); the expression stage reuses the mature cell network; the data are generated by the same model that is fitted, so identifiability results are best cases.</li>
</ul>
"""

# ---------------------------------------------------------------- rendering
counter = [0]
EQ = {}


def render_display(m):
    tex, key = m.group(1).strip(), (m.group(2) or "").strip()
    if tex == "{{MATRIX}}":
        tex = matrix_tex()
        key = ""
    counter[0] += 1 if key else 0
    num = counter[0]
    if key:
        EQ[key] = num
    mml = convert(tex, display="block")
    tag = f'<span class="no" id="eq-{key}">({num})</span>' if key else ""
    return f'<div class="eq"><div class="m">{mml}</div>{tag}</div>'


def render_inline(m):
    return convert(m.group(1).strip(), display="inline")


body = BODY
NUM["x6_sim"] = f"{x6:.2f}"
for k, v in NUM.items():
    body = body.replace("{{" + k + "}}", v)
body = re.sub(r"\$\$(.+?)\$\$(?:[ \t]*@@([A-Za-z]+))?", render_display, body, flags=re.S)
body = re.sub(r"(?<!\\)\$(.+?)\$", render_inline, body)
body = re.sub(r"@@([A-Za-z]+)", lambda m: f'<a href="#eq-{m.group(1)}">({EQ[m.group(1)]})</a>', body)
rows = "".join(f'<tr><td>{html.escape(l)}</td><td>{html.escape(s)}</td><td class="{"ok" if ok else "bad"}">{"pass" if ok else "FAIL"}</td></tr>'
               for l, s, ok, _ in CHECKS)
body = body.replace("{{CHECKTABLE}}", f"<table><tr><th>Check</th><th>Result</th><th></th></tr>{rows}</table>")
assert "{{" not in body, re.findall(r"\{\{[^}]*\}\}", body)[:5]
assert all(ok for *_, ok, _ in CHECKS), [c[0] for c in CHECKS if not c[2]]

CSS = """
:root{--bg:#fbfbfa;--panel:#f1f4f5;--ink:#1e2430;--muted:#5b6472;--line:#d5dbe1;--teal:#1d6a73;--teal-t:#e3eff0;--red:#9e1b32;--red-t:#f7e8eb;--warn:#fff6e0;--warn-line:#e3b04b;--prop:#eef5f6}
@media (prefers-color-scheme:dark){:root{--bg:#14181f;--panel:#1c222b;--ink:#e7eaee;--muted:#9aa4b2;--line:#323b47;--teal:#6fc0c9;--teal-t:#1a2f33;--red:#f08a9b;--red-t:#35191f;--warn:#33290f;--warn-line:#a37a1e;--prop:#18262a}}
*{box-sizing:border-box}
body{margin:0;background:var(--bg);color:var(--ink);font:16px/1.65 -apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,Helvetica,Arial,sans-serif}
main{max-width:980px;margin:0 auto;padding:28px 20px 80px}
h1{font:700 2rem/1.2 Georgia,"Times New Roman",serif;margin:0 0 6px}
h2{font:700 1.4rem/1.25 Georgia,"Times New Roman",serif;margin:2.3em 0 .5em;padding-top:.4em;border-top:1px solid var(--line)}
h3{font-size:1.02rem;margin:1.6em 0 .3em;color:var(--teal)}
p{margin:.6em 0}
.sub{color:var(--muted);margin:0 0 18px}
.banner{background:var(--warn);border:1px solid var(--warn-line);border-radius:10px;padding:12px 16px;margin:16px 0}
.box{background:var(--teal-t);border-radius:12px;padding:16px 20px;margin:18px 0}
.box h3{margin-top:0}
ul{margin:.4em 0 .6em 1.2em;padding:0}
li{margin:.3em 0}
table{border-collapse:collapse;width:100%;margin:12px 0;font-size:.93rem}
th,td{border-bottom:1px solid var(--line);padding:7px 9px;text-align:left;vertical-align:top}
th{background:var(--panel);font-weight:600}
code{font:.88em ui-monospace,SFMono-Regular,Menlo,Consolas,monospace;background:var(--panel);padding:1px 5px;border-radius:4px}
math{font-size:1.1em}
.eq{display:grid;grid-template-columns:1fr auto;align-items:center;gap:12px;margin:.9em 0;overflow-x:auto}
.eq .m{overflow-x:auto;padding:2px 0}
.eq .no{color:var(--muted);font-size:.95rem;white-space:nowrap}
.prop{background:var(--prop);border-left:4px solid var(--teal);border-radius:6px;padding:10px 16px;margin:16px 0}
.proof{color:var(--ink);border-left:2px solid var(--line);padding:2px 0 2px 14px;margin:8px 0 16px;font-size:.96rem}
details{margin:.6em 0}
summary{cursor:pointer;font-weight:600}
a{color:var(--teal);text-decoration:none}
.ok{color:var(--teal);font-weight:700}
.bad{color:var(--red);font-weight:700}
"""

page = (f'<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">'
        f"<title>Mathematical formulation of the coupling demonstration</title><style>{CSS}</style></head>"
        f"<body><main>{body}</main></body></html>")
out = ROOT / "docs" / "mathematics.html"
out.parent.mkdir(exist_ok=True)
out.write_text(page, encoding="utf-8")
print("wrote", out, f"{len(page) / 1e3:.0f} kB;", len(EQ), "numbered equations;", len(CHECKS), "checks, all pass")
