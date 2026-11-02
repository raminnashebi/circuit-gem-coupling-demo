#!/usr/bin/env python3
"""Build docs/report.html from results/results.json and results/figures/ (run src/coupling_demo.py first)."""
import base64
import html
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
R = json.loads((ROOT / "results" / "results.json").read_text())


def f2(x, d=2):
    return f"{x:.{d}f}"


def img(name, alt, caption):
    b64 = base64.b64encode((ROOT / "results" / "figures" / name).read_bytes()).decode()
    return (f'<figure><img alt="{html.escape(alt)}" src="data:image/png;base64,{b64}">'
            f"<figcaption>{caption}</figcaption></figure>")


L1 = R["level1_capacity"]
W = R["level1_window"]
S = R["sampling"]
PB = R["prediction_bands"]
T = R["level2_truth_design"]
EN = R["ensemble"]
ID = R["identifiability"]
CK = R["checks"]
DF = R["data_fits"]
CP = R["circuit_parameters"]
BP = R["base_parameters"]
cap_n = {a: L1["NADPH"][a]["network_capacity"] for a in ("0", "60", "120")}
cap_a = {a: L1["ATP"][a]["network_capacity"] for a in ("0", "60", "120")}
mv_n, mv_a = L1["NADPH"]["0"]["marginal"], L1["ATP"]["0"]["marginal"]
scen = list(EN["best_protected_fraction_per_scenario"].keys())
cm = ID["contraction_toxin_data"]
cm2 = ID["contraction_plus_enzyme_amount"]

rx_rows = "".join(
    f"<tr><td><code>{html.escape(r['name'])}</code></td><td>{html.escape(' + '.join(f'{abs(c):.3g} {m}' if abs(c) != 1 else m for m, c in r['stoichiometry'].items() if c < 0)) or '(source)'}"
    f" &rarr; {html.escape(' + '.join(f'{abs(c):.3g} {m}' if abs(c) != 1 else m for m, c in r['stoichiometry'].items() if c > 0)) or '(sink)'}</td>"
    f"<td>{html.escape(r['gene_rule'])}</td></tr>" for r in R["network"]["reactions"])

CSS = """
:root{--bg:#fbfbfa;--panel:#f1f4f5;--ink:#1e2430;--muted:#5b6472;--line:#d5dbe1;--teal:#1d6a73;--teal-t:#e3eff0;--red:#9e1b32;--red-t:#f7e8eb;--warn:#fff6e0;--warn-line:#e3b04b}
@media (prefers-color-scheme:dark){:root:not([data-theme="light"]){--bg:#14181f;--panel:#1c222b;--ink:#e7eaee;--muted:#9aa4b2;--line:#323b47;--teal:#6fc0c9;--teal-t:#1a2f33;--red:#f08a9b;--red-t:#35191f;--warn:#33290f;--warn-line:#a37a1e}}
:root[data-theme="dark"]{--bg:#14181f;--panel:#1c222b;--ink:#e7eaee;--muted:#9aa4b2;--line:#323b47;--teal:#6fc0c9;--teal-t:#1a2f33;--red:#f08a9b;--red-t:#35191f;--warn:#33290f;--warn-line:#a37a1e}
*{box-sizing:border-box}
body{margin:0;background:var(--bg);color:var(--ink);font:16px/1.6 -apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,Helvetica,Arial,sans-serif}
main{max-width:980px;margin:0 auto;padding:28px 20px 80px}
h1{font:700 2rem/1.2 Georgia,"Times New Roman",serif;margin:0 0 6px}
h2{font:700 1.45rem/1.25 Georgia,"Times New Roman",serif;margin:2.4em 0 .5em;padding-top:.4em;border-top:1px solid var(--line)}
h3{font-size:1.05rem;margin:1.6em 0 .3em;color:var(--teal)}
p{margin:.6em 0}
.sub{color:var(--muted);margin:0 0 18px}
.banner{background:var(--warn);border:1px solid var(--warn-line);border-radius:10px;padding:12px 16px;margin:16px 0}
.box{background:var(--teal-t);border-radius:12px;padding:16px 20px;margin:18px 0}
.box h3{margin-top:0}
.red{background:var(--red-t)}
.red h3{color:var(--red)}
ul,ol{margin:.4em 0 .6em 1.2em;padding:0}
li{margin:.25em 0}
table{border-collapse:collapse;width:100%;margin:12px 0;font-size:.93rem}
th,td{border-bottom:1px solid var(--line);padding:7px 9px;text-align:left;vertical-align:top}
th{background:var(--panel);font-weight:600}
code{font:.88em ui-monospace,SFMono-Regular,Menlo,Consolas,monospace;background:var(--panel);padding:1px 5px;border-radius:4px}
pre{background:var(--panel);padding:12px 14px;border-radius:8px;overflow-x:auto;font:.88rem/1.5 ui-monospace,Menlo,Consolas,monospace}
figure{margin:18px 0}
img{max-width:100%;height:auto;border:1px solid var(--line);border-radius:8px;background:#fff}
figcaption{font-size:.9rem;color:var(--muted);margin-top:6px}
.num{font-weight:700;color:var(--teal)}
.bad{font-weight:700;color:var(--red)}
.cols{display:grid;grid-template-columns:1fr 1fr;gap:14px}
@media(max-width:700px){.cols{grid-template-columns:1fr}}
.tag{display:inline-block;font-size:.78rem;font-weight:700;letter-spacing:.05em;text-transform:uppercase;color:var(--teal);margin-right:.5em}
details{margin:.6em 0}
summary{cursor:pointer;font-weight:600}
svg text{font-family:-apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,Helvetica,Arial,sans-serif}
"""

DIAGRAM = """
<svg viewBox="-40 0 940 400" role="img" aria-label="Diagram: the circuit model and the metabolic network exchange the enzyme bound, the delivered flux and the toxin level" style="width:100%;height:auto">
<defs><marker id="ar" viewBox="0 0 10 10" refX="8" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse"><path d="M0 0L10 5L0 10z" fill="currentColor"/></marker></defs>
<g fill="none" stroke-width="2">
<rect x="20" y="20" width="560" height="140" rx="14" style="fill:var(--red-t);stroke:var(--red)"/>
<rect x="20" y="230" width="560" height="140" rx="14" style="fill:var(--teal-t);stroke:var(--teal)"/>
</g>
<text x="40" y="48" font-size="16" font-weight="700" style="fill:var(--red)">Circuit model: ordinary differential equations</text>
<g font-size="14" text-anchor="middle" style="fill:var(--ink)">
<rect x="45" y="75" width="130" height="50" rx="8" style="fill:var(--bg);stroke:var(--red)" stroke-width="1.5"/><text x="110" y="105">Sensor (toxin level)</text>
<rect x="235" y="75" width="130" height="50" rx="8" style="fill:var(--bg);stroke:var(--red)" stroke-width="1.5"/><text x="300" y="105">Regulator R</text>
<rect x="425" y="75" width="130" height="50" rx="8" style="fill:var(--bg);stroke:var(--red)" stroke-width="1.5"/><text x="490" y="105">Effector enzyme E</text>
</g>
<g stroke-width="2" style="stroke:var(--red);color:var(--red)"><line x1="177" y1="100" x2="231" y2="100" marker-end="url(#ar)"/><line x1="367" y1="100" x2="421" y2="100" marker-end="url(#ar)"/></g>
<text x="40" y="258" font-size="16" font-weight="700" style="fill:var(--teal)">Metabolic network: linear programming at each step</text>
<g font-size="14" style="fill:var(--ink)">
<text x="45" y="296">Steady state S v = 0, bounds on every flux</text>
<text x="45" y="320">Maintenance tasks: pump ATP, basal NADPH (no growth)</text>
<text x="45" y="344">Capacities fall with cell age (G6PD, hexokinase)</text>
</g>
<g stroke-width="2.5" fill="none">
<line x1="470" y1="130" x2="470" y2="224" style="stroke:var(--red);color:var(--red)" marker-end="url(#ar)"/>
<line x1="530" y1="224" x2="530" y2="130" style="stroke:var(--teal);color:var(--teal)" marker-end="url(#ar)"/>
<path d="M20 330 H-10 V100 H41" style="stroke:var(--muted);color:var(--muted)" stroke-dasharray="6 5" marker-end="url(#ar)"/>
</g>
<g font-size="13.5" style="fill:var(--ink)">
<text x="610" y="165" font-weight="700" style="fill:var(--red)">Circuit to network</text>
<text x="610" y="185">enzyme level sets the flux bound</text>
<text x="610" y="204">v(DETOX) &le; kcat &times; E(t)</text>
<text x="610" y="232" font-weight="700" style="fill:var(--teal)">Network to circuit</text>
<text x="610" y="252">flux the cell can really deliver,</text>
<text x="610" y="271">spare ATP and NADPH</text>
<text x="610" y="299" font-weight="700" style="fill:var(--muted)">Feedback</text>
<text x="610" y="319">delivered flux changes the toxin</text>
<text x="610" y="338">level, which the sensor reads</text>
</g>
<text transform="translate(-20,215) rotate(-90)" text-anchor="middle" font-size="12.5" style="fill:var(--muted)">toxin level X</text>
<text x="450" y="394" font-size="12.5" text-anchor="middle" style="fill:var(--muted)">Two stages: expression in the erythroblast (circuit on, costs ATP), then the mature red cell (no new protein, enzyme slowly decays)</text>
</svg>
"""

BODY = f"""
<h1>Coupling a gene circuit to a metabolic model</h1>
<p class="sub">A toy demonstration of how a gene circuit model can be coupled to a metabolic network model, with invented numbers and simulated data.</p>
<div class="banner"><strong>Read this first.</strong> The cell, the circuit and every "measurement" here are made up or simulated by the script that wrote this page. The demonstration shows how the method works and what kind of answers it gives. It says nothing about real red cells, and no number here should be quoted as a biological result.</div>

<div class="box"><h3>In one minute</h3>
<ul>
<li><strong>The idea.</strong> Treat the gene circuit and the cell as two models that swap a few variables. The circuit tells the cell how much enzyme there is. The cell tells the circuit how much flux it can really deliver.</li>
<li><strong>What the demonstration does.</strong> A small red cell style network, an enzyme circuit that removes a toxin using NADPH, five simulated datasets, and six steps from a one-off feasibility check to a search over {EN['n_designs']:,} designs.</li>
<li><strong>What came out.</strong> The network's capacity sets an upper bound on what any design can achieve. The toxin stays controlled for about <span class="num">{T['toxin_loses_control_at_age_d']:.0f} days</span>, but the cell's spare redox capacity falls below the safety threshold at about <span class="num">{T['reserve_below_threshold_at_age_d']:.0f} days</span>. Expression strength and kcat cannot be told apart from the toxin data, yet the prediction that matters (protected life) barely changes across the designs the data allow.</li>
</ul></div>

<h2>1. The idea</h2>
{DIAGRAM}
<p>Two things differ from a standard metabolic model. First, a mature red cell does not grow, so there is no biomass objective. The network is asked to meet maintenance tasks (pump ATP and a basal NADPH demand) and then to say how much extra the circuit can draw. Second, the circuit is made in one cell stage and acts in another, so the same coupling is run in two regimes.</p>

<h2>2. The method in six steps</h2>
<table>
<tr><th>#</th><th>Step</th><th>What happens</th><th>Question it answers</th></tr>
<tr><td>1</td><td>Add the effector to the network</td><td>A new reaction <code>DETOX</code> is tied to the transgene by a gene rule. Its upper bound is the enzyme capacity, <code>v &le; kcat &times; E</code>.</td><td>How does the circuit enter the model?</td></tr>
<tr><td>2</td><td>Ask the network three questions</td><td>Linear programming gives the <em>network capacity</em> (largest feasible DETOX flux, the maximum in flux variability analysis), <em>marginal values</em> (what one more unit of each resource would buy, a shadow price) and the <em>redox reserve</em> left after the effector runs.</td><td>Is the effector feasible? What limits it? Is there a safe window?</td></tr>
<tr><td>3</td><td>Constrain with data</td><td>Exchange rate intervals, a tracer fraction and enzyme activities bound the network. Flux sampling replaces a single objective.</td><td>What do the data allow? Which measurement tightens a prediction?</td></tr>
<tr><td>4</td><td>Couple in time</td><td>Each hour the circuit ODEs propose a flux bound, the linear program returns what the cell can deliver, and the delivered flux feeds back into the toxin balance.</td><td>When does the design fail, and why?</td></tr>
<tr><td>5</td><td>Search designs</td><td>{EN['n_designs']:,} Sobol designs are simulated with a table of exact linear program results.</td><td>Which designs hold up across uncertainty in the host?</td></tr>
<tr><td>6</td><td>Identifiability</td><td>Designs that reproduce a simulated toxin time course are kept (rejection sampling). Parameter ranges show what the data fix.</td><td>Which new measurement should the lab make?</td></tr>
</table>
<p>Step 5 relies on a short argument. The set of feasible DETOX fluxes is an interval from zero to the network capacity, because the feasible region is convex and zero is feasible. So the flux the cell delivers is exactly <code>min(enzyme bound, network capacity)</code>, and the network capacity depends only on cell age and expression cost. That lets the design search use a table of exact results instead of a linear program per step. The checks in section 7 confirm that the table and the direct linear program agree.</p>

<h2>3. The toy cell</h2>
<p>{R['checks']['n_metabolites']} metabolites and {R['checks']['n_reactions']} reactions, in arbitrary flux units. A young cell takes up about 0.4 glucose units. Steady state means <code>S v = 0</code>. ATP/ADP and NADPH/NADP are conserved pairs. The pentose phosphate route (PPP) is lumped with its recycling and carbon dioxide is not tracked. Gene symbols are illustrative.</p>
<table><tr><th>Reaction</th><th>Stoichiometry</th><th>Gene rule</th></tr>{rx_rows}</table>
<p><strong>Maintenance tasks instead of growth.</strong> The pump reaction <code>ATPASE</code> must carry at least {BP['m_atp']} and the redox reaction <code>GSH_REDOX</code> at least {BP['m_nadph']}. <strong>Ageing.</strong> G6PD capacity halves every {BP['half_lives_days']['G6PD']:.0f} days and hexokinase every {BP['half_lives_days']['HK']:.0f} days. <strong>Network facts that follow.</strong> NADPH comes only from the PPP, so the network capacity for DETOX is exactly <code>2 &times; G6PD capacity &minus; basal NADPH demand</code>. At age 0 that is 2 &times; {BP['ppp_cap']} &minus; {BP['m_nadph']} = {f2(cap_n['0'])}, which the linear program reproduces.</p>

<h2>4. The circuit model</h2>
<p>A toxin <em>X</em> arrives at a constant load and is also cleared slowly by other routes. The sensor reads <em>X</em>, a regulator <em>R</em> follows, and the effector enzyme <em>E</em> removes <em>X</em> using NADPH.</p>
<pre>dR/dt = on(t) &middot; a_R &middot; h(X; K_X) &minus; d_R &middot; R                      h(u; K) = u&sup2; / (K&sup2; + u&sup2;)
dE/dt = on(t) &middot; a_E &middot; h(R; K_R) &minus; d_E(t) &middot; E
dX/dt = L &minus; v_DETOX &minus; k_c &middot; X
bound:  v_DETOX &le; kcat &middot; E &middot; X / (Km + X)     then the network returns what it can deliver
cost:   ATP drawn by expression = c_R &middot; a_R &middot; h(X; K_X) + c_E &middot; a_E &middot; h(R; K_R)   (expression stage only)</pre>
<ul>
<li><strong>Expression stage (6 days, stands in for the erythroblast).</strong> <code>on = 1</code>. Making proteins costs ATP, which the network must supply. The same toy network stands in for the erythroblast. In the real project this stage needs its own model, with growth and mitochondria.</li>
<li><strong>Mature stage (120 days).</strong> <code>on = 0</code>. There are no ribosomes, so no new protein: <em>R</em> disappears in about two days and <em>E</em> decays with a half-life of {R['circuit_parameters']['true_design']['hl']:.0f} days.</li>
<li><strong>Targets (invented).</strong> Toxin level at most {CP['X_target']:.0f} (it would be {CP['L']/CP['kc']:.0f} without the circuit). That needs a DETOX flux of at least <span class="num">{f2(W['need_flux'])}</span>. Protection also needs a redox reserve of at least <span class="num">{W['reserve_min_fold']}&times;</span> the basal NADPH flux.</li>
<li><strong>True design used to make the data:</strong> expression strength a<sub>E</sub> = {CP['true_design']['a_E']}, kcat = {CP['true_design']['kcat']}, sensor threshold K<sub>X</sub> = {CP['true_design']['K_X']}, enzyme half-life {CP['true_design']['hl']:.0f} d.</li>
</ul>

<h2>5. The synthetic data</h2>
<table>
<tr><th>Data</th><th>What was "measured"</th><th>How it was made</th><th>Used for</th></tr>
<tr><td>D1 exchange rates</td><td>Glucose uptake and lactate release, 7 cell ages</td><td>Parsimonious solution of the true network, 4 replicates, 6% noise</td><td>Bounds on the network (95% intervals)</td></tr>
<tr><td>D2 tracer fraction</td><td>Share of glucose through the PPP (as from a [1,2-<sup>13</sup>C<sub>2</sub>]glucose tracer)</td><td>True fraction, 4 replicates, 15% noise</td><td>Ratio constraint on the PPP</td></tr>
<tr><td>D3 enzyme activity</td><td>G6PD and hexokinase activity against cell age</td><td>True decay, 3 replicates, 10% noise</td><td>Capacities and their half-lives</td></tr>
<tr><td>D4 toxin time course</td><td>Toxin level on 11 days (2 in the expression stage, 9 in the mature cell)</td><td>True design, 3 replicates, 12% noise</td><td>Accepting or rejecting designs</td></tr>
<tr><td>D5 enzyme amount</td><td>Effector level at the end of the expression stage</td><td>True value, 3 replicates, 15% noise</td><td>The extra measurement tested in step 6</td></tr>
</table>

<h2>6. Results</h2>
<h3>6.1 Step 2: what limits the effector</h3>
{img('fig1_capacity.png', 'Delivered DETOX flux against enzyme capacity for three cell ages and two effector types', 'Delivered flux rises with the enzyme capacity until the network cannot supply more, then plateaus. Left: an effector that uses NADPH. Right: one that uses ATP. Each plateau is the network capacity at that cell age.')}
<p><strong>NADPH effector.</strong> The network capacity falls from <span class="num">{f2(cap_n['0'])}</span> at age 0 to <span class="num">{f2(cap_n['60'])}</span> at 60 days and <span class="num">{f2(cap_n['120'])}</span> at 120 days, as G6PD declines. The marginal values name the limit: one more unit of G6PD capacity buys <span class="num">{f2(mv_n['G6PD capacity (PPP)'],1)}</span> units of detox flux, one unit less basal NADPH demand buys <span class="num">{f2(mv_n['NADPH supply (less basal demand)'],1)}</span>, and hexokinase, ATP and glucose are worth zero.</p>
<p><strong>ATP effector.</strong> A different host constraint binds. The network capacity is {f2(cap_a['0'])} at age 0 and {f2(cap_a['120'])} at 120 days. Here hexokinase capacity is worth <span class="num">{f2(mv_a['hexokinase capacity'],1)}</span> per unit and G6PD nothing. The same pipeline tells you which cofactor limits which kind of effector.</p>
<p class="sub">Reading note: the dashed line is the flux needed to hold the toxin at its target. At 120 days even a perfectly expressed enzyme can no longer reach it.</p>

<h3>6.2 Step 2: the safe operating window</h3>
{img('fig2_window.png', 'Operating window of the effector against cell age and enzyme capacity', 'For each enzyme level and cell age, the network gives the delivered flux and the redox reserve left. Grey: too little detox. Red: detox is enough but the spare redox capacity is below the threshold. Blue: both are met.')}
<p>The window is bounded from below by the detox needed and from above by the reserve that must be kept. It narrows as the cell ages and closes at about <span class="num">{W['closing_age_days']:.0f} days</span>. After about 107 days no enzyme level reaches the needed flux at all. This is the kind of figure that tells an experimental team how strongly to express an enzyme, and in which cells it will work.</p>

<h3>6.3 Step 3: what each dataset tells us</h3>
{img('fig3_data.png', 'Synthetic exchange data, flux sampling with and without the tracer, and predicted network capacity for DETOX with and without the enzyme assay', 'Left: simulated exchange rates with 95% intervals (lines are the true values; they are flat because the demand of the toy cell does not change). Middle: objective free flux sampling at cell age 60 days. Right: the predicted network capacity for DETOX, from data that fix fluxes only (grey) and after adding the enzyme activity assay (teal).')}
<ul>
<li><strong>Exchange rates alone cannot fix the PPP split.</strong> Sampling the flux space gives a PPP share of glucose between {f2(S['exchange_only'][0])} and {f2(S['exchange_only'][2])} (95% range). Adding the tracer fraction narrows it to {f2(S['with_tracer'][0])} to {f2(S['with_tracer'][2])}. The true value, {f2(S['true_ppp_fraction'],3)}, is inside. The feasible range of the PPP flux itself shrinks from {f2(S['fva_ppp_flux_exchange_only'][0],2)}&ndash;{f2(S['fva_ppp_flux_exchange_only'][1],2)} to {f2(S['fva_ppp_flux_with_tracer'][0],2)}&ndash;{f2(S['fva_ppp_flux_with_tracer'][1],2)}.</li>
<li><strong>Fluxes tell you what the cell does now, enzyme levels tell you what it can do under load.</strong> The network capacity for the effector depends on G6PD capacity, which exchange rates and a basal tracer fraction do not measure. At age 60 the predicted network capacity is {f2(PB['age_60_without_enzyme_assay_95'][0])} to {f2(PB['age_60_without_enzyme_assay_95'][1])} without the enzyme assay (my invented prior) and <span class="num">{f2(PB['age_60_with_enzyme_assay_95'][0],3)} to {f2(PB['age_60_with_enzyme_assay_95'][1],3)}</span> with it. The true value is {f2(PB['age_60_truth'],3)}. The fitted G6PD half-life is {f2(DF['G6PD_half_life_d'],1)} days (95% interval {f2(DF['G6PD_half_life_95'][0],0)}&ndash;{f2(DF['G6PD_half_life_95'][1],0)}) against a true {DF['G6PD_half_life_true']:.0f}.</li>
</ul>
<p>This is the practical reason for the "what to measure" list: in a real project, it is the enzyme activity or abundance by cell age, not another round of exchange rates, that narrows the prediction.</p>

<h3>6.4 Step 4: coupling in time</h3>
{img('fig4_dynamics.png', 'Time course of regulator, DETOX flux, toxin level and redox reserve in the expression stage and the mature cell', 'The coupled run for the true design. Left column: the six day expression stage. Right column: 120 days of the mature cell. From top: regulator, DETOX flux (enzyme bound, network capacity, delivered flux), toxin level, spare redox capacity.')}
<ol>
<li><strong>Expression stage.</strong> The sensor reads the toxin and switches the circuit on. The toxin falls from {CP['L']/CP['kc']:.0f} to about {f2(T['x_at_end_of_expression'],1)}. The ATP drawn for making proteins peaks at {f2(T['max_expression_cost'])} flux units and the cell never has less than {f2(T['min_atp_spare_during_expression'])} spare. No crisis occurs.</li>
<li><strong>Maturation.</strong> With no ribosomes, the regulator falls to zero within about two days and the circuit can no longer respond. Only the enzyme pool is left and it decays slowly. This is the two regime picture from the slides.</li>
<li><strong>The cost to the host.</strong> The detox flux uses NADPH, so the redox reserve falls as G6PD declines. It crosses the <span class="num">{W['reserve_min_fold']}&times;</span> threshold at about <span class="bad">{T['reserve_below_threshold_at_age_d']:.0f} days</span>. The toxin stays under its target until about <span class="num">{T['toxin_loses_control_at_age_d']:.0f} days</span>, when the network itself becomes the limit (from about {T['metabolism_becomes_limiting_at_age_d']:.0f} days).</li>
</ol>
<p>The lesson is the point of coupling: judged by the circuit output alone, the design looks successful for more than three months. Only the coupled model shows that the safety margin of the host is gone after about seven weeks.</p>

<h3>6.5 Steps 5 and 6: designs and identifiability</h3>
{img('fig5_designs.png', 'Plausible designs on the expression and kcat plane, narrowing of parameter ranges by the data, and design score', 'Left: all designs (grey) and those that reproduce the simulated toxin time course (teal) lie along a line; adding the enzyme amount keeps a short segment (orange). Middle: how much each range narrows. Right: protected fraction of life, worst of three G6PD decline scenarios.')}
<ul>
<li><strong>The network caps every design.</strong> Of {EN['n_designs']:,} designs, {EN['designs_failing_from_expression_cost']*100:.0f}% fail because expression costs more ATP than the cell has to spare. For the rest, the best protected fraction of mature life is {f2(EN['best_protected_fraction_per_scenario'][scen[0]])}, {f2(EN['best_protected_fraction_per_scenario'][scen[1]])} and {f2(EN['best_protected_fraction_per_scenario'][scen[2]])} in the fast, baseline and slow G6PD scenarios. These equal the network's own limits ({f2(EN['network_bound_protected_fraction'][scen[0]])}, {f2(EN['network_bound_protected_fraction'][scen[1]])}, {f2(EN['network_bound_protected_fraction'][scen[2]])}) to within a day. Tuning the circuit cannot beat the host.</li>
<li><strong>Expression and kcat are not identifiable separately.</strong> {ID['n_accepted']} designs reproduce the toxin time course, and they lie on a line (correlation of the logs {f2(ID['corr_log_aE_log_kcat_accepted'],2)}). The data narrow the <em>product</em> a<sub>E</sub> &times; kcat by {cm['a_E × kcat']*100:.0f}% but expression by only {cm['a_E']*100:.0f}% and kcat by {cm['kcat']*100:.0f}%. The reason is structural: in the mature cell and in the toxin balance, enzyme level only acts through kcat &times; E.</li>
<li><strong>One more measurement resolves it.</strong> Adding the simulated enzyme amount at the end of expression narrows expression to {cm2['a_E']*100:.0f}% and kcat to {cm2['kcat']*100:.0f}% (from {cm['a_E']*100:.0f}% and {cm['kcat']*100:.0f}%), and the ATP burden during expression from {f2(ID['expression_burden_accepted_designs_5_50_95'][0])}&ndash;{f2(ID['expression_burden_accepted_designs_5_50_95'][2])} to {f2(ID['expression_burden_accepted_with_enzyme_amount_5_50_95'][0])}&ndash;{f2(ID['expression_burden_accepted_with_enzyme_amount_5_50_95'][2])} (true value {f2(ID['expression_burden_truth_design'])}). This burden is where expression and kcat differ: the cheaper design is the one with higher kcat and lower expression.</li>
<li><strong>The prediction is identifiable even though the parameters are not.</strong> Across all designs consistent with the toxin data, the protected fraction of life lies between {f2(ID['protected_fraction_of_accepted_designs_baseline'][0])} and {f2(ID['protected_fraction_of_accepted_designs_baseline'][2])} (true design {f2(ID['protected_fraction_truth_design'])}). This is the same separation of robust predictions from non-identifiable parameters as in the bioenergetics ensembles.</li>
</ul>

<h2>7. Checks</h2>
<table>
<tr><th>Check</th><th>Result</th></tr>
<tr><td>ATP/ADP and NADPH/NADP pools conserved; steady state residual of a solution</td><td>{'yes' if CK['adenylate_pool_conserved'] and CK['nadp_pool_conserved'] else 'NO'}; residual {CK['baseline_steady_state_residual']:.0e}</td></tr>
<tr><td>Table of exact linear program results against direct linear programs (network capacity by age and by expression cost)</td><td>largest difference {CK['table_vs_exact_lp_max_abs_diff']:.1e}; {CK['cost_table_vs_exact_lp_max_abs_diff']:.1e}</td></tr>
<tr><td>Coupled run with a linear program at every hour against the same run from the table (toxin level)</td><td>largest relative difference {CK['lp_run_vs_table_run_max_rel_diff_X']:.0e}</td></tr>
<tr><td>Redox reserve from the closed form <code>1 + (capacity &minus; v) / basal</code> against the linear program</td><td>largest difference {CK['reserve_closed_form_vs_lp_max_abs_diff']:.1e}</td></tr>
<tr><td>Time step 0.5 h against 0.125 h (toxin level)</td><td>largest relative difference {CK['dt_0.5_vs_0.125_max_rel_diff_X']*100:.1f}%</td></tr>
<tr><td>Flux sampler: constraint violation; difference between two independent chains (in standard deviations)</td><td>{S['max_steady_state_violation']:.0e}; {S['between_chain_mean_diff_over_sd']:.2f}</td></tr>
<tr><td>True network capacity inside the 95% band that uses the enzyme assay (ages 0 to 120 d)</td><td>{PB['truth_inside_post_band_fraction']*100:.0f}% of ages</td></tr>
<tr><td>Fitted G6PD half-life against the truth used to make the data</td><td>{f2(DF['G6PD_half_life_d'],1)} d (95% {f2(DF['G6PD_half_life_95'][0],0)}&ndash;{f2(DF['G6PD_half_life_95'][1],0)}) against {DF['G6PD_half_life_true']:.0f} d</td></tr>
</table>

<h2>8. Limits of this demonstration</h2>
<ul>
<li><strong>Invented cell.</strong> The network is a teaching toy. It leaves out NADH balance, 2,3-BPG, ion and volume dynamics, carbon dioxide and compartments. Its numbers (capacities, demands, thresholds) were chosen to be sensible, not measured.</li>
<li><strong>Same model makes and fits the data.</strong> The identifiability results are therefore best cases. Real data add model error and noise that is not lognormal.</li>
<li><strong>Quasi steady state.</strong> A flux balance model assumes metabolism settles much faster than the hour to day changes in the circuit. That suits fluxes, but it gives no cofactor ratios (NADPH/NADP&#8314;, GSH/GSSG). A real application would need a small kinetic core for those. The redox reserve here is a flux margin standing in for them.</li>
<li><strong>Simple circuit.</strong> One toxin, one effector, no cell to cell variation, no protein burden beyond ATP, no loss of function from damage.</li>
<li><strong>Simple inference.</strong> The design search accepts designs whose predictions fall within a tolerance of the data ({ID['n_accepted']} of {EN['n_designs']:,}). This shows structure, not a full Bayesian posterior.</li>
<li><strong>One network for two stages.</strong> The expression stage reuses the mature cell network.</li>
</ul>

<h2>9. Questions and next steps</h2>
<h3>Questions it may prompt</h3>
<ul>
<li><em>Is steady state valid at each step?</em> For flux questions yes, since central metabolism relaxes in seconds or minutes. For cofactor ratios I would add a kinetic core of central red cell metabolism.</li>
<li><em>Why a linear program per step and not a kinetic model?</em> The network is large and its kinetics are partly unknown. Constraints need no rate laws. Where dynamics matter, a kinetic model is used for the core only.</li>
<li><em>How would this scale?</em> Replace the toy matrix with RBC-GEM for the mature cell, with stage extensions from Human-GEM such as mitochondrial reactions for the reticulocyte. The functions need only a stoichiometric matrix, bounds and gene rules, so the same code runs on a model loaded with COBRApy.</li>
<li><em>What if the real circuit is not an enzyme?</em> Any circuit output that changes a flux bound, a demand reaction or a regulatory switch enters the same way.</li>
</ul>
<h3>Next steps for a real application</h3>
<ol>
<li>Swap in a real stage specific red cell model and literature capacities (enzyme activities by cell age).</li>
<li>Add a kinetic core for central metabolism to get NADPH/NADP&#8314; and ATP ratios.</li>
<li>Use measured metabolomics as exchange intervals, as in step 3, and test identifiability before the experiments.</li>
<li>Give the erythroblast stage its own model with growth and mitochondria.</li>
</ol>

<h2>10. Files and how to reproduce</h2>
<ul>
<li><code>src/coupling_demo.py</code>: the whole demonstration; writes <code>results/results.json</code> (every number quoted above) and <code>results/figures/</code>.</li>
<li><code>src/build_report.py</code>: builds this page, <code>docs/report.html</code>.</li>
<li><code>src/build_math.py</code>: builds <code>docs/mathematics.html</code>, the equations, derivations and numerical checks behind the demonstration.</li>
<li><code>tests/</code>: checks the internal consistency of the results and every headline number.</li>
</ul>
<pre>pip install -r requirements.txt
make all    # results, figures and both reports
make test   # the checks</pre>
<details><summary>Glossary</summary>
<p><strong>Metabolic network (GEM).</strong> The list of reactions a cell can run, as a stoichiometric matrix <code>S</code>. <strong>Flux.</strong> Rate of a reaction. <strong>Linear programming.</strong> Optimising a linear goal under linear constraints. <strong>Network capacity.</strong> Largest feasible flux of one reaction (the maximum from FVA). <strong>Marginal value (shadow price).</strong> Gain in the goal from relaxing one constraint by one unit. <strong>FVA.</strong> Flux variability analysis, the minimum and maximum of each flux. <strong>Flux sampling.</strong> Drawing random feasible flux states without choosing an objective. <strong>G6PD, PPP.</strong> The enzyme and pathway that make NADPH from glucose. <strong>Effector.</strong> The enzyme the circuit produces. <strong>Identifiable.</strong> A parameter is identifiable if the data can fix it.</p>
</details>
"""

page = (f'<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">'
        f"<title>Gene circuit and metabolic model coupling</title><style>{CSS}</style></head><body><main>{BODY}</main></body></html>")
out = ROOT / "docs" / "report.html"
out.parent.mkdir(exist_ok=True)
out.write_text(page, encoding="utf-8")
print("wrote", out, f"{len(page) / 1e6:.2f} MB")
