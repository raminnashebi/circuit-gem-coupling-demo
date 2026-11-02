"""Checks on results/results.json: internal consistency and the headline numbers in the README and reports."""
import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
R = json.loads((ROOT / "results" / "results.json").read_text())
C = R["checks"]


def test_network_is_balanced():
    assert C["adenylate_pool_conserved"] and C["nadp_pool_conserved"]
    assert C["baseline_steady_state_residual"] < 1e-9
    assert (C["n_reactions"], C["n_metabolites"], C["rank_S"]) == (12, 10, 8)


def test_lookup_tables_match_exact_linear_programs():
    assert C["table_vs_exact_lp_max_abs_diff"] < 1e-6
    assert C["cost_table_vs_exact_lp_max_abs_diff"] < 1e-6
    assert C["lp_run_vs_table_run_max_rel_diff_X"] < 1e-3


def test_closed_form_reserve_matches_linear_program():
    assert C["reserve_closed_form_vs_lp_max_abs_diff"] < 1e-9


def test_time_step_has_converged():
    assert C["dt_0.5_vs_0.125_max_rel_diff_X"] < 0.05


def test_capacity_and_its_limiting_resource():
    nadph = R["level1_capacity"]["NADPH"]["0"]
    assert nadph["network_capacity"] == pytest.approx(0.68, abs=1e-6)
    assert nadph["marginal"]["G6PD capacity (PPP)"] == pytest.approx(2.0, abs=1e-6)


def test_safe_window_closes_at_about_day_58():
    assert round(R["level1_window"]["closing_age_days"]) == 58


def test_reserve_fails_long_before_toxin_control():
    truth = R["level2_truth_design"]
    assert truth["reserve_below_threshold_at_age_d"] == pytest.approx(52, abs=1)
    assert truth["toxin_loses_control_at_age_d"] == pytest.approx(110, abs=1)
    assert truth["reserve_below_threshold_at_age_d"] < truth["toxin_loses_control_at_age_d"]


def test_data_fix_the_product_not_the_parts():
    ident = R["identifiability"]
    assert R["ensemble"]["n_designs"] == 16384
    assert ident["n_accepted"] == 76
    assert ident["corr_log_aE_log_kcat_accepted"] < -0.98
    toxin_only = ident["contraction_toxin_data"]
    assert toxin_only["a_E × kcat"] == pytest.approx(0.94, abs=0.01)
    assert toxin_only["a_E"] == pytest.approx(0.26, abs=0.01)
    assert toxin_only["kcat"] == pytest.approx(0.04, abs=0.01)
    with_enzyme = ident["contraction_plus_enzyme_amount"]
    assert with_enzyme["a_E"] == pytest.approx(0.79, abs=0.02)
    assert with_enzyme["kcat"] == pytest.approx(0.80, abs=0.02)


def test_prediction_is_robust_across_fitted_designs():
    low, _, high = R["identifiability"]["protected_fraction_of_accepted_designs_baseline"]
    assert 0.43 <= low <= high <= 0.45


def test_window_edge_matches_its_closed_form():
    assert C["window_edge_grid_vs_formula_max_abs_diff"] < 0.01  # within one grid step


def test_slide_series_agree_with_the_tables():
    s = R["level1_capacity_series"]
    for key in ("NADPH", "ATP"):
        for i, age in ((0, "0"), (30, "60"), (60, "120")):
            assert s[key][i] == pytest.approx(R["level1_capacity"][key][age]["network_capacity"], abs=1e-6)
    assert len(R["identifiability"]["accepted_log10"]["a_E"]) == R["identifiability"]["n_accepted"]


def test_outputs_exist():
    for name in ("fig1_capacity", "fig2_window", "fig3_data", "fig4_dynamics", "fig5_designs"):
        assert (ROOT / "results" / "figures" / f"{name}.png").stat().st_size > 0
    for name in ("step1_limits", "step2_window", "step4_in_time", "step6_measure"):
        assert (ROOT / "results" / "slides" / f"{name}.png").stat().st_size > 0
    for name in ("d1_network", "d2_coupling", "d3_loop"):
        assert (ROOT / "results" / "method" / f"{name}.png").stat().st_size > 0
    for name in ("report.html", "mathematics.html"):
        assert (ROOT / "docs" / name).stat().st_size > 0
