import pandas as pd
from pathlib import Path

from tools import run_bayesian_analysis as bayes

REPO_ROOT = Path(__file__).resolve().parents[1]
MULTI_STAN_TEXT = (REPO_ROOT / "stan_models" / "VIDRA.stan").read_text(encoding="utf-8")


def test_az_source_mapping_is_retained_for_step3():
    assert bayes.COMBINATION_SLOPE[4] == [1, 2]


def test_gene_output_paths_are_partitioned_by_gene():
    assert bayes.gene_output_partition("/tmp/results", "ENSG1") == "/tmp/results/as_gene=ENSG1"
    assert bayes.gene_output_file("/tmp/results", "ENSG1") == "/tmp/results/as_gene=ENSG1/results.parquet"


def test_list_gene_partition_names_reads_local_partition_dirs(tmp_path):
    (tmp_path / "as_gene=ENSG000001").mkdir()
    (tmp_path / "as_gene=ENSG000002").mkdir()
    (tmp_path / "_temporary").mkdir()

    assert bayes.list_gene_partition_names(str(tmp_path)) == ["ENSG000001", "ENSG000002"]


class _FakeRunset:
    def __init__(self, stdout_files=None, stderr_files=None):
        self.stdout_files = stdout_files or []
        self.stderr_files = stderr_files or []


class _FakeFit:
    def __init__(self, runset):
        self.runset = runset


def test_fit_has_nonconvergence_warning_detects_warning_text(tmp_path):
    stdout = tmp_path / "stdout.txt"
    stdout.write_text("The algorithm may not have converged.\n", encoding="utf-8")
    fit = _FakeFit(_FakeRunset(stdout_files=[str(stdout)]))
    assert bayes._fit_has_nonconvergence_warning(fit) is True


def test_fit_has_nonconvergence_warning_returns_false_when_warning_absent(tmp_path):
    stdout = tmp_path / "stdout.txt"
    stdout.write_text("All good.\n", encoding="utf-8")
    fit = _FakeFit(_FakeRunset(stdout_files=[str(stdout)]))
    assert bayes._fit_has_nonconvergence_warning(fit) is False


def test_output_schema_includes_has_burden():
    assert "has_burden" in bayes.ALL_OUTPUT_COLS


def test_scalar_burden_inputs_return_no_burden_for_zero_bO():
    df = pd.DataFrame({"bO": [0.0, 0.0], "bOse": [2.0, 2.0]})

    assert bayes._extract_scalar_burden_inputs(df) == {
        "has_burden": 0,
        "bO": 0.0,
        "bOse": 2.0,
    }


def test_scalar_burden_inputs_collapse_consistent_repeated_burden():
    df = pd.DataFrame({"bO": [0.75, 0.75, 0.75], "bOse": [0.2, 0.2, 0.2]})

    assert bayes._extract_scalar_burden_inputs(df) == {
        "has_burden": 1,
        "bO": 0.75,
        "bOse": 0.2,
    }


def test_scalar_burden_inputs_keep_non_significant_burden_as_burden():
    df = pd.DataFrame({"bO": [0.2, 0.2], "bOse": [0.2, 0.2]})

    assert bayes._extract_scalar_burden_inputs(df) == {
        "has_burden": 1,
        "bO": 0.2,
        "bOse": 0.2,
    }


def test_scalar_burden_inputs_reject_conflicting_nonzero_burden():
    df = pd.DataFrame({"bO": [0.75, 0.25], "bOse": [0.2, 0.2]})

    try:
        bayes._extract_scalar_burden_inputs(df)
    except ValueError as exc:
        assert "Inconsistent non-zero burden values" in str(exc)
    else:
        raise AssertionError("Expected inconsistent burden values to fail")


def test_as_error_report_sets_has_burden_from_bO():
    df = pd.DataFrame(
        [
            {
                "as_disease": "EFO_0000001",
                "variant": "1_100_A_G",
                "GsourceLab": 0,
                "GqtlLab": 0,
                "bO": 0.75,
            }
        ]
    )

    result = bayes.as_error_report(df, gene="ENSG000001", error_msg="boom")
    assert bool(result.loc[0, "has_burden"]) is True


def test_as_error_report_sets_has_burden_false_without_burden():
    df = pd.DataFrame(
        [
            {
                "as_disease": "EFO_0000001",
                "variant": "1_100_A_G",
                "GsourceLab": 0,
                "GqtlLab": 0,
                "bO": 0.0,
            }
        ]
    )

    result = bayes.as_error_report(df, gene="ENSG000001", error_msg="boom")
    assert bool(result.loc[0, "has_burden"]) is False


def test_apply_annotation_transforms_coerces_invalid_numeric_payloads_to_defaults():
    df = pd.DataFrame(
        [
            {
                "variant": "1_100_A_G",
                "as_blosum62": "not-a-number",
                "as_loftool": None,
                "as_clinicalSignificance": "2",
            }
        ]
    )

    transformed = bayes._apply_annotation_transforms(df.copy())

    assert transformed.loc[0, "as_blosum62"] == 1.0
    assert transformed.loc[0, "as_loftool"] == 0.0
    assert transformed.loc[0, "as_clinicalSignificance"] == 1.0


def test_multi_stan_uses_scalar_burden_likelihood():
    assert "int<lower=0, upper=1> has_burden;" in MULTI_STAN_TEXT
    assert "real bO;" in MULTI_STAN_TEXT
    assert "real<lower=0> bOse;" in MULTI_STAN_TEXT
    assert "bO ~ normal(intercept_random, bOse);" in MULTI_STAN_TEXT
    assert "vector[N] bO" not in MULTI_STAN_TEXT


def test_multi_stan_has_shared_intercept_for_logor_sources():
    assert "real intercept_random;" in MULTI_STAN_TEXT
    assert "intercept_random ~ normal(0, 10);" in MULTI_STAN_TEXT
    assert "intercept_random + slope_random[4] * protein_prior[n]" in MULTI_STAN_TEXT
    assert "intercept_random + protein_prior[n] * slope_random[3]" in MULTI_STAN_TEXT
    assert "slope_random[5] ~ normal( -1, 5 )" not in MULTI_STAN_TEXT
    assert "real intercept_burden_logor;" not in MULTI_STAN_TEXT
    assert "real intercept_coding_logor;" not in MULTI_STAN_TEXT


def test_multi_stan_adapts_shared_intercept_for_clinvar_through_sigmoid():
    assert "intercept_clinvar_baseline" not in MULTI_STAN_TEXT
    assert "inv_logit(intercept_random) + slope_random[5] * protein_prior[n]" in MULTI_STAN_TEXT
