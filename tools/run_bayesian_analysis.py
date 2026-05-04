#!/usr/bin/env python3
from __future__ import annotations

"""VIDRA Step 3: run the Bayesian models for one gene at a time.

This script is the only stage still orchestrated by Nextflow in the hybrid
VIDRA_5 pipeline. Step 1 and Step 2 are run externally and produce:

  - `vidra_analysis_ready/as_gene=<gene>/...`
  - `variant_annotations/annotations.parquet`

Each Step 3 task receives exactly one Step 1 gene partition plus the shared
annotation parquet, runs the Stan models with pandas/CmdStanPy, and writes
only that gene's result partition:

  `vidra_results/as_gene=<gene>/results.parquet`
"""

import atexit
import argparse
import filecmp
import fcntl
import glob
import os
import shutil
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.dataset as ds
import pyarrow.fs as pafs
import pyarrow.parquet as pq
from cmdstanpy import CmdStanModel, set_cmdstan_path

# =============================================================================
# Constants
# =============================================================================

_STAN_NONCONVERGENCE_TEXT = "The algorithm may not have converged."
_STAN_VB_SUMMARY = {
    "total": 0,
    "nonconverged": 0,
    "by_model": {},
}

# Maps (GsourceLab, GqtlLab) pairs → slope_random index in VIDRA.stan
COMBINATION_SLOPE = {1: [0, 0], 2: [0, 1], 3: [3, 2], 4: [1, 2], 5: [2, 2]}

ANNOTATION_DEFAULTS = {
    "as_blosum62": 1.0,
    "as_conservation": 1.0,
    "as_sift": 1.0,
    "as_polyphen": 1.0,
    "as_cadd": 1.0,
    "as_alphamissense": 1.0,
    "as_primateai": 0.0,
    "as_revel": 1.0,
    "as_clinicalSignificance": 0.0,
    "as_loftool": 0.0,
    "as_plddt": 1.0,
    "as_consequence": 1.0,
    "as_foldx": 1.0,
    "foldxDdq_raw": float("nan"),
}

ANNOTATION_COLS = [
    "as_blosum62",
    "as_conservation",
    "as_sift",
    "as_polyphen",
    "as_cadd",
    "as_alphamissense",
    "as_revel",
    "as_primateai",
    "as_loftool",
    "as_plddt",
    "as_consequence",
    "as_clinicalSignificance",
    "foldxDdq_raw",
    "most_severe_consequence",
]

LOF_CONSEQUENCES = {"stop_gained", "start_lost"}

GENE_MEAN_IMPUTE_COLS = [
    "as_revel",
    "as_alphamissense",
    "as_sift",
    "as_polyphen",
    "as_foldx",
    "as_cadd",
    "as_conservation",
]

INVERSION_COLS = [
    "as_revel",
    "as_polyphen",
    "as_cadd",
    "as_alphamissense",
    "as_plddt",
]

FULLRANK_MAX_N = 200
ADVI_MAX_ITER = 10000
ADVI_GRAD_SAMPLES = 20
ADVI_DRAWS = 1000

POSTERIOR_COLS = [
    "mean",
    "median",
    "pct_1",
    "pct_2_5",
    "pct_5",
    "pct_10",
    "pct_25",
    "pct_40",
    "pct_50",
    "pct_60",
    "pct_75",
    "pct_90",
    "pct_95",
    "pct_97_5",
    "pct_99",
    "pp_slope_pos",
    "pp_slope_neg",
]

META_COLS = ["gene", "as_disease", "parameter", "model", "n_variants", "source", "qtl"]
BURDEN_COLS = ["has_burden"]
ALL_OUTPUT_COLS = META_COLS + POSTERIOR_COLS + BURDEN_COLS

MODELS = {}


# =============================================================================
# Generic path helpers
# =============================================================================

def join_root_path(root: str, relative_path: str) -> str:
    base = str(root).rstrip("/")
    rel = str(relative_path).lstrip("/")
    return f"{base}/{rel}" if base else rel


def is_gs_path(path: str | Path | None) -> bool:
    return str(path or "").startswith("gs://")


def _strip_gs_prefix(path: str) -> str:
    return str(path)[len("gs://"):]


def _get_gcsfs():
    try:
        import gcsfs
    except ImportError as exc:
        raise RuntimeError("Working with gs:// paths requires gcsfs.") from exc
    return gcsfs.GCSFileSystem()


def path_exists(path: str | Path) -> bool:
    if is_gs_path(path):
        fs = _get_gcsfs()
        return fs.exists(_strip_gs_prefix(str(path)))
    return Path(path).exists()


def remove_path(path: str | Path) -> None:
    path = str(path)
    if not path_exists(path):
        return
    if is_gs_path(path):
        fs = _get_gcsfs()
        fs.rm(_strip_gs_prefix(path), recursive=True)
        return

    local = Path(path)
    if local.is_dir():
        shutil.rmtree(local)
    else:
        local.unlink()


def _get_arrow_filesystem(path_like: str):
    if is_gs_path(path_like):
        fs = _get_gcsfs()
        return pafs.PyFileSystem(pafs.FSSpecHandler(fs)), _strip_gs_prefix(path_like)
    return None, str(path_like)


def open_parquet_dataset(path_like: str):
    filesystem, stripped = _get_arrow_filesystem(path_like)
    if filesystem is not None:
        return ds.dataset(stripped, filesystem=filesystem, format="parquet")
    return ds.dataset(path_like, format="parquet")


def read_parquet_dataframe(path_like: str, columns: list[str] | None = None) -> pd.DataFrame:
    dataset = open_parquet_dataset(path_like)
    table = dataset.to_table(columns=columns) if columns else dataset.to_table()
    return table.to_pandas()


def _write_table_to_path(table: pa.Table, output_path: str) -> None:
    if is_gs_path(output_path):
        fs = _get_gcsfs()
        with fs.open(_strip_gs_prefix(output_path), "wb") as handle:
            pq.write_table(table, handle)
        return

    out = Path(output_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    pq.write_table(table, out)


def _write_text_file(path: str, text: str) -> None:
    if is_gs_path(path):
        fs = _get_gcsfs()
        with fs.open(_strip_gs_prefix(path), "w", encoding="utf-8") as handle:
            handle.write(text)
        return

    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(text, encoding="utf-8")


def list_gene_partition_names(root_dir: str) -> list[str]:
    if is_gs_path(root_dir):
        fs = _get_gcsfs()
        prefix = _strip_gs_prefix(root_dir.rstrip("/"))
        if not fs.exists(prefix):
            return []
        entries = fs.ls(prefix, detail=True)
        genes = []
        for entry in entries:
            name = str(entry.get("name", "")).rstrip("/").split("/")[-1]
            entry_type = str(entry.get("type", ""))
            if name.startswith("as_gene=") and entry_type in {"directory", "dir"}:
                genes.append(name.replace("as_gene=", "", 1))
        return sorted(set(genes))

    root = Path(root_dir)
    if not root.exists():
        return []
    return sorted(
        path.name.replace("as_gene=", "", 1)
        for path in root.iterdir()
        if path.is_dir() and path.name.startswith("as_gene=")
    )


def gene_output_partition(output_dir: str, gene: str) -> str:
    return join_root_path(output_dir, f"as_gene={gene}")


def gene_output_file(output_dir: str, gene: str) -> str:
    return join_root_path(gene_output_partition(output_dir, gene), "results.parquet")


# =============================================================================
# Stan model setup
# =============================================================================

def _fit_has_nonconvergence_warning(fit) -> bool:
    runset = getattr(fit, "runset", None)
    if runset is None:
        return False

    for attr in ("stdout_files", "stderr_files"):
        for path in getattr(runset, attr, []) or []:
            try:
                text = Path(path).read_text(encoding="utf-8", errors="ignore")
            except OSError:
                continue
            if _STAN_NONCONVERGENCE_TEXT in text:
                return True
    return False


def _record_stan_vb_status(model_name: str, fit) -> None:
    summary = _STAN_VB_SUMMARY
    summary["total"] += 1

    model_stats = summary["by_model"].setdefault(
        model_name,
        {"total": 0, "nonconverged": 0},
    )
    model_stats["total"] += 1

    if _fit_has_nonconvergence_warning(fit):
        summary["nonconverged"] += 1
        model_stats["nonconverged"] += 1


def _emit_stan_vb_summary() -> None:
    total = int(_STAN_VB_SUMMARY["total"] or 0)
    if total <= 0:
        return

    nonconverged = int(_STAN_VB_SUMMARY["nonconverged"] or 0)
    converged = total - nonconverged
    print(
        f"Stan VB convergence summary: total={total}, "
        f"converged={converged}, nonconverged={nonconverged}"
    )
    for model_name in sorted(_STAN_VB_SUMMARY["by_model"]):
        stats = _STAN_VB_SUMMARY["by_model"][model_name]
        print(
            f"  {model_name}: total={int(stats['total'] or 0)}, "
            f"nonconverged={int(stats['nonconverged'] or 0)}"
        )


atexit.register(_emit_stan_vb_summary)


def _stan_cache_dir():
    return os.environ.get("VIDRA_STAN_CACHE_DIR", "/tmp/vidra_stan_cache")


def _copy_if_changed(src, dst):
    if not os.path.isfile(dst):
        shutil.copy2(src, dst)
        return True
    try:
        same = filecmp.cmp(src, dst, shallow=False)
    except OSError:
        same = False
    if not same:
        shutil.copy2(src, dst)
        return True
    return False


def _configure_cmdstan_path():
    candidates = []
    env_path = os.environ.get("CMDSTAN", "").strip()
    if env_path:
        candidates.append(env_path)

    candidates.extend(["/cmdstan", "/opt/cmdstan"])
    candidates.extend(sorted(glob.glob("/cmdstan/cmdstan-*"), reverse=True))
    candidates.extend(sorted(glob.glob("/opt/cmdstan/cmdstan-*"), reverse=True))
    candidates.extend(sorted(glob.glob("/root/.cmdstan/cmdstan-*"), reverse=True))
    candidates.extend(sorted(glob.glob("/home/*/.cmdstan/cmdstan-*"), reverse=True))

    seen = set()
    for path in candidates:
        if not path or path in seen:
            continue
        seen.add(path)
        if not os.path.isdir(os.path.join(path, "bin")):
            continue
        set_cmdstan_path(path)
        os.environ["CMDSTAN"] = path
        print(f"Using CmdStan path: {path}")
        return path

    raise ValueError("No valid CmdStan installation found.")


def _ensure_cached_compiled_models(path_multi, path_single):
    cache_dir = _stan_cache_dir()
    os.makedirs(cache_dir, exist_ok=True)

    cached_multi = os.path.join(cache_dir, "VIDRA.stan")
    cached_single = os.path.join(cache_dir, "VIDRA_single_variant.stan")
    cached_multi_exe = os.path.join(cache_dir, "VIDRA")
    cached_single_exe = os.path.join(cache_dir, "VIDRA_single_variant")
    lock_path = os.path.join(cache_dir, ".compile.lock")

    with open(lock_path, "w", encoding="utf-8") as lockf:
        fcntl.flock(lockf, fcntl.LOCK_EX)

        multi_updated = _copy_if_changed(path_multi, cached_multi)
        single_updated = _copy_if_changed(path_single, cached_single)

        if multi_updated and os.path.exists(cached_multi_exe):
            os.remove(cached_multi_exe)
        if single_updated and os.path.exists(cached_single_exe):
            os.remove(cached_single_exe)

        multi_ready = os.path.isfile(cached_multi_exe) and os.access(cached_multi_exe, os.X_OK)
        single_ready = os.path.isfile(cached_single_exe) and os.access(cached_single_exe, os.X_OK)

        if not multi_ready:
            CmdStanModel(stan_file=cached_multi)
        if not single_ready:
            CmdStanModel(stan_file=cached_single)

        fcntl.flock(lockf, fcntl.LOCK_UN)

    return {
        "multi_stan": cached_multi,
        "single_stan": cached_single,
        "multi_exe": cached_multi_exe,
        "single_exe": cached_single_exe,
    }


def _find_stan_model_artifacts():
    candidate_dirs = [_stan_cache_dir()]
    for env_var in ("VIDRA_STAN_MODELS_DIR", "STAN_MODELS_DIR"):
        value = os.environ.get(env_var)
        if value:
            candidate_dirs.append(value)

    candidate_dirs.extend([
        "/opt/vidra/stan_models",
        "/container/stan_models",
        str(Path(__file__).resolve().parent.parent / "stan_models"),
    ])

    first_stan_only = None
    seen = set()
    for base in candidate_dirs:
        if not base or base in seen:
            continue
        seen.add(base)

        multi_stan = os.path.join(base, "VIDRA.stan")
        single_stan = os.path.join(base, "VIDRA_single_variant.stan")
        if not (os.path.isfile(multi_stan) and os.path.isfile(single_stan)):
            continue

        multi_exe = os.path.join(base, "VIDRA")
        single_exe = os.path.join(base, "VIDRA_single_variant")
        has_exe = (
            os.path.isfile(multi_exe)
            and os.path.isfile(single_exe)
            and os.access(multi_exe, os.X_OK)
            and os.access(single_exe, os.X_OK)
        )

        if has_exe:
            return {
                "mode": "exe",
                "multi_stan": multi_stan,
                "single_stan": single_stan,
                "multi_exe": multi_exe,
                "single_exe": single_exe,
            }

        if first_stan_only is None:
            first_stan_only = {
                "mode": "stan",
                "multi_stan": multi_stan,
                "single_stan": single_stan,
            }

    if first_stan_only is not None:
        return first_stan_only

    raise ValueError("Could not locate Stan model sources.")


def get_models():
    if MODELS:
        return MODELS

    artifacts = _find_stan_model_artifacts()
    if artifacts["mode"] == "exe":
        MODELS["multi"] = CmdStanModel(
            stan_file=artifacts["multi_stan"],
            exe_file=artifacts["multi_exe"],
        )
        MODELS["single"] = CmdStanModel(
            stan_file=artifacts["single_stan"],
            exe_file=artifacts["single_exe"],
        )
        return MODELS

    _configure_cmdstan_path()
    cached = _ensure_cached_compiled_models(
        artifacts["multi_stan"],
        artifacts["single_stan"],
    )
    MODELS["multi"] = CmdStanModel(
        stan_file=cached["multi_stan"],
        exe_file=cached["multi_exe"],
    )
    MODELS["single"] = CmdStanModel(
        stan_file=cached["single_stan"],
        exe_file=cached["single_exe"],
    )
    return MODELS


# =============================================================================
# Posterior helpers
# =============================================================================

def clean_posteriorForAs(serie):
    s = np.asarray(serie)
    return {
        "mean": float(np.mean(s)),
        "median": float(np.median(s)),
        "pct_1": float(np.percentile(s, 1)),
        "pct_2_5": float(np.percentile(s, 2.5)),
        "pct_5": float(np.percentile(s, 5)),
        "pct_10": float(np.percentile(s, 10)),
        "pct_25": float(np.percentile(s, 25)),
        "pct_40": float(np.percentile(s, 40)),
        "pct_50": float(np.percentile(s, 50)),
        "pct_60": float(np.percentile(s, 60)),
        "pct_75": float(np.percentile(s, 75)),
        "pct_90": float(np.percentile(s, 90)),
        "pct_95": float(np.percentile(s, 95)),
        "pct_97_5": float(np.percentile(s, 97.5)),
        "pct_99": float(np.percentile(s, 99)),
        "pp_slope_pos": float((s > 0).mean()),
        "pp_slope_neg": float((s < 0).mean()),
    }


def _get_slope_key(combination_slope, value_to_find):
    for key, value in combination_slope.items():
        if value == value_to_find:
            return key
    return None


def _safe_float(val, default=0.0):
    try:
        f = float(val)
        return default if not np.isfinite(f) else f
    except (TypeError, ValueError):
        return default


def _fill_annotation(df, col, default=0.0):
    if col not in df.columns:
        df[col] = default
    else:
        df[col] = df[col].fillna(default)
    return df


# =============================================================================
# Stan runners
# =============================================================================

def _get_scratch_dir():
    scratch = os.environ.get("SPARK_LOCAL_DIRS", "/tmp").split(",")[0]
    os.environ["TMPDIR"] = scratch
    return scratch


def _extract_scalar_burden_inputs(df: pd.DataFrame) -> dict:
    """
    Convert Step 1 burden columns into the scalar Stan contract.

    Step 1 joins the same gene/disease burden estimate onto every variant row.
    Passing it as a vector would count the same burden evidence once per variant,
    so Step 3 collapses it to one scalar and rejects conflicting non-zero values.
    """
    if "bO" not in df.columns:
        return {
            "has_burden": 0,
            "bO": 0.0,
            "bOse": 2.0,
        }

    burden = pd.DataFrame({
        "bO": pd.to_numeric(df["bO"], errors="coerce"),
        "bOse": pd.to_numeric(df["bOse"], errors="coerce") if "bOse" in df.columns else 2.0,
    })
    burden = burden.replace([np.inf, -np.inf], np.nan)
    burden = burden[burden["bO"].notna() & (burden["bO"].abs() > 1e-12)].copy()

    if burden.empty:
        return {
            "has_burden": 0,
            "bO": 0.0,
            "bOse": 2.0,
        }

    burden["bOse"] = burden["bOse"].where(burden["bOse"].notna() & (burden["bOse"] > 0.0), 2.0)
    first_logor = float(burden["bO"].iloc[0])
    first_se = max(float(burden["bOse"].iloc[0]), 1e-6)

    inconsistent = (
        ~np.isclose(burden["bO"].to_numpy(dtype=float), first_logor, rtol=1e-9, atol=1e-12)
        | ~np.isclose(burden["bOse"].to_numpy(dtype=float), first_se, rtol=1e-9, atol=1e-12)
    )
    if bool(inconsistent.any()):
        sample = burden.drop_duplicates().head(5).to_dict(orient="records")
        raise ValueError(f"Inconsistent non-zero burden values within one gene/disease block: {sample}")

    return {
        "has_burden": 1,
        "bO": first_logor,
        "bOse": first_se,
    }


def AS_singleVars(df, gene, h1=0.1):
    model = get_models()["single"]
    if len(df) == 0:
        return None

    row = df.iloc[0]
    df_dict = {
        "h1": float(h1),
        "N": len(df),
        "numG1": _safe_float(row["GsourceLab"]),
        "numG2": _safe_float(row["GqtlLab"]),
        "xc": _safe_float(row["xc"], 0.0),
        "xcse": max(_safe_float(row["xcse"], 0.1), 1e-6),
        "yOR": _safe_float(row["yc"]),
        "yORse": max(_safe_float(row["ycse"], 0.14), 1e-6),
        "as_blosum62": _safe_float(row.get("as_blosum62", 1.0)),
        "as_conservation": _safe_float(row.get("as_conservation", 1.0)),
        "as_sift": _safe_float(row.get("as_sift", 1.0)),
        "as_polyphen": _safe_float(row.get("as_polyphen", 1.0)),
        "as_clinicalSignificance": _safe_float(row.get("as_clinicalSignificance", 0.0)),
        "as_cadd": _safe_float(row.get("as_cadd", 1.0)),
        "as_alphamissense": _safe_float(row.get("as_alphamissense", 1.0)),
        "as_revel": _safe_float(row.get("as_revel", 1.0)),
        "as_primateai": _safe_float(row.get("as_primateai", 0.0)),
    }

    # Burden-informed intercepts are multi-variant-only in the current model.
    # The single-variant Stan model does not consume bO/bOse and remains a
    # legacy fallback until a dedicated single-variant intercept redesign.
    scratch = _get_scratch_dir()
    with tempfile.TemporaryDirectory(prefix="stan_sv_", dir=scratch) as tmpdir:
        fit = model.variational(
            data=df_dict,
            seed=412,
            algorithm="fullrank",
            iter=ADVI_MAX_ITER,
            grad_samples=ADVI_GRAD_SAMPLES,
            draws=ADVI_DRAWS,
            require_converged=False,
            show_console=False,
            refresh=1000,
            output_dir=tmpdir,
        )
        _record_stan_vb_status("single_variant", fit)
        slope_posteriors = fit.stan_variable("slope", mean=False)

    res = clean_posteriorForAs(slope_posteriors)
    res["gene"] = gene
    res["as_disease"] = str(df["as_disease"].iloc[0])
    res["n_variants"] = int(df["variant"].nunique())
    res["source"] = str(int(df_dict["numG1"]))
    res["qtl"] = str(int(df_dict["numG2"]))
    res["model"] = "single_variant"
    res["parameter"] = "slope"
    res["has_burden"] = False
    return pd.DataFrame([res])


def AS_multiVars(df, gene):
    model = get_models()["multi"]
    nu = max(int(df["variant"].nunique()) - 1, 1)
    burden_inputs = _extract_scalar_burden_inputs(df)
    df_dict = {
        "nu": nu,
        "N": len(df),
        **burden_inputs,
        "numG1": df["GsourceLab"].astype(float).tolist(),
        "numG2": df["GqtlLab"].astype(float).tolist(),
        "xc": df["xc"].astype(float).replace([np.inf, -np.inf], np.nan).fillna(1.0).tolist(),
        "xcse": df["xcse"].astype(float).replace([np.inf, -np.inf], np.nan).fillna(1.0).clip(lower=1e-6).tolist(),
        "yOR": df["yc"].astype(float).replace([np.inf, -np.inf], np.nan).fillna(0.0).tolist(),
        "yORse": df["ycse"].astype(float).replace([np.inf, -np.inf], np.nan).fillna(0.14).clip(lower=1e-6).tolist(),
        "as_blosum62": df["as_blosum62"].astype(float).tolist(),
        "as_conservation": df["as_conservation"].astype(float).tolist(),
        "as_sift": df["as_sift"].astype(float).tolist(),
        "as_polyphen": df["as_polyphen"].astype(float).tolist(),
        "as_cadd": df["as_cadd"].astype(float).tolist(),
        "as_alphamissense": df["as_alphamissense"].astype(float).tolist(),
        "as_revel": df["as_revel"].astype(float).tolist(),
        "as_foldx": df["as_foldx"].astype(float).tolist(),
        "as_consequence": df["as_consequence"].astype(float).tolist(),
        "as_clinicalSignificance": df["as_clinicalSignificance"].astype(float).tolist(),
        "as_primateai": df["as_primateai"].astype(float).tolist(),
    }

    scratch = _get_scratch_dir()
    stan_tmpdir = tempfile.mkdtemp(prefix="stan_mv_", dir=scratch)
    n_unique_variants = int(df["variant"].nunique())
    algorithm = "fullrank" if n_unique_variants <= FULLRANK_MAX_N else "meanfield"

    fit = model.variational(
        data=df_dict,
        seed=412,
        algorithm=algorithm,
        iter=ADVI_MAX_ITER,
        grad_samples=ADVI_GRAD_SAMPLES,
        draws=ADVI_DRAWS,
        require_converged=False,
        show_console=False,
        refresh=1000,
        output_dir=stan_tmpdir,
    )
    _record_stan_vb_status(f"multiple_variant_{algorithm}", fit)

    posteriors = fit.variational_sample_pd
    slope_samples = posteriors["slope"].values
    intercept_samples = posteriors["intercept_random"].values
    sr_cols = [c for c in posteriors.columns if c.startswith("slope_random[")]
    slope_random_samples = posteriors[sorted(sr_cols)].values

    param_samples = {
        "slope": slope_samples,
        "intercept": intercept_samples,
    }
    for i in range(slope_random_samples.shape[1]):
        param_samples[f"slope_random[{i+1}]"] = slope_random_samples[:, i]

    combination_observed = [
        list(pair)
        for pair in df[["GsourceLab", "GqtlLab"]].drop_duplicates().itertuples(index=False)
    ]
    slope_to_meta_analyse = []
    for combination in combination_observed:
        key = _get_slope_key(COMBINATION_SLOPE, combination)
        if key is not None:
            slope_to_meta_analyse.append(key)

    list_posteriors = []
    for slope_idx in slope_to_meta_analyse:
        col_name = f"slope_random[{slope_idx}]"
        if col_name in posteriors.columns:
            list_posteriors += posteriors[col_name].tolist()

    rows = []
    for param_name, samples in param_samples.items():
        stats = clean_posteriorForAs(samples)
        stats["parameter"] = param_name
        rows.append(stats)

    meta_stats = clean_posteriorForAs(pd.Series(list_posteriors or slope_samples))
    meta_stats["parameter"] = "meta_slope"
    rows.append(meta_stats)

    output = pd.DataFrame(rows)
    output["gene"] = gene
    output["as_disease"] = str(df["as_disease"].iloc[0])
    output["n_variants"] = int(df["variant"].nunique())
    output["source"] = str(sorted(set(int(x) for x in df["GsourceLab"])))
    output["qtl"] = str(sorted(set(int(x) for x in df["GqtlLab"])))
    output["model"] = f"multiple_variant_{algorithm}"
    output["has_burden"] = bool(burden_inputs["has_burden"])

    shutil.rmtree(stan_tmpdir, ignore_errors=True)
    return output


def as_error_report(df, gene, error_msg=""):
    try:
        has_burden = bool(_extract_scalar_burden_inputs(df)["has_burden"])
    except Exception:
        bO_vals = pd.to_numeric(df["bO"], errors="coerce") if "bO" in df.columns else pd.Series(dtype=float)
        has_burden = bool((bO_vals.replace([np.inf, -np.inf], np.nan).abs() > 1e-12).any())
    row = {col: None for col in ALL_OUTPUT_COLS}
    row.update(
        {
            "gene": gene,
            "as_disease": str(df["as_disease"].iloc[0]),
            "n_variants": int(df["variant"].nunique()),
            "model": "error_report",
            "parameter": str(error_msg)[:500] if error_msg else "error",
            "source": str(sorted(set(int(x) for x in df["GsourceLab"]))),
            "qtl": str(sorted(set(int(x) for x in df["GqtlLab"]))),
            "has_burden": has_burden,
        }
    )
    return pd.DataFrame([row], columns=ALL_OUTPUT_COLS)


def fitASmodels(disease_df, gene, h1=0.1):
    n_unique_variants = disease_df["variant"].nunique()
    try:
        result = AS_multiVars(disease_df, gene) if n_unique_variants > 1 else AS_singleVars(disease_df, gene, h1)
        if result is None or result.empty:
            return as_error_report(disease_df, gene, "empty result")
        return result
    except Exception as exc:  # pragma: no cover - exercised through error rows
        import traceback

        error_msg = f"{type(exc).__name__}: {exc} || {traceback.format_exc()[-400:]}"
        return as_error_report(disease_df, gene, error_msg)


# =============================================================================
# Per-gene preprocessing and disease fitting
# =============================================================================

def preprocess_gene(gene_df):
    if gene_df.empty:
        return gene_df

    if "foldxDdq_raw" in gene_df.columns:
        foldx_raw = pd.to_numeric(gene_df["foldxDdq_raw"], errors="coerce")
        max_score = foldx_raw.max()
        if pd.notna(max_score) and max_score != 0:
            inverted = max_score - foldx_raw
            max_inv = inverted.max()
            gene_df["as_foldx"] = (1.0 - (inverted / max_inv)).fillna(1.0) if max_inv > 0 else 1.0
        else:
            gene_df["as_foldx"] = 1.0

    missense_cols = set(GENE_MEAN_IMPUTE_COLS)
    for col, default in ANNOTATION_DEFAULTS.items():
        if col in missense_cols:
            if col not in gene_df.columns:
                gene_df[col] = float("nan")
        else:
            _fill_annotation(gene_df, col, default)

    for col in GENE_MEAN_IMPUTE_COLS:
        if col in gene_df.columns:
            gene_mean = gene_df[col].mean(skipna=True)
            if pd.notna(gene_mean):
                gene_df[col] = gene_df[col].fillna(gene_mean)

    if "most_severe_consequence" in gene_df.columns:
        lof_mask = gene_df["most_severe_consequence"].isin(LOF_CONSEQUENCES)
        if lof_mask.any():
            for col in GENE_MEAN_IMPUTE_COLS:
                if col in gene_df.columns:
                    gene_df.loc[lof_mask, col] = 0.0

    for col in GENE_MEAN_IMPUTE_COLS:
        if col in gene_df.columns:
            gene_df[col] = gene_df[col].fillna(ANNOTATION_DEFAULTS.get(col, 1.0))

    numeric_cols = ["xc", "xcse", "yc", "ycse", "bO", "bOse", "GsourceLab", "GqtlLab"]
    for col in numeric_cols:
        if col in gene_df.columns:
            gene_df[col] = gene_df[col].fillna(0.0)

    gene_df["GsourceLab"] = gene_df["GsourceLab"].astype(int)
    gene_df["GqtlLab"] = gene_df["GqtlLab"].astype(int)

    # Legacy deduplication from the older Step 3 path.
    # Step 1 now performs source-aware deduplication for common, coding, ClinVar,
    # and AZ rare-variant rows, so Step 3 should consume the prepared rows as-is
    # instead of dropping duplicates again here.
    #
    # gene_df = gene_df.sort_values("yc")
    # gene_df = gene_df.drop_duplicates(
    #     subset=["variant", "as_disease", "GsourceLab", "GqtlLab"],
    #     keep="last",
    # )
    gene_df = gene_df.groupby(["as_disease", "GsourceLab", "GqtlLab"], group_keys=False).filter(
        lambda x: not ((len(x["variant"]) == 1) and (x["GsourceLab"] == 3).all())
    )
    return gene_df


def process_disease(disease_df, h1=0.1):
    if disease_df.empty:
        return pd.DataFrame(columns=ALL_OUTPUT_COLS)

    gene = str(disease_df["as_gene"].iloc[0])
    try:
        get_models()
    except Exception as exc:  # pragma: no cover - exercised through error rows
        import traceback

        err = f"ModelLoad|{type(exc).__name__}: {exc}|{traceback.format_exc()[-500:]}"
        return as_error_report(disease_df.iloc[[0]], gene, err)

    result = fitASmodels(disease_df, gene, h1)
    if result is None or result.empty:
        return pd.DataFrame(columns=ALL_OUTPUT_COLS)

    for col in ALL_OUTPUT_COLS:
        if col not in result.columns:
            result[col] = None
    return result[ALL_OUTPUT_COLS]


def _merge_annotations(analysis_df: pd.DataFrame, annotations_df: pd.DataFrame) -> pd.DataFrame:
    df = analysis_df.copy()
    cols_to_drop = [c for c in ANNOTATION_COLS if c in df.columns and c != "as_clinicalSignificance"]
    if "as_clinicalSignificance" in df.columns:
        df = df.rename(columns={"as_clinicalSignificance": "_base_as_clinicalSignificance"})
    if cols_to_drop:
        df = df.drop(columns=cols_to_drop)

    merged = df.merge(annotations_df, on="variant", how="left")
    if "_base_as_clinicalSignificance" in merged.columns:
        base = merged["_base_as_clinicalSignificance"].fillna(0.0)
        current = merged["as_clinicalSignificance"].fillna(0.0)
        merged["as_clinicalSignificance"] = np.where(base > 0.0, base, current)
        merged = merged.drop(columns=["_base_as_clinicalSignificance"])
    return merged


def _apply_annotation_transforms(df: pd.DataFrame) -> pd.DataFrame:
    for col in INVERSION_COLS:
        if col in df.columns:
            df[col] = 1.0 - pd.to_numeric(df[col], errors="coerce")

    if "as_consequence" in df.columns:
        df["as_consequence"] = 1.0 - pd.to_numeric(df["as_consequence"], errors="coerce") / 12.0

    if "as_clinicalSignificance" in df.columns:
        cs_max = pd.to_numeric(df["as_clinicalSignificance"], errors="coerce").max()
        if pd.notna(cs_max) and cs_max > 1.0:
            df["as_clinicalSignificance"] = pd.to_numeric(df["as_clinicalSignificance"], errors="coerce") / float(cs_max)

    deferred = set(GENE_MEAN_IMPUTE_COLS)
    for col, default in ANNOTATION_DEFAULTS.items():
        if col in df.columns and col not in deferred and not (isinstance(default, float) and np.isnan(default)):
            # Step 2 encodes the supported annotation fields as numerics. Coerce
            # any unexpected string payloads to NaN so we can fall back to the
            # same per-column defaults instead of depending on removed pandas
            # "errors='ignore'" behavior.
            df[col] = pd.to_numeric(df[col], errors="coerce")
            df[col] = df[col].fillna(default)

    return df


def run_gene_analysis(*, analysis_dir: str, annotations_path: str, output_dir: str, gene: str, h1: float) -> pd.DataFrame:
    print("--- VIDRA Bayesian Analysis Per Gene ---")
    print(f"Gene:        {gene}")
    print(f"Analysis:    {analysis_dir}")
    print(f"Annotations: {annotations_path}")
    print(f"Output dir:  {output_dir}")
    print(f"h1:          {h1}")

    analysis_df = read_parquet_dataframe(analysis_dir)
    if analysis_df.empty:
        remove_path(gene_output_partition(output_dir, gene))
        return pd.DataFrame(columns=ALL_OUTPUT_COLS)

    if "as_gene" not in analysis_df.columns:
        analysis_df["as_gene"] = gene

    annotations_df = read_parquet_dataframe(annotations_path)
    merged = _merge_annotations(analysis_df, annotations_df)
    merged = _apply_annotation_transforms(merged)
    preprocessed = preprocess_gene(merged)

    if preprocessed.empty:
        remove_path(gene_output_partition(output_dir, gene))
        return pd.DataFrame(columns=ALL_OUTPUT_COLS)

    disease_frames = []
    for _, disease_df in preprocessed.groupby("as_disease", sort=False):
        result = process_disease(disease_df.copy(), h1=h1)
        if result is not None and not result.empty:
            disease_frames.append(result)

    if not disease_frames:
        remove_path(gene_output_partition(output_dir, gene))
        return pd.DataFrame(columns=ALL_OUTPUT_COLS)

    final_results = pd.concat(disease_frames, ignore_index=True)
    final_results = final_results[ALL_OUTPUT_COLS]
    output_partition = gene_output_partition(output_dir, gene)
    output_file = gene_output_file(output_dir, gene)
    remove_path(output_partition)
    _write_table_to_path(pa.Table.from_pandas(final_results, preserve_index=False), output_file)
    print(f"Wrote {len(final_results)} result rows to {output_file}")
    return final_results


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run VIDRA Step 3 for one gene.")
    parser.add_argument("--analysis_dir", required=True, help="Path to one as_gene=<gene> Step 1 partition")
    annotations_group = parser.add_mutually_exclusive_group(required=True)
    annotations_group.add_argument("--annotations_file", help="Path to variant_annotations/annotations.parquet")
    annotations_group.add_argument("--annotations_dir", help="Directory containing the annotation parquet")
    parser.add_argument("--output_dir", required=True, help="Result root directory containing as_gene=<gene> partitions")
    parser.add_argument("--gene", required=True, help="ENSG identifier for the current task")
    parser.add_argument("--h1", type=float, default=0.1, help="Hyperparameter passed through to the Stan models")
    return parser.parse_args()


def main(args: argparse.Namespace) -> None:
    annotations_path = args.annotations_file or args.annotations_dir
    run_gene_analysis(
        analysis_dir=args.analysis_dir,
        annotations_path=annotations_path,
        output_dir=args.output_dir,
        gene=args.gene,
        h1=float(args.h1),
    )


if __name__ == "__main__":
    main(parse_args())
