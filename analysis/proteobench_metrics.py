"""ProteoBench quant-accuracy metrics for our 20 benchmark analyses.

This module wraps ProteoBench's `QuantScoresHYE` + `QuantDatapointHYE`
pipeline so we can compute the same per-replicate-threshold metrics
ProteoBench publishes (median_abs_epsilon_global, ROC AUC, per-species
log2 fold-change, CV_median, ...) on **our** `diann_report.pr_matrix.tsv`
files, without having to either resubmit to ProteoBench or recompute the
metric formulas ourselves.

Two adapters are needed because quantmsdiann's pr_matrix.tsv format
differs from the DIA-NN main report ProteoBench's parsers expect:

1. **Wide → long melt.** pr_matrix.tsv stores one column per run; we
   melt to a long-format DataFrame whose columns match what
   `ParseSettingsBuilder` expects (`Modified.Sequence`, `Protein.Ids`,
   `Precursor.Charge`, `Run`, `Precursor.Normalised`).
2. **Species re-annotation.** quantmsdiann strips species suffixes from
   the Protein.Ids column (bare `Q96P70` instead of `Q96P70;Q96P70_HUMAN`).
   ProteoBench's species detection works by substring match on the
   Protein.Ids string for `_HUMAN` / `_YEAST` / `_ECOLI`. We rebuild a
   UniProt-accession → species map once from three SwissProt FASTA
   streams cached under `data/quantmsdiann_benchmarks/uniprot/` and use
   it to add the species suffix to each accession before passing to
   ProteoBench.

Per-analysis results are cached as a single JSON under
`data/quantmsdiann_benchmarks/proteobench_metrics/<dataset>_<version>.json`
so the figure rebuild is offline thereafter. 20 invocations total when
the cache is cold.
"""
from __future__ import annotations

import json
import re
from functools import lru_cache
from pathlib import Path
from typing import Iterable

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = REPO_ROOT / "data" / "quantmsdiann_benchmarks"
UNIPROT_DIR = DATA_DIR / "uniprot"
METRICS_CACHE_DIR = DATA_DIR / "proteobench_metrics"

# Per-dataset → ProteoBench module mapping. Matches the spec mapping in
# the benchmarks design doc.
DATASET_TO_MODULE = {
    "PXD049412":             "quant_lfq_DIA_ion_singlecell",
    "PXD062685":             "quant_lfq_DIA_ion_diaPASEF",
    "PXD070049":             "quant_lfq_DIA_ion_ZenoTOF",
    "ProteoBench_Module_7":  "quant_lfq_DIA_ion_Astral",
}

# Per-module → ProteoBench parse-settings subdir (relative to the
# proteobench package's `io/parsing/io_parse_settings/` root).
MODULE_TO_PARSE_DIR = {
    "quant_lfq_DIA_ion_singlecell": ("Quant", "lfq", "DIA", "ion", "singlecell"),
    "quant_lfq_DIA_ion_diaPASEF":   ("Quant", "lfq", "DIA", "ion", "diaPASEF"),
    "quant_lfq_DIA_ion_ZenoTOF":    ("Quant", "lfq", "DIA", "ion", "ZenoTOF"),
    "quant_lfq_DIA_ion_Astral":     ("Quant", "lfq", "DIA", "ion", "Astral"),
}

# Columns that are NOT per-run intensity columns in a DIA-NN pr_matrix.tsv.
# Everything else is treated as a run column when melting.
PR_MATRIX_META_COLS = {
    "Protein.Group", "Protein.Ids", "Protein.Names", "Genes",
    "First.Protein.Description", "Proteotypic",
    "Stripped.Sequence", "Modified.Sequence",
    "Precursor.Charge", "Precursor.Id",
}

# Recognised raw-file extensions appended by DIA-NN to the run column
# headers. Stripped to derive the bare run identifier that matches
# ProteoBench's condition_mapper keys.
_RUN_EXT_RE = re.compile(r"\.(raw|d|mzML)$", re.IGNORECASE)


# ---------------------------------------------------------------------------
# UniProt species map
# ---------------------------------------------------------------------------

_SPECIES_FASTAS = {
    "HUMAN": UNIPROT_DIR / "human_reviewed.fasta",
    "YEAST": UNIPROT_DIR / "yeast_reviewed.fasta",
    "ECOLI": UNIPROT_DIR / "ecoli_reviewed.fasta",
}


def _parse_accessions_from_fasta(fasta_path: Path) -> set[str]:
    """Return the set of UniProt primary accessions from a SwissProt FASTA.
    Each `>sp|ACC|NAME_SPECIES ...` header contributes the ACC token."""
    out: set[str] = set()
    with open(fasta_path, encoding="utf-8") as fh:
        for line in fh:
            if not line.startswith(">"):
                continue
            # Format: >sp|ACC|NAME_SPECIES ...
            parts = line[1:].split("|", 2)
            if len(parts) >= 2:
                out.add(parts[1].strip())
    return out


@lru_cache(maxsize=1)
def build_species_map() -> dict[str, str]:
    """Build a UniProt accession → species string map from the three
    cached SwissProt FASTAs. Returns a dict keyed by accession. Raises
    FileNotFoundError if a FASTA is missing — those files must be
    fetched once via the README instructions."""
    out: dict[str, str] = {}
    for species, path in _SPECIES_FASTAS.items():
        if not path.exists():
            raise FileNotFoundError(
                f"Missing SwissProt cache: {path}. Fetch via\n"
                f"  curl -s -o {path} "
                f"'https://rest.uniprot.org/uniprotkb/stream?query=reviewed:true"
                f"+AND+organism_id:<id>&format=fasta&compressed=false'"
            )
        for acc in _parse_accessions_from_fasta(path):
            out[acc] = species
    return out


def annotate_species_suffix(protein_ids: str, species_map: dict[str, str]) -> str:
    """Re-annotate a semicolon-separated `Protein.Ids` cell from
    quantmsdiann's bare-accession output (`Q96P70;P09417`) to the
    species-suffixed form ProteoBench expects (`Q96P70_HUMAN;P09417_HUMAN`).

    Accessions not in the species map are passed through unchanged
    (e.g. contaminants, decoys). The exact suffix style matches
    ProteoBench's species_mapper rules — a simple `_HUMAN` / `_YEAST` /
    `_ECOLI` token anywhere in the cell makes `str.contains("_HUMAN")`
    return True."""
    if not protein_ids or protein_ids == "nan":
        return protein_ids
    out_tokens: list[str] = []
    for acc in str(protein_ids).split(";"):
        acc_clean = acc.strip()
        species = species_map.get(acc_clean)
        if species:
            out_tokens.append(f"{acc_clean}_{species}")
        else:
            out_tokens.append(acc_clean)
    return ";".join(out_tokens)


# ---------------------------------------------------------------------------
# Wide-to-long melt of pr_matrix.tsv
# ---------------------------------------------------------------------------

def melt_pr_matrix(matrix_path: Path) -> pd.DataFrame:
    """Read a DIA-NN `diann_report.pr_matrix.tsv` and melt it to the
    long-format DataFrame ProteoBench's `convert_to_standard_format`
    consumes.

    Returns a DataFrame with columns
    `[Modified.Sequence, Protein.Ids, Precursor.Charge, Run,
    Precursor.Normalised]`. NaN intensities are dropped (a precursor
    that wasn't quantified in a given run contributes no row). The
    `Run` column has the `.raw` / `.d` / `.mzML` extension stripped to
    match ProteoBench's bare-filename condition_mapper keys.
    """
    wide = pd.read_csv(matrix_path, sep="\t")
    run_cols = [c for c in wide.columns if c not in PR_MATRIX_META_COLS]
    long_df = wide.melt(
        id_vars=["Modified.Sequence", "Protein.Ids", "Precursor.Charge"],
        value_vars=run_cols,
        var_name="Run",
        value_name="Precursor.Normalised",
    ).dropna(subset=["Precursor.Normalised"])
    long_df["Run"] = long_df["Run"].str.replace(_RUN_EXT_RE, "", regex=True)
    return long_df.reset_index(drop=True)


# ---------------------------------------------------------------------------
# Metric computation per analysis
# ---------------------------------------------------------------------------

def _proteobench_parse_settings_dir(module_id: str) -> Path:
    """Return the absolute path to the per-module parse-settings dir
    bundled with the installed `proteobench` package."""
    import proteobench
    return (
        Path(proteobench.__file__).parent / "io" / "parsing"
        / "io_parse_settings" / Path(*MODULE_TO_PARSE_DIR[module_id])
    )


def compute_proteobench_metrics(
    matrix_path: Path,
    module_id: str,
    *,
    software_version: str = "",
    proteobench_version_pin: str | None = None,
) -> dict:
    """Run ProteoBench's metric computation on one quantmsdiann
    `pr_matrix.tsv`. Returns the per-replicate-threshold `results`
    dict (keys `"1"..."6"`) plus headline aliases used by the figures.

    The function is intentionally side-effect-free — it does NOT write
    a cache file; the caller (`cached_proteobench_metrics`) handles
    caching so unit tests can exercise this path in isolation."""
    from proteobench.io.parsing.parse_settings import ParseSettingsBuilder
    from proteobench.score.quantscoresHYE import QuantScoresHYE
    from proteobench.datapoint.quant_datapoint import QuantDatapointHYE

    species_map = build_species_map()
    long_df = melt_pr_matrix(matrix_path)
    long_df["Protein.Ids"] = long_df["Protein.Ids"].astype(str).apply(
        lambda s: annotate_species_suffix(s, species_map)
    )

    parse_settings_dir = _proteobench_parse_settings_dir(module_id)
    parser = ParseSettingsBuilder(
        parse_settings_dir=str(parse_settings_dir),
        module_id=module_id,
    ).build_parser("DIA-NN")
    standard_format, replicate_to_raw = parser.convert_to_standard_format(long_df)

    score = QuantScoresHYE(
        "precursor ion",
        parser.species_expected_ratio(),
        parser.species_dict(),
    )
    intermediate = score.generate_intermediate(standard_format, replicate_to_raw)

    user_input = {
        "software_version": software_version,
        "search_engine": "DIA-NN",
        "search_engine_version": software_version,
        "ident_fdr_psm": 0.01,
        "ident_fdr_peptide": None,
        "ident_fdr_protein": 0.01,
        "enable_match_between_runs": False,
        "precursor_mass_tolerance": "",
        "fragment_mass_tolerance": "",
        "enzyme": "Trypsin/P",
        "allowed_miscleavages": 1,
        "min_peptide_length": 7,
        "max_peptide_length": 50,
    }
    datapoint = QuantDatapointHYE.generate_datapoint(
        intermediate, "DIA-NN", user_input,
    )

    # `datapoint.results` is keyed by int. JSON-serialise the keys as
    # strings for parity with ProteoBench's published submission JSONs.
    results_str_keys = {
        str(k): {kk: _jsonable(vv) for kk, vv in v.items()}
        for k, v in datapoint.results.items()
    }
    return {
        "matrix_path": str(matrix_path),
        "module_id": module_id,
        "software_version": software_version,
        "proteobench_version": (
            proteobench_version_pin
            or getattr(__import__("proteobench"), "__version__", "")
        ),
        "n_long_rows": int(len(long_df)),
        "n_standard_rows": int(len(standard_format)),
        "n_intermediate_rows": int(len(intermediate)),
        "replicate_to_raw": {
            k: list(v) for k, v in replicate_to_raw.items()
        },
        "results": results_str_keys,
    }


def _jsonable(value):
    """Coerce a metric value to a JSON-serialisable scalar. NumPy
    scalars get cast to Python types; lists pass through; everything
    else is stringified as a fallback."""
    if value is None:
        return None
    if hasattr(value, "item"):
        try:
            return value.item()
        except (TypeError, ValueError):
            pass
    if isinstance(value, (int, float, str, bool)):
        return value
    if isinstance(value, (list, tuple)):
        return [_jsonable(v) for v in value]
    return str(value)


# ---------------------------------------------------------------------------
# Cache layer
# ---------------------------------------------------------------------------

def metrics_cache_path(dataset: str, version: str) -> Path:
    return METRICS_CACHE_DIR / f"{dataset}_{version}.json"


def cached_proteobench_metrics(
    dataset: str,
    version: str,
    *,
    fetch: bool = True,
) -> dict:
    """Cache-or-compute wrapper. Returns the metrics dict either from
    the on-disk JSON cache or by computing it fresh.

    `fetch=False` skips computation entirely — the cache file must
    already exist (raises FileNotFoundError otherwise). Used by the
    figure rebuild path so it never silently invokes ProteoBench's
    parser on a CI runner that lacks the SwissProt FASTAs."""
    cache_path = metrics_cache_path(dataset, version)
    if cache_path.exists():
        with open(cache_path, encoding="utf-8") as fh:
            return json.load(fh)
    if not fetch:
        raise FileNotFoundError(
            f"Cache miss for {dataset}/{version} at {cache_path} and "
            f"fetch=False. Run with fetch=True (default) to populate."
        )
    module_id = DATASET_TO_MODULE[dataset]
    matrix = (
        DATA_DIR / dataset / version / "diann_report.pr_matrix.tsv"
    )
    payload = compute_proteobench_metrics(
        matrix, module_id, software_version=version,
    )
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    with open(cache_path, "w", encoding="utf-8") as fh:
        json.dump(payload, fh, indent=2)
    return payload


# ---------------------------------------------------------------------------
# Long-format extractor for the figure layer
# ---------------------------------------------------------------------------

def iter_metric_rows(
    payloads: Iterable[tuple[str, str, dict]],
) -> list[dict]:
    """Yield long-format rows from a sequence of (dataset, version,
    payload) tuples. Each row is one (dataset, version, threshold,
    metric) combination. Headline metrics + per-species variants are
    emitted; the figure consumer selects what it needs."""
    rows: list[dict] = []
    headline_metrics = [
        "nr_prec",
        "nr_prec_HUMAN", "nr_prec_YEAST", "nr_prec_ECOLI",
        "median_abs_epsilon_global", "mean_abs_epsilon_global",
        "median_abs_epsilon_eq_species", "mean_abs_epsilon_eq_species",
        "median_abs_epsilon_HUMAN", "median_abs_epsilon_YEAST", "median_abs_epsilon_ECOLI",
        "mean_log2_empirical_HUMAN", "mean_log2_empirical_YEAST",
        "mean_log2_empirical_ECOLI",
        "median_log2_empirical_HUMAN", "median_log2_empirical_YEAST",
        "median_log2_empirical_ECOLI",
        "CV_median", "CV_q75", "CV_q90", "CV_q95",
        "roc_auc", "roc_auc_directional",
        "variance_epsilon_global",
    ]
    for dataset, version, payload in payloads:
        results = payload.get("results", {})
        for threshold_str, metric_dict in results.items():
            try:
                threshold = int(threshold_str)
            except ValueError:
                continue
            for metric in headline_metrics:
                value = metric_dict.get(metric)
                if value is None:
                    continue
                rows.append({
                    "dataset": dataset,
                    "version": version,
                    "threshold": threshold,
                    "metric": metric,
                    "value": value,
                })
    return rows
