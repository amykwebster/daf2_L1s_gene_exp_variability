#!/usr/bin/env python3
"""
Aging analysis for muscle and neuron gene sets from AnnData object.

"""

from __future__ import annotations

import argparse
import json
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile
import textwrap
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

OUTPUT_ROOT_DEFAULT = Path("output")
OUTPUT_ROOT_DEFAULT.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str((OUTPUT_ROOT_DEFAULT / ".mplconfig").resolve()))

import anndata as ad
import matplotlib
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats
from statsmodels.stats.multitest import multipletests
from sklearn.decomposition import PCA

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def log(message: str) -> None:
    print(f"[INFO] {message}", flush=True)


def warn(message: str) -> None:
    print(f"[WARN] {message}", flush=True)


def ensure_dir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def write_text(path: Path, text: str) -> None:
    path.write_text(text, encoding="utf-8")


def age_sort_key(value: object) -> Tuple[int, float, str]:
    if pd.isna(value):
        return (2, float("inf"), "nan")
    text = str(value).strip()
    numbers = re.findall(r"[-+]?\d*\.?\d+", text)
    if numbers:
        return (0, float(numbers[0]), text)
    return (1, float("inf"), text)


def unique_preview(series: pd.Series, limit: int = 20) -> List[str]:
    values = [str(x) for x in series.dropna().astype(str).unique().tolist()]
    values = sorted(values, key=age_sort_key)
    if len(values) > limit:
        return values[:limit] + [f"... ({len(values) - limit} more)"]
    return values


def score_candidate_columns(columns: Sequence[str], keywords: Sequence[str]) -> List[str]:
    scored = []
    for column in columns:
        name = column.lower()
        score = sum(1 for keyword in keywords if keyword in name)
        if score > 0:
            scored.append((score, column))
    scored.sort(key=lambda item: (-item[0], item[1]))
    return [column for _, column in scored]


def infer_sample_from_obs_names(obs_names: Sequence[str]) -> pd.Series:
    inferred: List[str] = []
    for name in obs_names:
        text = str(name)
        if text.count("-") >= 2:
            inferred.append(text.split("-", 2)[-1])
        elif "_" in text:
            inferred.append(text.rsplit("_", 1)[0])
        else:
            inferred.append(text)
    return pd.Series(inferred, index=obs_names, name="inferred_sample_id", dtype="string")


def detect_metadata_fields(obs: pd.DataFrame) -> Tuple[Dict[str, List[str]], Dict[str, str], str]:
    candidates = {
        "age": score_candidate_columns(
            obs.columns,
            ["age", "timepoint", "time_point", "time", "day", "stage", "tp"],
        ),
        "annotation": score_candidate_columns(
            obs.columns,
            ["annotate", "annotation", "cell_type", "celltype", "cluster", "label", "identity"],
        ),
        "sample": score_candidate_columns(
            obs.columns,
            ["sample", "rep", "replicate", "batch", "library", "donor", "animal", "plate"],
        ),
    }

    choices: Dict[str, str] = {}
    if candidates["age"]:
        choices["age"] = candidates["age"][0]
    elif "timepoint" in obs.columns:
        choices["age"] = "timepoint"

    if candidates["annotation"]:
        choices["annotation"] = candidates["annotation"][0]
    elif "annotate_name" in obs.columns:
        choices["annotation"] = "annotate_name"

    inferred_sample = infer_sample_from_obs_names(obs.index)
    sample_choice = None
    if candidates["sample"]:
        sample_choice = candidates["sample"][0]
    elif inferred_sample.nunique(dropna=False) > 1:
        sample_choice = inferred_sample.name
    if sample_choice is not None:
        choices["sample"] = sample_choice

    return candidates, choices, inferred_sample.name


def summarize_dataset(
    adata: ad.AnnData,
    obs: pd.DataFrame,
    var: pd.DataFrame,
    output_dir: Path,
) -> Tuple[Dict[str, List[str]], Dict[str, str], pd.Series]:
    log("Inspecting dataset structure and metadata candidates.")
    candidates, choices, inferred_sample_name = detect_metadata_fields(obs)
    inferred_sample = infer_sample_from_obs_names(obs.index)

    summary_lines = [
        "Dataset summary",
        "===============",
        f"Cells: {adata.n_obs}",
        f"Genes: {adata.n_vars}",
        "",
        "obs columns:",
        ", ".join(map(str, obs.columns.tolist())) if len(obs.columns) else "(none)",
        "",
        "var columns:",
        ", ".join(map(str, var.columns.tolist())) if len(var.columns) else "(none)",
        "",
        "First rows of obs:",
        obs.head().to_string(),
        "",
        "Candidate metadata fields:",
    ]

    for key, columns in candidates.items():
        summary_lines.append(f"- {key}: {columns if columns else ['(none detected)']}")
        for column in columns[:3]:
            summary_lines.append(f"  * {column}: {unique_preview(obs[column])}")
    summary_lines.append(f"- sample (inferred from obs_names): {unique_preview(inferred_sample)}")
    summary_lines.append("")
    summary_lines.append("Chosen metadata fields:")
    for key in ["age", "annotation", "sample"]:
        value = choices.get(key, "(not found)")
        summary_lines.append(f"- {key}: {value}")

    write_text(output_dir / "dataset_structure_summary.txt", "\n".join(summary_lines) + "\n")

    metadata_candidates = {
        key: [
            {
                "column": column,
                "preview_unique_values": unique_preview(obs[column]),
            }
            for column in columns[:5]
        ]
        for key, columns in candidates.items()
    }
    metadata_candidates["sample_inferred_from_obs_names"] = [
        {
            "column": inferred_sample_name,
            "preview_unique_values": unique_preview(inferred_sample),
        }
    ]
    (output_dir / "metadata_candidates.json").write_text(
        json.dumps(
            {"candidates": metadata_candidates, "chosen_fields": choices},
            indent=2,
        ),
        encoding="utf-8",
    )

    return candidates, choices, inferred_sample


def choose_gene_names(var: pd.DataFrame, n_vars: int) -> pd.Index:
    if "gene_names" in var.columns:
        gene_names = var["gene_names"].astype(str)
        if gene_names.nunique() == n_vars:
            return pd.Index(gene_names.tolist(), name="gene")
    return pd.Index(var.index.astype(str).tolist(), name="gene")


def identify_tissue_labels(
    obs: pd.DataFrame,
    annotation_col: str,
    tissue_name: str,
) -> Tuple[pd.Series, pd.DataFrame]:
    annotations = obs[annotation_col].astype(str)
    if tissue_name == "muscle":
        pattern = r"muscle"
    elif tissue_name == "neuron":
        pattern = r"neuron|neuronal"
    else:
        raise ValueError(f"Unsupported tissue: {tissue_name}")

    mask = annotations.str.contains(pattern, case=False, regex=True, na=False)
    label_counts = (
        obs.loc[mask, annotation_col]
        .astype(str)
        .value_counts()
        .rename_axis("cell_type")
        .reset_index(name="n_cells")
    )
    label_counts["tissue_group"] = tissue_name
    return mask, label_counts


def sorted_age_values(series: pd.Series) -> List[str]:
    values = [str(x) for x in series.dropna().astype(str).unique().tolist()]
    return sorted(values, key=age_sort_key)


def build_pseudobulk_counts(
    adata: ad.AnnData,
    tissue_obs: pd.DataFrame,
    cell_indices: np.ndarray,
    sample_col: str,
    age_col: str,
    gene_names: pd.Index,
    min_cells_per_sample: int,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    log(f"Building pseudobulk counts for {len(tissue_obs)} cells across {tissue_obs[sample_col].nunique()} samples.")
    sample_table = (
        tissue_obs.groupby(sample_col, dropna=False)
        .agg(
            age=(age_col, lambda x: x.dropna().astype(str).iloc[0] if x.dropna().shape[0] else "NA"),
            n_cells=("__cell_id__", "size"),
        )
        .reset_index()
        .rename(columns={sample_col: "sample_id"})
    )
    sample_table["include_in_pseudobulk"] = sample_table["n_cells"] >= min_cells_per_sample

    included_samples = sample_table.loc[sample_table["include_in_pseudobulk"], "sample_id"].astype(str).tolist()
    if not included_samples:
        warn("No pseudobulk samples passed the minimum cell threshold.")
        empty = pd.DataFrame(columns=gene_names, index=pd.Index([], name="sample_id"))
        return empty, sample_table

    counts_by_sample: Dict[str, np.ndarray] = {}
    sample_lookup = tissue_obs[sample_col].astype(str).to_numpy()

    for sample_id in included_samples:
        local_idx = np.flatnonzero(sample_lookup == sample_id)
        global_idx = cell_indices[local_idx]
        matrix = adata.X[global_idx, :]
        summed = np.asarray(matrix.sum(axis=0)).ravel()
        counts_by_sample[sample_id] = summed.astype(np.int64, copy=False)

    counts_df = pd.DataFrame.from_dict(counts_by_sample, orient="index")
    counts_df.index.name = "sample_id"
    counts_df.columns = gene_names
    counts_df = counts_df.loc[included_samples]
    return counts_df, sample_table


def choose_comparison_ages(
    sample_table: pd.DataFrame,
    all_age_values: Sequence[str],
    min_replicates_per_group: int,
) -> Dict[str, Optional[str]]:
    result = {
        "youngest_overall": all_age_values[0] if all_age_values else None,
        "oldest_overall": all_age_values[-1] if all_age_values else None,
        "young_for_de": None,
        "old_for_de": None,
        "comparison_note": "",
    }

    usable = sample_table.loc[sample_table["include_in_pseudobulk"]].copy()
    if usable.empty:
        result["comparison_note"] = "No pseudobulk samples passed filtering."
        return result

    usable["age"] = usable["age"].astype(str)
    age_to_n = usable.groupby("age")["sample_id"].nunique().to_dict()

    for age in all_age_values:
        if age_to_n.get(age, 0) >= min_replicates_per_group:
            result["young_for_de"] = age
            break

    for age in reversed(all_age_values):
        if age_to_n.get(age, 0) >= min_replicates_per_group:
            result["old_for_de"] = age
            break

    if result["young_for_de"] is None and all_age_values:
        result["young_for_de"] = all_age_values[0]
    if result["old_for_de"] is None and all_age_values:
        result["old_for_de"] = all_age_values[-1]

    if (
        result["young_for_de"] == result["youngest_overall"]
        and result["old_for_de"] == result["oldest_overall"]
    ):
        result["comparison_note"] = "Using the youngest and oldest overall ages for DE."
    else:
        result["comparison_note"] = (
            "Oldest overall age was not statistically usable after pseudobulk filtering "
            f"(or lacked replicates). Using {result['young_for_de']} vs {result['old_for_de']} for DE."
        )
    return result


def run_edger_de(
    counts: pd.DataFrame,
    sample_meta: pd.DataFrame,
    out_csv: Path,
) -> Optional[pd.DataFrame]:
    if shutil.which("Rscript") is None:
        return None

    check = subprocess.run(
        ["Rscript", "-e", "quit(status=ifelse(requireNamespace('edgeR', quietly=TRUE), 0, 1))"],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    if check.returncode != 0:
        return None

    with tempfile.TemporaryDirectory(prefix="edger_", dir=str(out_csv.parent)) as tmp_dir:
        tmp_path = Path(tmp_dir)
        counts_path = tmp_path / "counts.csv"
        meta_path = tmp_path / "sample_metadata.csv"
        script_path = tmp_path / "run_edger.R"
        counts.T.to_csv(counts_path)
        sample_meta.to_csv(meta_path, index=False)

        r_script = textwrap.dedent(
            """
            suppressPackageStartupMessages({
              library(edgeR)
            })

            args <- commandArgs(trailingOnly = TRUE)
            counts_path <- args[1]
            meta_path <- args[2]
            out_path <- args[3]

            counts <- read.csv(counts_path, row.names = 1, check.names = FALSE)
            meta <- read.csv(meta_path, stringsAsFactors = FALSE, check.names = FALSE)
            rownames(meta) <- meta$sample_id
            meta <- meta[colnames(counts), , drop = FALSE]
            group <- factor(meta$group, levels = c("young", "old"))

            y <- DGEList(counts = counts)
            keep <- filterByExpr(y, group = group)
            y <- y[keep, , keep.lib.sizes = FALSE]
            y <- calcNormFactors(y)
            design <- model.matrix(~ group)
            y <- estimateDisp(y, design, robust = TRUE)
            fit <- glmQLFit(y, design, robust = TRUE)
            qlf <- glmQLFTest(fit, coef = "groupold")
            tt <- topTags(qlf, n = Inf, sort.by = "none")$table
            tt$gene <- rownames(tt)
            tt$average_expression <- rowMeans(cpm(y, log = TRUE, prior.count = 2))
            tt <- tt[, c("gene", "logFC", "PValue", "FDR", "average_expression")]
            write.csv(tt, out_path, row.names = FALSE)
            """
        ).strip()
        write_text(script_path, r_script + "\n")

        run = subprocess.run(
            ["Rscript", str(script_path), str(counts_path), str(meta_path), str(out_csv)],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            env={**os.environ, "LC_ALL": "C"},
        )
        if run.returncode != 0:
            warn("edgeR failed; falling back to Python DE.")
            warn(run.stderr.strip() or "edgeR returned a non-zero exit code.")
            return None

    results = pd.read_csv(out_csv)
    results = results.rename(
        columns={
            "logFC": "log2FC",
            "PValue": "pvalue",
            "FDR": "fdr",
        }
    )
    return results


def run_fallback_de(counts: pd.DataFrame, sample_meta: pd.DataFrame) -> pd.DataFrame:
    log("Running fallback DE on pseudobulk logCPM profiles.")
    counts = counts.astype(float)
    lib_sizes = counts.sum(axis=1)
    log_cpm = np.log2(((counts + 0.5).div(lib_sizes + 1.0, axis=0) * 1e6) + 1.0)

    group = sample_meta.set_index("sample_id").loc[counts.index, "group"]
    young = log_cpm.loc[group == "young"]
    old = log_cpm.loc[group == "old"]

    average_expression = log_cpm.mean(axis=0)
    log2fc = old.mean(axis=0) - young.mean(axis=0)

    if min(len(young), len(old)) >= 2:
        stat, pvalue = stats.ttest_ind(
            old.to_numpy(),
            young.to_numpy(),
            axis=0,
            equal_var=False,
            nan_policy="omit",
        )
        pvalue = np.asarray(pvalue, dtype=float)
        pvalue[~np.isfinite(pvalue)] = 1.0
        fdr = multipletests(pvalue, method="fdr_bh")[1]
    else:
        warn("Insufficient replicate counts for a statistical test; saving effect-size-only ranking.")
        stat = np.full(log2fc.shape[0], np.nan)
        pvalue = np.full(log2fc.shape[0], np.nan)
        fdr = np.full(log2fc.shape[0], np.nan)

    results = pd.DataFrame(
        {
            "gene": counts.columns,
            "log2FC": log2fc.to_numpy(),
            "pvalue": pvalue,
            "fdr": fdr,
            "average_expression": average_expression.to_numpy(),
            "statistic": stat,
        }
    )
    return results


def run_de_analysis(
    counts: pd.DataFrame,
    sample_table: pd.DataFrame,
    young_age: str,
    old_age: str,
    tissue_dir: Path,
) -> Tuple[pd.DataFrame, str]:
    usable = sample_table.loc[sample_table["include_in_pseudobulk"]].copy()
    usable["age"] = usable["age"].astype(str)
    de_meta = usable.loc[usable["age"].isin([young_age, old_age]), ["sample_id", "age", "n_cells"]].copy()
    de_meta["group"] = np.where(de_meta["age"] == young_age, "young", "old")

    counts_de = counts.loc[de_meta["sample_id"]].copy()
    gene_filter = (counts_de.sum(axis=0) >= 10) & ((counts_de > 0).sum(axis=0) >= 2)
    counts_de = counts_de.loc[:, gene_filter]

    if counts_de.empty:
        warn("Gene filtering removed all genes. Returning an empty DE table.")
        empty = pd.DataFrame(columns=["gene", "log2FC", "pvalue", "fdr", "average_expression"])
        return empty, "none"

    de_meta.to_csv(tissue_dir / "pseudobulk_samples_used_for_de.csv", index=False)
    counts_de.to_csv(tissue_dir / "pseudobulk_counts_used_for_de.csv.gz", compression="gzip")

    results_path = tissue_dir / "de_results_full.csv"
    edger_results = None
    method = "fallback_logcpm_ttest"
    if len(de_meta["group"].unique()) == 2:
        edger_results = run_edger_de(counts_de, de_meta, results_path)
    if edger_results is not None:
        results = edger_results
        method = "edgeR_glmQLF"
    else:
        results = run_fallback_de(counts_de, de_meta)
        results.to_csv(results_path, index=False)

    if "fdr" not in results.columns:
        results["fdr"] = np.nan
    if "pvalue" not in results.columns:
        results["pvalue"] = np.nan
    if "average_expression" not in results.columns:
        results["average_expression"] = np.nan

    # Preserve all filtered genes and provide a signed ranking metric.
    pvalue_for_rank = results["pvalue"].fillna(1.0).clip(lower=1e-300)
    if np.allclose(pvalue_for_rank.to_numpy(), 1.0, equal_nan=False):
        rank_score = results["log2FC"].fillna(0.0)
    else:
        rank_score = results["log2FC"].fillna(0.0) * -np.log10(pvalue_for_rank)
    results["rank_score"] = rank_score
    results = results.sort_values(["fdr", "pvalue", "log2FC"], ascending=[True, True, False], na_position="last")
    results.to_csv(results_path, index=False)
    return results, method


def save_gene_sets(
    results: pd.DataFrame,
    tissue_name: str,
    tissue_dir: Path,
    gene_set_summary_rows: List[Dict[str, object]],
    fdr_threshold: float,
    min_abs_log2fc: float,
) -> None:
    gene_set_dir = ensure_dir(tissue_dir / "gene_sets")
    ranked_dir = ensure_dir(tissue_dir / "ranked_lists")

    valid_fdr = results["fdr"].notna() if "fdr" in results.columns else pd.Series(False, index=results.index)
    up_mask = valid_fdr & (results["fdr"] < fdr_threshold) & (results["log2FC"] > min_abs_log2fc)
    down_mask = valid_fdr & (results["fdr"] < fdr_threshold) & (results["log2FC"] < -min_abs_log2fc)

    direction_map = {
        "Up": results.loc[up_mask].sort_values("log2FC", ascending=False),
        "Down": results.loc[down_mask].sort_values("log2FC", ascending=True),
    }

    for direction, subset in direction_map.items():
        set_name = f"{tissue_name.capitalize()}_Age_{direction}"
        subset.to_csv(gene_set_dir / f"{set_name}.csv", index=False)
        write_text(gene_set_dir / f"{set_name}.txt", "\n".join(subset["gene"].astype(str).tolist()) + ("\n" if len(subset) else ""))
        gene_set_summary_rows.append(
            {
                "gene_set": set_name,
                "n_genes": int(len(subset)),
                "fdr_threshold": fdr_threshold,
                "min_abs_log2fc": min_abs_log2fc,
            }
        )

    ranked = results.sort_values("rank_score", ascending=False).copy()
    ranked.to_csv(ranked_dir / f"{tissue_name}_aging_ranked_genes.csv", index=False)
    ranked[["gene", "rank_score", "log2FC", "pvalue", "fdr"]].to_csv(
        ranked_dir / f"{tissue_name}_aging_ranked_genes_for_gsea.tsv",
        sep="\t",
        index=False,
    )


def plot_bar_counts(
    counts_series: pd.Series,
    title: str,
    ylabel: str,
    out_path: Path,
) -> None:
    if counts_series.empty:
        return
    ordered = counts_series.sort_index(key=lambda idx: [age_sort_key(x) for x in idx])
    plt.figure(figsize=(6, 4))
    sns.barplot(x=ordered.index.astype(str), y=ordered.values, color="#4C78A8")
    plt.title(title)
    plt.xlabel("Age")
    plt.ylabel(ylabel)
    plt.tight_layout()
    plt.savefig(out_path, dpi=200)
    plt.close()


def plot_volcano(results: pd.DataFrame, title: str, out_path: Path) -> None:
    if results.empty or "pvalue" not in results.columns:
        return
    plot_df = results.copy()
    plot_df["neg_log10_pvalue"] = -np.log10(plot_df["pvalue"].fillna(1.0).clip(lower=1e-300))
    plot_df["significant"] = plot_df["fdr"].fillna(1.0) < 0.05

    plt.figure(figsize=(6, 5))
    sns.scatterplot(
        data=plot_df,
        x="log2FC",
        y="neg_log10_pvalue",
        hue="significant",
        palette={True: "#D62728", False: "#7F7F7F"},
        s=18,
        linewidth=0,
        legend=False,
    )
    plt.axvline(0.0, color="black", linestyle="--", linewidth=0.8)
    plt.axhline(-math.log10(0.05), color="black", linestyle=":", linewidth=0.8)
    plt.title(title)
    plt.xlabel("log2 fold change (old vs young)")
    plt.ylabel("-log10 p-value")
    plt.tight_layout()
    plt.savefig(out_path, dpi=200)
    plt.close()


def plot_heatmap(
    counts: pd.DataFrame,
    sample_table: pd.DataFrame,
    results: pd.DataFrame,
    title: str,
    out_path: Path,
) -> None:
    if counts.empty or results.empty:
        return
    sample_meta = sample_table.loc[sample_table["include_in_pseudobulk"], ["sample_id", "age"]].copy()
    sample_meta = sample_meta.set_index("sample_id").loc[counts.index]

    up_genes = (
        results.loc[(results["fdr"].fillna(1.0) < 0.05) & (results["log2FC"] > 0)]
        .sort_values(["fdr", "log2FC"], ascending=[True, False])
        .head(50)["gene"]
        .astype(str)
        .tolist()
    )
    down_genes = (
        results.loc[(results["fdr"].fillna(1.0) < 0.05) & (results["log2FC"] < 0)]
        .sort_values(["fdr", "log2FC"], ascending=[True, True])
        .head(50)["gene"]
        .astype(str)
        .tolist()
    )
    selected_genes = up_genes + [gene for gene in down_genes if gene not in up_genes]
    if not selected_genes:
        warn(f"Skipping heatmap for {title}: no significant genes passed thresholds.")
        return

    heat = counts.loc[:, [gene for gene in selected_genes if gene in counts.columns]].astype(float)
    lib_sizes = heat.sum(axis=1)
    log_cpm = np.log2(((heat + 0.5).div(lib_sizes + 1.0, axis=0) * 1e6) + 1.0)
    z = log_cpm.sub(log_cpm.mean(axis=0), axis=1).div(log_cpm.std(axis=0).replace(0, np.nan), axis=1)
    z = z.fillna(0.0).T

    ordered_samples = sample_meta.sort_values("age", key=lambda s: s.map(age_sort_key)).index.tolist()
    z = z.loc[:, ordered_samples]

    plt.figure(figsize=(max(6, len(ordered_samples) * 0.7), max(8, z.shape[0] * 0.12)))
    sns.heatmap(z, cmap="vlag", center=0, xticklabels=True, yticklabels=True)
    plt.title(title)
    plt.xlabel("Pseudobulk sample")
    plt.ylabel("Gene")
    plt.tight_layout()
    plt.savefig(out_path, dpi=200)
    plt.close()


def plot_pca(
    counts: pd.DataFrame,
    sample_table: pd.DataFrame,
    title: str,
    out_path: Path,
) -> None:
    if counts.shape[0] < 2 or counts.shape[1] < 2:
        return
    lib_sizes = counts.sum(axis=1)
    log_cpm = np.log2(((counts + 0.5).div(lib_sizes + 1.0, axis=0) * 1e6) + 1.0)
    variable = log_cpm.var(axis=0).sort_values(ascending=False).head(min(1000, log_cpm.shape[1])).index
    matrix = log_cpm.loc[:, variable]
    pca = PCA(n_components=2)
    coords = pca.fit_transform(matrix)

    sample_meta = sample_table.loc[sample_table["include_in_pseudobulk"], ["sample_id", "age"]].copy()
    sample_meta = sample_meta.set_index("sample_id").loc[counts.index]
    plot_df = pd.DataFrame(coords, columns=["PC1", "PC2"], index=counts.index).join(sample_meta)

    plt.figure(figsize=(6, 5))
    sns.scatterplot(data=plot_df, x="PC1", y="PC2", hue="age", s=90)
    for sample_id, row in plot_df.iterrows():
        plt.text(row["PC1"], row["PC2"], str(sample_id), fontsize=8)
    plt.title(title)
    plt.tight_layout()
    plt.savefig(out_path, dpi=200)
    plt.close()


def plot_qc_for_tissue(
    tissue_name: str,
    tissue_obs: pd.DataFrame,
    sample_table: pd.DataFrame,
    counts: pd.DataFrame,
    results: pd.DataFrame,
    qc_dir: Path,
) -> None:
    plot_bar_counts(
        tissue_obs["age_value"].value_counts(),
        title=f"{tissue_name.capitalize()} cells by age",
        ylabel="Cells",
        out_path=qc_dir / f"{tissue_name}_cell_counts_by_age.png",
    )
    plot_bar_counts(
        sample_table.loc[sample_table["include_in_pseudobulk"], "age"].value_counts(),
        title=f"{tissue_name.capitalize()} pseudobulk samples by age",
        ylabel="Pseudobulk samples",
        out_path=qc_dir / f"{tissue_name}_pseudobulk_counts_by_age.png",
    )
    plot_volcano(
        results,
        title=f"{tissue_name.capitalize()} aging volcano",
        out_path=qc_dir / f"{tissue_name}_aging_volcano.png",
    )
    plot_heatmap(
        counts,
        sample_table,
        results,
        title=f"{tissue_name.capitalize()} top age-associated genes",
        out_path=qc_dir / f"{tissue_name}_aging_heatmap.png",
    )
    plot_pca(
        counts,
        sample_table,
        title=f"{tissue_name.capitalize()} pseudobulk PCA",
        out_path=qc_dir / f"{tissue_name}_pseudobulk_pca.png",
    )


def main() -> int:
    parser = argparse.ArgumentParser(description="First-pass muscle/neuron aging analysis from AnnData.")
    parser.add_argument("--input", default="ad_worm_aging.h5ad", help="Input AnnData file.")
    parser.add_argument("--output", default="output", help="Output directory.")
    parser.add_argument("--min-cells-per-sample", type=int, default=20, help="Minimum cells required per pseudobulk sample.")
    parser.add_argument("--min-replicates-per-group", type=int, default=2, help="Minimum pseudobulk replicates needed to use an age in DE.")
    parser.add_argument("--fdr-threshold", type=float, default=0.05, help="FDR cutoff for gene sets.")
    parser.add_argument("--min-abs-log2fc", type=float, default=0.0, help="Absolute log2FC cutoff for gene sets.")
    args = parser.parse_args()

    input_path = Path(args.input)
    output_root = ensure_dir(Path(args.output))
    metadata_dir = ensure_dir(output_root / "metadata_summary")
    muscle_dir = ensure_dir(output_root / "muscle")
    neuron_dir = ensure_dir(output_root / "neuron")
    qc_dir = ensure_dir(output_root / "qc")

    if not input_path.exists():
        raise FileNotFoundError(f"Input file not found: {input_path}")

    sns.set_theme(style="whitegrid")

    log(f"Opening {input_path} in backed mode.")
    adata = ad.read_h5ad(input_path, backed="r")
    obs = adata.obs.copy()
    var = adata.var.copy()

    _, metadata_choices, inferred_sample = summarize_dataset(adata, obs, var, metadata_dir)

    if "age" not in metadata_choices:
        raise RuntimeError("Could not detect an age/timepoint column in adata.obs.")
    if "annotation" not in metadata_choices:
        raise RuntimeError("Could not detect a cell annotation column in adata.obs.")

    age_col = metadata_choices["age"]
    annotation_col = metadata_choices["annotation"]
    sample_col = metadata_choices.get("sample", inferred_sample.name)

    if sample_col == inferred_sample.name:
        obs[sample_col] = inferred_sample.astype(str).values
    elif sample_col not in obs.columns:
        obs[sample_col] = inferred_sample.astype(str).values
        warn(f"Selected sample column `{sample_col}` was missing; using inferred sample IDs instead.")
        sample_col = inferred_sample.name
        metadata_choices["sample"] = sample_col

    obs["age_value"] = obs[age_col].astype(str)
    obs["annotation_value"] = obs[annotation_col].astype(str)
    obs["__cell_id__"] = obs.index.astype(str)

    all_age_values = sorted_age_values(obs["age_value"])
    age_counts = (
        obs["age_value"].value_counts()
        .rename_axis("age")
        .reset_index(name="n_cells")
        .sort_values("age", key=lambda s: s.map(age_sort_key))
    )
    age_counts.to_csv(metadata_dir / "all_age_values_and_cell_counts.csv", index=False)

    gene_names = choose_gene_names(var, adata.n_vars)

    gene_set_summary_rows: List[Dict[str, object]] = []
    comparison_summary: Dict[str, Dict[str, Optional[str]]] = {}
    de_method_summary: Dict[str, str] = {}

    for tissue_name, tissue_dir in [("muscle", muscle_dir), ("neuron", neuron_dir)]:
        log(f"Processing tissue group: {tissue_name}")
        mask, label_table = identify_tissue_labels(obs, annotation_col, tissue_name)
        label_table.to_csv(tissue_dir / f"{tissue_name}_cell_types.csv", index=False)

        tissue_obs = obs.loc[mask].copy()
        if tissue_obs.empty:
            warn(f"No cells matched the {tissue_name} definition.")
            comparison_summary[tissue_name] = {
                "youngest_overall": all_age_values[0] if all_age_values else None,
                "oldest_overall": all_age_values[-1] if all_age_values else None,
                "young_for_de": None,
                "old_for_de": None,
                "comparison_note": f"No {tissue_name} cells were found.",
            }
            de_method_summary[tissue_name] = "not_run"
            continue

        tissue_obs.to_csv(tissue_dir / f"{tissue_name}_per_cell_metadata.csv", index=True)
        tissue_age_counts = (
            tissue_obs["age_value"].value_counts()
            .rename_axis("age")
            .reset_index(name="n_cells")
            .sort_values("age", key=lambda s: s.map(age_sort_key))
        )
        tissue_age_counts.to_csv(tissue_dir / f"{tissue_name}_cell_counts_by_age.csv", index=False)

        cell_indices = np.flatnonzero(mask.to_numpy())
        counts, sample_table = build_pseudobulk_counts(
            adata=adata,
            tissue_obs=tissue_obs,
            cell_indices=cell_indices,
            sample_col=sample_col,
            age_col="age_value",
            gene_names=gene_names,
            min_cells_per_sample=args.min_cells_per_sample,
        )
        sample_table = sample_table.sort_values("age", key=lambda s: s.map(age_sort_key))
        sample_table.to_csv(tissue_dir / f"{tissue_name}_pseudobulk_samples.csv", index=False)

        comparison = choose_comparison_ages(sample_table, all_age_values, args.min_replicates_per_group)
        comparison_summary[tissue_name] = comparison
        write_text(tissue_dir / "comparison_summary.json", json.dumps(comparison, indent=2) + "\n")

        if counts.empty:
            warn(f"No usable pseudobulk counts for {tissue_name}.")
            empty_results = pd.DataFrame(columns=["gene", "log2FC", "pvalue", "fdr", "average_expression", "rank_score"])
            empty_results.to_csv(tissue_dir / "de_results_full.csv", index=False)
            de_method_summary[tissue_name] = "not_run"
            continue

        counts.to_csv(tissue_dir / f"{tissue_name}_pseudobulk_counts_all_genes.csv.gz", compression="gzip")

        young_age = comparison.get("young_for_de")
        old_age = comparison.get("old_for_de")
        if not young_age or not old_age or young_age == old_age:
            warn(f"Could not identify a valid age comparison for {tissue_name}.")
            empty_results = pd.DataFrame(columns=["gene", "log2FC", "pvalue", "fdr", "average_expression", "rank_score"])
            empty_results.to_csv(tissue_dir / "de_results_full.csv", index=False)
            de_method_summary[tissue_name] = "not_run"
            continue

        results, de_method = run_de_analysis(
            counts=counts,
            sample_table=sample_table,
            young_age=young_age,
            old_age=old_age,
            tissue_dir=tissue_dir,
        )
        de_method_summary[tissue_name] = de_method
        save_gene_sets(
            results=results,
            tissue_name=tissue_name,
            tissue_dir=tissue_dir,
            gene_set_summary_rows=gene_set_summary_rows,
            fdr_threshold=args.fdr_threshold,
            min_abs_log2fc=args.min_abs_log2fc,
        )
        plot_qc_for_tissue(
            tissue_name=tissue_name,
            tissue_obs=tissue_obs,
            sample_table=sample_table,
            counts=counts,
            results=results,
            qc_dir=qc_dir,
        )

    pd.DataFrame(gene_set_summary_rows).to_csv(output_root / "metadata_summary" / "gene_set_summary.csv", index=False)
    write_text(
        metadata_dir / "analysis_choices.json",
        json.dumps(
            {
                "metadata_choices": metadata_choices,
                "comparison_summary": comparison_summary,
                "de_method_summary": de_method_summary,
                "parameters": {
                    "min_cells_per_sample": args.min_cells_per_sample,
                    "min_replicates_per_group": args.min_replicates_per_group,
                    "fdr_threshold": args.fdr_threshold,
                    "min_abs_log2fc": args.min_abs_log2fc,
                },
            },
            indent=2,
        )
        + "\n",
    )

    save_readme(output_root, input_path, metadata_choices, comparison_summary, de_method_summary)
    log("Analysis complete.")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:  # pragma: no cover - entrypoint guard
        print(f"[ERROR] {exc}", file=sys.stderr)
        raise
