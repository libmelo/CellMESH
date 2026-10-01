"""Reproduce the methodological probes in ALGORITHM_REVIEW_PD_2026-09-22.md.

Run from the package checkout:
    python docs/review_pd_scoring_2026_09_22.py > /tmp/pd_review.json
This script uses synthetic expression and reads the packaged priors only.
It does not change scoring code or existing data.
"""
from pathlib import Path
import hashlib
import importlib.metadata
import json
import platform
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import numpy as np
import pandas as pd
from anndata import AnnData

from cellmesh import compute_metabolite_availability, load_cell_mesh_database, run_cell_mesh
from cellmesh._numerics import _positive_reference_scores, _reaction_activity
from cellmesh.score import _get_reaction_gene_sets, _normalize_enzyme_metabolite


def adata(values, genes, labels, samples=None):
    obs = pd.DataFrame({"cell_type": labels}, index=[f"c{i}" for i in range(len(labels))])
    if samples is not None:
        obs["sample"] = samples
    return AnnData(np.asarray(values, dtype=float), obs=obs, var=pd.DataFrame(index=genes))


def enzyme(rows):
    return pd.DataFrame([
        {"metabolite": "Met", "hmdb_id": "HMDB0000001", "gene": g, "role": role, "reaction": rxn}
        for g, role, rxn in rows
    ])


SENSOR = pd.DataFrame([{
    "metabolite": "Met", "hmdb_id": "HMDB0000001", "sensor_gene": "R",
    "sensor_type": "Cell surface receptor",
}])
PRIOR = enzyme([("P", "production", "p"), ("D", "degradation", "d")])


def event(result, sender="A", receiver="A"):
    return result.events.query("sender == @sender and receiver == @receiver").iloc[0]


def scalar(value):
    return float(value) if np.isfinite(value) else None


def main():
    out = {
        "date": "2026-09-22",
        "git_head": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
        "environment": {"python": platform.python_version(), **{
            p: importlib.metadata.version(p) for p in ("numpy", "pandas", "scipy", "anndata")
        }},
        "source_sha256": {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
                          for p in sorted((ROOT / "cellmesh").glob("*.py"))},
        "proposal_sha256": hashlib.sha256((ROOT.parent / "mCCC_PD_normalization_scoring_review.md").read_bytes()).hexdigest(),
    }

    # Identical production and sensor evidence; alter only consumption coverage.
    missing = []
    for state in ("supported", "prior_missing", "prior_gene_unavailable", "prior_no_expression"):
        genes = ["P", "D", "R"] if state != "prior_gene_unavailable" else ["P", "R"]
        row = [1, 0 if state == "prior_no_expression" else 1, 1] if "D" in genes else [1, 1]
        a = adata([row] * 4, genes, ["A", "A", "B", "B"])
        prior = PRIOR if state != "prior_missing" else PRIOR.query("role == 'production'")
        result = run_cell_mesh(a, prior, SENSOR, min_cells=1)
        ar = result.availability_results
        missing.append({"state": ar["metadata"].iloc[0]["consumption_status"],
                        "C_ref": scalar(ar["C_ref"].iloc[0]),
                        "sender_A": float(ar["availability"].iloc[0]["A"]),
                        "event_A_A": float(event(result)["cell_mesh_score"])})
    out["consumption_coverage"] = missing

    # A tenfold change in production alone is removed by within-sample references.
    a = adata([[1, 1, 1]] * 4 + [[10, 1, 1]] * 4, ["P", "D", "R"],
              ["A", "A", "B", "B"] * 2, ["S1"] * 4 + ["S2"] * 4)
    result = run_cell_mesh(a, PRIOR, SENSOR, sample_key="sample", sample_mode="sample_aware", min_cells=1)
    out["within_sample_rescaling"] = [{
        "sample": sample, "P_A": float(ar["P"].iloc[0]["A"]),
        "P_ref": float(ar["P_ref"].iloc[0]), "P_score_A": float(ar["P_score"].iloc[0]["A"]),
        "sender_A": float(ar["availability"].iloc[0]["A"]),
        "event_A_A": float(result.sample_events.query("sample == @sample and sender == 'A' and receiver == 'A'")["cell_mesh_score"].iloc[0]),
    } for sample, ar in result.availability_results["availability_by_sample"].items()]

    out["positive_reference_discontinuity"] = []
    for values in ([1.0, 0.0], [1.0, 1e-12]):
        scores, ref = _positive_reference_scores(np.array(values), "mean", "audit")
        out["positive_reference_discontinuity"].append({"values": values, "reference": ref, "scores": scores.tolist()})

    out["saturation_vs_relative_ratio"] = []
    for p, d in [(4, 2), (1, .5), (.2, .1), (100, 10), (2, 1)]:
        ps, ds = p / (1 + p), d / (1 + d)
        out["saturation_vs_relative_ratio"].append({"p": p, "d": d, "R_rel": p / (p + d),
            "saturated_ratio": ps / (ps + ds), "sender_base": ps * ps / (ps + ds)})

    # A measured zero subunit and an unmeasured subunit are numerically different.
    complex_prior = enzyme([("G1;G2", "production", "r")])
    out["partial_reaction_coverage"] = []
    for genes, values in [(["G1", "G2"], [[100, 0], [1, 1]]), (["G1"], [[100], [1]])]:
        ar = compute_metabolite_availability(adata(values, genes, ["A", "B"]), complex_prior,
                                            sender_abundance_exponent=0, min_cells=1)
        out["partial_reaction_coverage"].append({"measured_genes": genes, "P_A": float(ar["P"].iloc[0]["A"]),
            "P_B": float(ar["P"].iloc[0]["B"]), "sender_A": float(ar["availability"].iloc[0]["A"])})
    x = np.array([[100, 0], [0, 100]], dtype=float)
    out["aggregation_noncommutativity"] = {
        "activity_of_group_mean": float(_reaction_activity(x.mean(axis=0, keepdims=True))[0]),
        "mean_of_cell_activities": float(_reaction_activity(x).mean()),
        "cells_coexpressing_both_genes": int(np.all(x > 0, axis=1).sum()),
    }

    subset_prior = enzyme([("G1", "production", "r1"), ("G1;G2", "production", "r2")])
    out["strict_subset_reaction_selection"] = _get_reaction_gene_sets(_normalize_enzyme_metabolite(subset_prior)).to_dict("records")

    out["composition_only_change"] = []
    for n_a in (100, 40):
        result = run_cell_mesh(adata([[1, 1, 1]] * 200, ["P", "D", "R"],
                                    ["A"] * n_a + ["B"] * (200 - n_a)), PRIOR, SENSOR, min_cells=1)
        out["composition_only_change"].append({"fraction_A": n_a / 200, "event_A_A": float(event(result)["cell_mesh_score"]),
                                               "event_B_B": float(event(result, "B", "B")["cell_mesh_score"])})

    out["local_expression_support"] = []
    for positive_cells in (1, 100):
        x = np.zeros((200, 3))
        x[:, 1] = 1
        x[:positive_cells, 0] = x[:positive_cells, 2] = 1
        result = run_cell_mesh(adata(x, ["P", "D", "R"], ["A"] * 100 + ["B"] * 100), PRIOR, SENSOR)
        ev = event(result)
        out["local_expression_support"].append({"positive_cells_A": positive_cells,
            "group_size_A": 100, "raw_counts_each_signal_A": positive_cells,
            "sensor_expr_frac": float(ev["sensor_expr_frac"]), "passes_min_cells": bool(ev["passes_min_cells"]),
            "event_A_A": float(ev["cell_mesh_score"])})

    # Current packaged prior coverage, with no expression-based filtering.
    ep, sp = load_cell_mesh_database()
    all_ids, shared = set(ep.hmdb_id), set(ep.hmdb_id) & set(sp.hmdb_id)
    role_ids = {role: set(ep.loc[ep.role.eq(role), "hmdb_id"]) for role in ("production", "degradation", "export")}
    producing_shared = shared & role_ids["production"]
    parsed = _normalize_enzyme_metabolite(ep)
    grouped = parsed.groupby(["hmdb_id", "reaction", "direction"]).gene.apply(lambda x: frozenset(x))
    unique_sets = {}
    for (hmdb, rxn, direction), genes in grouped.items():
        unique_sets.setdefault((hmdb, direction), set()).add(genes)
    n_unique = sum(map(len, unique_sets.values()))
    n_subsets = sum(sum(any(g < h for h in group) for g in group) for group in unique_sets.values())
    selected = _get_reaction_gene_sets(parsed)
    out["default_prior_audit"] = {
        "input_provenance": {"enzyme": ep.attrs["input_provenance"], "sensor": sp.attrs["input_provenance"]},
        "enzyme_hmdb_count": len(all_ids), "shared_hmdb_count": len(shared),
        "shared_with_production": len(producing_shared),
        "shared_production_without_consumption": len(producing_shared - role_ids["degradation"]),
        "shared_with_export": len(shared & role_ids["export"]),
        "reaction_groups_before_set_selection": len(grouped),
        "unique_direction_gene_sets": n_unique, "strict_subsets_omitted": n_subsets,
        "retained_reaction_sets": len(selected),
        "retained_multigene_sets": int(selected.genes.map(len).gt(1).sum()),
    }
    print(json.dumps(out, ensure_ascii=False, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
