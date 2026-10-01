"""Descriptive reports built from captured scoring values, not a second scorer.

Coverage describes the input gene panel, not biological function or confidence.
The reports do not change the tested event universe, scores, or p-values.
"""
from __future__ import annotations

import html
import json
import os
from pathlib import Path
import tempfile

import numpy as np
import pandas as pd

from .positivity import INHERIT

from .database import _normalize_hmdb_id, _valid_hmdb_mask
from .score import _normalize_enzyme_metabolite, _get_reaction_gene_sets


REPORT_SCHEMA_VERSION = 1
EVENT_KEYS = ["sender", "receiver", "hmdb_id", "sensor_gene"]


def _units(result):
    availability = result.availability_results
    if not isinstance(availability, dict):
        raise ValueError("Explanation requires availability_results retained by run_cell_mesh")
    if "availability_by_sample" in availability:
        return list(availability["availability_by_sample"].items())
    return [(None, availability)]


def _priors(result):
    if not isinstance(result.priors, dict) or not {"enzyme", "interaction"} <= result.priors.keys():
        raise ValueError("Explanation requires prior snapshots; rerun analysis with this version")
    return result.priors["enzyme"], result.priors["interaction"]


def _definitions(result):
    enzyme, _ = _priors(result)
    parsed = _normalize_enzyme_metabolite(enzyme)
    retained = _get_reaction_gene_sets(parsed)
    kept = {(str(r.hmdb_id), str(r.reaction), str(r.direction)): frozenset(r.genes)
            for r in retained.itertuples()}
    rows = []
    for key, part in parsed.groupby(["hmdb_id", "reaction", "direction"], sort=False):
        genes = sorted(set(part["gene"]))
        candidates = [(k, v) for k, v in kept.items() if k[0] == key[0] and k[2] == key[2]]
        matched = [k[1] for k, v in candidates if frozenset(genes) == v]
        supersets = [k[1] for k, v in candidates if frozenset(genes) < v]
        status = ("retained" if key in kept else
                  "duplicate_gene_set" if matched else "strict_subset")
        rows.append(dict(hmdb_id=key[0], reaction=key[1], direction=key[2],
                         genes=genes, included=status == "retained", selection_status=status,
                         represented_by=matched if matched else supersets))
    return pd.DataFrame(rows, columns=["hmdb_id", "reaction", "direction", "genes",
                                      "included", "selection_status", "represented_by"])


def _value(table, hmdb, cell_type=None):
    if table is None or len(table) == 0:
        return np.nan
    ids = table.index.get_level_values("hmdb_id")
    selected = table.loc[ids == hmdb]
    if selected.empty:
        return np.nan
    if isinstance(selected, pd.Series):
        return float(selected.iloc[0])
    if cell_type not in selected.columns:
        return np.nan
    return float(selected.iloc[0][cell_type])


def _sum_available(values):
    # Reporting an unweighted sum must not invalidate a valid weighted score.
    with np.errstate(over="ignore", invalid="ignore"):
        value = np.sum(np.asarray(values, dtype=np.longdouble))
    return float(value) if np.isfinite(value) and value <= np.finfo(float).max else np.nan


def _empty_concat(parts, columns=()):
    return pd.concat(parts, ignore_index=True) if parts else pd.DataFrame(columns=columns)


def build_report_tables(result):
    """Materialize reporting tables once, keyed by sample/type rather than event pairs."""
    enzyme, interaction = _priors(result)
    definitions = _definitions(result)
    retained = definitions.loc[definitions["included"]]
    all_hmdb = sorted(set(enzyme.loc[_valid_hmdb_mask(enzyme.hmdb_id), "hmdb_id"]) |
                      set(interaction.loc[_valid_hmdb_mask(interaction.hmdb_id), "hmdb_id"]))
    needed_genes = sorted(set(enzyme["gene"]) | set(interaction["sensor_gene"]))
    groups = {(m, d): list(g.itertuples()) for (m, d), g in
              retained.groupby(["hmdb_id", "direction"], sort=False)}
    annotated_counts = definitions.groupby(["hmdb_id", "direction"], sort=False).size().to_dict()
    reactions, genes, coverage, senders, receivers, references, production = [], [], [], [], [], [], []
    for sample, unit in _units(result):
        pb, fractions = unit.get("pseudobulk"), unit.get("expr_frac")
        if pb is None or fractions is None:
            raise ValueError("Explanation requires pseudobulk and expression-fraction intermediates")
        raw, weighted = unit.get("reaction_scores_unweighted"), unit.get("reaction_scores")
        if raw is None or weighted is None:
            raise ValueError("Reaction activities were not captured; rerun analysis with this version")
        panel = set(pb.columns)
        for cell_type in pb.index:
            for gene in needed_genes:
                present = gene in panel
                genes.append(dict(sample=sample, cell_type=cell_type, gene=gene,
                                  input_present=present,
                                  mean_expression=float(pb.at[cell_type, gene]) if present else np.nan,
                                  expr_frac=float(fractions.at[cell_type, gene]) if present else np.nan))
            for row in definitions.itertuples():
                key = (row.hmdb_id, row.reaction, row.direction)
                present = sorted(set(row.genes) & panel)
                evaluable = row.included and bool(present)
                reactions.append(dict(sample=sample, cell_type=cell_type, hmdb_id=row.hmdb_id,
                    reaction=row.reaction, direction=row.direction, genes=row.genes,
                    input_genes=present, missing_genes=sorted(set(row.genes) - panel),
                    included=row.included, selection_status=row.selection_status,
                    represented_by=row.represented_by, evaluable=evaluable,
                    activity_unweighted=float(raw.at[key, cell_type]) if evaluable else np.nan,
                    activity_weighted=float(weighted.at[key, cell_type]) if evaluable else np.nan))
            for hmdb in all_hmdb:
                for direction, name in [("product", "P"), ("substrate", "C"), ("exporter", "E")]:
                    group = groups.get((hmdb, direction), [])
                    sets = [r.genes for r in group]
                    universe = set().union(*sets) if sets else set()
                    counts = [len(set(g) & panel) for g in sets]
                    activities = [float(raw.at[(hmdb, r.reaction, direction), cell_type])
                                  for r in group]
                    state = ("prior_missing" if not sets else "input_genes_unavailable" if not any(counts)
                             else "positive_expression" if any(x > 0 for x in activities)
                             else "measured_zero")
                    coverage.append(dict(sample=sample, cell_type=cell_type, hmdb_id=hmdb, direction=name,
                        prior_present=bool(sets), n_prior_genes=len(universe),
                        n_input_genes=len(universe & panel),
                        input_gene_fraction=len(universe & panel) / len(universe) if universe else np.nan,
                        missing_genes=sorted(universe - panel),
                        n_annotated_reactions=annotated_counts.get((hmdb, direction), 0),
                        n_reactions=len(sets),
                        n_fully_covered=sum(c == len(g) for c, g in zip(counts, sets)),
                        n_partially_covered=sum(0 < c < len(g) for c, g in zip(counts, sets)),
                        n_uncovered=sum(c == 0 for c in counts),
                        n_positive_reactions=sum(x > 0 for x in activities), expression_state=state))
                availability = _value(unit.get("availability"), hmdb, cell_type)
                if not np.isnan(availability):
                    entry = dict(sample=sample, hmdb_id=hmdb, sender=cell_type,
                        cell_fraction=float(unit["cell_fractions"].loc[cell_type]),
                        abundance_weight=float(unit["sender_abundance_weights"].loc[cell_type]),
                        sender_score=availability)
                    for name, direction in [("P", "product"), ("C", "substrate"), ("E", "exporter")]:
                        group = groups.get((hmdb, direction), [])
                        entry[name + "_unweighted"] = _sum_available(
                            [raw.at[(hmdb, r.reaction, direction), cell_type] for r in group])
                        entry[name] = _value(unit.get(name), hmdb, cell_type)
                        entry[name + "_score"] = _value(unit.get(name + "_score"), hmdb, cell_type)
                        entry[name + "_ref"] = _value(unit.get(name + "_ref"), hmdb)
                    for name in ["base_availability", "E_effective", "E_factor"]:
                        entry[name] = _value(unit.get(name), hmdb, cell_type)
                    metadata = unit["metadata"]
                    selected = metadata.loc[metadata.index.get_level_values("hmdb_id") == hmdb]
                    for name in ["production_status", "consumption_status", "export_status"]:
                        entry[name] = selected.iloc[0][name]
                    senders.append(entry)
        for name in ["P", "C", "E"]:
            for idx, row in unit[name].iterrows():
                positive = row.index[row > 0].tolist()
                references.append(dict(sample=sample, feature_kind=name, feature_id=idx[1],
                    reference=_value(unit[name + "_ref"], idx[1]),
                    method=unit["pce_reference"], cell_types=positive,
                    n_celltypes=len(positive), positive_values=row.loc[positive].tolist()))
        receiver_components = unit.get("receiver_components", pd.DataFrame())
        if not receiver_components.empty:
            for gene, group in receiver_components.groupby("sensor_gene", sort=False):
                first = group.iloc[0]
                positive = group.loc[group.mean_expression > 0]
                references.append(dict(sample=sample, feature_kind="receiver", feature_id=gene,
                    reference=first.reference, method=first.reference_method,
                    cell_types=positive.receiver.tolist(), n_celltypes=len(positive),
                    positive_values=positive.mean_expression.tolist()))
        # Preserve unmeasured receptors from the original prior as explicit NAs.
        for row in interaction.drop_duplicates(["hmdb_id", "sensor_gene"]).itertuples():
            for cell_type in pb.index:
                entry = dict(sample=sample, hmdb_id=row.hmdb_id, sensor_gene=row.sensor_gene,
                             receiver=cell_type, input_present=row.sensor_gene in panel)
                match = (receiver_components.loc[(receiver_components.sensor_gene == row.sensor_gene) &
                                                  (receiver_components.receiver == cell_type)]
                         if not receiver_components.empty else pd.DataFrame())
                if not match.empty:
                    entry.update(match.iloc[0].to_dict())
                receivers.append(entry)
        diagnostic = unit.get("production_diagnostics")
        if isinstance(diagnostic, pd.DataFrame):
            diagnostic = diagnostic.reset_index()
            diagnostic.insert(0, "sample", sample)
            production.append(diagnostic)
    coverage_frame = pd.DataFrame(coverage)
    summary = []
    if not coverage_frame.empty:
        for direction, group in coverage_frame.groupby("direction", sort=False):
            annotated = group.loc[group.prior_present]
            summary.append(dict(direction=direction, n_units=len(group),
                n_units_with_prior=len(annotated),
                n_units_with_input_genes=int((annotated.n_input_genes > 0).sum()),
                n_units_with_positive_expression=int((group.expression_state == "positive_expression").sum()),
                denominator="sample x observed cell type x prior-union HMDB; retained reactions"))
    diagnostics = result.sample_events if result.sample_events is not None else result.events
    status = diagnostics.copy()
    status.attrs = {}
    if "cell_mesh_score" in status:
        status["evaluation_status"] = np.where(status.cell_mesh_score.isna(), "not_evaluable",
                                        np.where(status.cell_mesh_score == 0, "evaluated_zero", "evaluated_positive"))
        status["reason"] = np.where(status.cell_mesh_score.isna(), "sample_endpoint_unobserved", "")
    return {
        "explanation": {"sender_components": pd.DataFrame(senders),
                        "receiver_components": pd.DataFrame(receivers),
                        "reaction_details": pd.DataFrame(reactions),
                        "gene_details": pd.DataFrame(genes),
                        "reference_context": pd.DataFrame(references)},
        "coverage": {"summary": pd.DataFrame(summary), "details": coverage_frame},
        "diagnostics": {"production_status": _empty_concat(production), "event_status": status},
    }


def coverage_report(result):
    """Return coverage counts and explicit denominators; absent priors have NA coverage."""
    return build_report_tables(result)["coverage"]


def _select(frame, **criteria):
    if frame.empty:
        return frame.copy()
    mask = pd.Series(True, index=frame.index)
    for key, value in criteria.items():
        mask &= frame[key].isna() if value is None else frame[key].eq(value)
    return frame.loc[mask].copy()


def _explain(result, tables, *, hmdb_id, sensor_gene, sender, receiver, sample=None, observed=None, aggregate_events=None):
    hmdb = _normalize_hmdb_id(hmdb_id)
    if hmdb is None or pd.isna(hmdb):
        raise ValueError("hmdb_id must be a valid HMDB identifier")
    sensor_gene, sender, receiver = (str(x).strip() for x in (sensor_gene, sender, receiver))
    units = _units(result)
    if sample is not None:
        sample = str(sample).strip()
        units = [(s, u) for s, u in units if s == sample]
        if not units:
            raise ValueError(f"Unknown sample {sample!r}; pooled results have no sample-level decomposition")
    _, interaction = _priors(result)
    prior = interaction.loc[(interaction.hmdb_id == hmdb) & (interaction.sensor_gene == sensor_gene)]
    if observed is None:
        observed = result.sample_events if result.sample_events is not None else result.events
    records = []
    for unit_sample, unit in units:
        criteria = dict(hmdb_id=hmdb, sensor_gene=sensor_gene, sender=sender, receiver=receiver)
        if "sample" in observed:
            criteria["sample"] = unit_sample
        match = _select(observed, **criteria)
        reasons = []
        if prior.empty:
            reasons.append("interaction_prior_missing")
        if not result.parameters.get("allow_self", True) and sender == receiver:
            reasons.append("self_communication_disabled")
        pb = unit["pseudobulk"]
        if sender not in pb.index:
            reasons.append("sender_unobserved")
        if receiver not in pb.index:
            reasons.append("receiver_unobserved")
        if sensor_gene not in pb.columns:
            reasons.append("sensor_gene_unavailable")
        cov = _select(tables["coverage"]["details"], sample=unit_sample, hmdb_id=hmdb,
                      direction="P")
        if cov.empty or not cov.iloc[0].prior_present:
            reasons.append("production_prior_missing")
        elif cov.iloc[0].n_input_genes == 0:
            reasons.append("production_genes_unavailable")
        score = float(match.iloc[0].cell_mesh_score) if not match.empty else np.nan
        status = "not_evaluable" if np.isnan(score) else "evaluated_zero" if score == 0 else "evaluated_positive"
        if score == 0:
            receiver_detail = _select(tables["explanation"]["receiver_components"], sample=unit_sample,
                                     hmdb_id=hmdb, sensor_gene=sensor_gene, receiver=receiver)
            if not receiver_detail.empty:
                if bool(receiver_detail.iloc[0].get("gate_applied", False)):
                    reasons.append("receiver_expression_fraction_gate")
                if receiver_detail.iloc[0].get("mean_expression", np.nan) == 0:
                    reasons.append("receiver_measured_zero")
            if _value(unit.get("P"), hmdb, sender) == 0:
                reasons.append("zero_production_capacity")
            if not reasons:
                reasons.append("zero_score_see_components_and_numerical_diagnostics")
        entry = match.iloc[0].to_dict() if not match.empty else criteria.copy()
        entry.update(sample=unit_sample, cell_mesh_score=score, evaluation_status=status,
                     reasons=reasons if reasons else ([] if not np.isnan(score) else ["event_not_generated"]))
        records.append(entry)
    explanation = tables["explanation"]
    sample_names = [s for s, _ in units]

    def selected(name, **criteria):
        frame = _select(explanation[name], **criteria)
        return frame.loc[frame["sample"].isin(sample_names)].copy() if not frame.empty else frame

    reactions = selected("reaction_details", hmdb_id=hmdb, cell_type=sender)
    needed = {sensor_gene}
    for names in reactions.get("genes", []):
        needed.update(names)
    gene_table = explanation["gene_details"]
    gene_table = gene_table.loc[gene_table["sample"].isin(sample_names) &
                                gene_table.cell_type.isin([sender, receiver]) & gene_table.gene.isin(needed)].copy()
    cov = _select(tables["coverage"]["details"], hmdb_id=hmdb, cell_type=sender)
    cov = cov.loc[cov["sample"].isin(sample_names)] if not cov.empty else cov
    aggregate = _select(result.events if aggregate_events is None else aggregate_events, hmdb_id=hmdb, sensor_gene=sensor_gene, sender=sender, receiver=receiver)
    aggregate.attrs = {}
    enzyme, _ = _priors(result)
    return dict(identity=dict(hmdb_id=hmdb, sensor_gene=sensor_gene, sender=sender, receiver=receiver),
        aggregate_event=aggregate, sample_events=pd.DataFrame(records),
        sender_components=selected("sender_components", hmdb_id=hmdb, sender=sender),
        receiver_components=selected("receiver_components", hmdb_id=hmdb, sensor_gene=sensor_gene, receiver=receiver),
        reaction_details=reactions, gene_details=gene_table, coverage=cov,
        reference_context=explanation["reference_context"].loc[
            explanation["reference_context"]["sample"].isin(sample_names) &
            (((explanation["reference_context"].feature_kind != "receiver") &
              explanation["reference_context"].feature_id.eq(hmdb)) |
             ((explanation["reference_context"].feature_kind == "receiver") &
              explanation["reference_context"].feature_id.eq(sensor_gene)))].copy(),
        enzyme_prior=enzyme.loc[enzyme.hmdb_id == hmdb].copy(), interaction_prior=prior.copy(),
        notes=["Coverage describes genes present in the input matrix, not reaction function or confidence.",
               "Sample-aware aggregate scores are medians of event scores, not recomputed from component medians.",
               "Positive-score prevalence is not a significance or replication rate.",
               "Measured zero expression does not establish biological absence; see numerical diagnostics."])


def explain_event(result, *, hmdb_id, sensor_gene, sender, receiver, sample=None):
    return _explain(result, build_report_tables(result), hmdb_id=hmdb_id, sensor_gene=sensor_gene,
                    sender=sender, receiver=receiver, sample=sample)


def _table(frame, limit=100):
    if frame.empty:
        return "<p>无对应记录。</p>"
    shown = frame.head(limit)
    scope = f"共 {len(frame):,} 行" + (f"，预览前 {len(shown):,} 行" if len(frame) > len(shown) else "")
    return (f"<p class='report-note'>{scope}。</p>" +
            "<div class='table'>" + shown.to_html(index=False, escape=True, na_rep="NA", border=0) + "</div>")


def _panel_rows(detail):
    specs = [
        ("P/C/E normalized scores", "sender_components", ["P_score", "C_score", "E_score"]),
        ("Sender: base and exporter adjustment", "sender_components",
         ["base_availability", "E_effective", "E_factor", "sender_score"]),
        ("Receiver: before / after expression gate", "receiver_components", ["score_before_gate", "sensor_score"]),
        ("Event scores by sample", "sample_events", ["cell_mesh_score"]),
    ]
    for title, name, columns in specs:
        frame = detail[name]
        values = []
        for row in frame.head(10).to_dict("records"):
            sample = "pooled" if row["sample"] is None else str(row["sample"])
            for column in columns:
                values.append((f"{sample} / {column}", row.get(column, np.nan), columns.index(column)))
        yield title, columns, values


def _score_svg(values, title):
    """Inline SVG of captured values, including missingness; no dependencies."""
    pieces = [f'<svg xmlns="http://www.w3.org/2000/svg" role="img" aria-label="{html.escape(title)}" '
              f'viewBox="0 0 720 {max(90, 26 * len(values) + 50)}">',
              f'<text x="5" y="20" font-size="15" font-weight="bold">{html.escape(title)}</text>']
    if not values:
        pieces.append('<text x="5" y="55">NA：没有可计算的组成项</text>')
    colors = ["#277b89", "#c16b3c", "#6864a1", "#397e53"]
    for i, (label, value, color) in enumerate(values):
        y = 50 + 26 * i
        pieces.append(f'<text x="5" y="{y}" font-size="12"><title>{html.escape(label)}</title>{html.escape(label[:38])}</text>')
        if pd.isna(value):
            pieces.append(f'<text x="290" y="{y}" font-size="13">NA：无法计算</text>')
        else:
            pieces.append(f'<rect x="290" y="{y-13}" width="{340*float(value):.4f}" height="16" fill="{colors[color]}"/>')
            pieces.append(f'<text x="{300+340*float(value):.4f}" y="{y}" font-size="13">{float(value):.6g}</text>')
    return "".join(pieces) + "</svg>"


def plot_event_explanation(result, *, hmdb_id, sensor_gene, sender, receiver, sample=None, figsize=(15, 11)):
    """Return (figure, axes) for four component panels; save PNG/PDF with figure.savefig.

    Matplotlib is optional and imported only here. The first ten sample units are
    shown per panel; complete values remain available through explain_event().
    """
    import matplotlib.pyplot as plt
    detail = explain_event(result, hmdb_id=hmdb_id, sensor_gene=sensor_gene,
                           sender=sender, receiver=receiver, sample=sample)
    fig, axes = plt.subplots(2, 2, figsize=figsize, constrained_layout=True)
    colors = ["#277b89", "#c16b3c", "#6864a1", "#397e53"]
    for ax, (title, _, values) in zip(axes.flat, _panel_rows(detail)):
        for i, (label, value, color) in enumerate(values):
            if pd.isna(value):
                ax.text(0.02, i, "NA — not evaluable", va="center")
            else:
                ax.barh(i, value, color=colors[color], height=0.65)
                ax.text(float(value) + 0.015, i, f"{float(value):.4g}", va="center", fontsize=9)
        if not values:
            ax.text(0.04, 0.5, "NA — no computable component", transform=ax.transAxes)
        ax.set_yticks(range(len(values)), labels=[row[0] for row in values], fontsize=8)
        ax.invert_yaxis()
        ax.set_xlim(0, 1.18)
        ax.set_xticks([0, 0.5, 1])
        ax.set_title(title, loc="left", fontsize=11)
        ax.spines[["top", "right"]].set_visible(False)
    fig.suptitle(f"{sender} → {hmdb_id} / {sensor_gene} → {receiver}\n"
                 "First 10 sample units; E_effective can be the neutral missing-prior value", fontsize=12)
    return fig, axes


def _input_archive_html(result, input_storage):
    archive = result.input_archive if input_storage is None else input_storage["archive"]
    provenance = result.input_provenance
    parts = []
    if archive:
        parts.append("<p>表达矩阵归档：<code>" + html.escape(str(archive["filename"])) +
                     "</code>，矩阵位于 <code>.X</code>。这里只展示归档记录；调用 <code>load_input()</code> 时校验文件。</p>")
    else:
        parts.append("<p>本结果未附带细胞级输入文件。结果 H5 保存汇总信息与解释数据。</p>")
    if provenance:
        parts.append("<p>输入来源：" + html.escape(str(provenance["source"])) +
                     f"；维度：{int(provenance['n_obs']):,} 个细胞 × {int(provenance['n_vars']):,} 个基因。"
                     "保存的是调用时提供的表达值；不推定为原始 counts。</p>")
        parts.append("<details><summary>输入来源与归档信息</summary><pre>" + html.escape(
            json.dumps({"expression_input": provenance, "input_archive": archive},
                       ensure_ascii=False, indent=2, default=str)) + "</pre></details>")
    return "".join(parts)


def render_html_report(result, *, events=None, selection=None, detail_n=INHERIT, plot_top_n=20,
                       tables=None, input_storage=None, positive_col=None, positive_threshold=None,
                       positive_qc_only=None, positive_min_score=INHERIT, delivery=None, csv_paths=None,
                       adata=None, general_report=None, navigation_report=None, event_links=None,
                       selection_evidence_path=None):
    """Render a single tier; file orchestration lives in result_io."""
    from .selection import report_selection
    from ._report_plots import validate_plot_top_n
    from .positivity import resolve_positive, delivery_record
    from ._report_pages import render_overview, render_event
    plot_top_n = validate_plot_top_n(plot_top_n)
    spec = resolve_positive(result, positive_col=positive_col, positive_threshold=positive_threshold,
                            positive_qc_only=positive_qc_only, positive_min_score=positive_min_score,
                            selection=selection) if delivery is None else delivery['positive']
    _, requests = report_selection(result, selection, events, detail_n, positive=spec)
    delivery = delivery_record(result, spec, selection, requests, INHERIT, plot_top_n) if delivery is None else delivery
    if selection is None and requests:
        if len(requests) != 1:
            raise ValueError("Multiple event reports require save() or to_html() file orchestration")
        return render_event(result, requests[0], tables=tables, spec=spec, adata=adata,
                            general_report=general_report, navigation_report=navigation_report)
    return render_overview(result, selection=selection, tables=tables, delivery=delivery,
                           csv_paths=csv_paths, plot_top_n=plot_top_n, input_storage=input_storage,
                           event_links=event_links, general_report=general_report,
                           navigation_report=navigation_report, adata=adata,
                           selection_evidence_path=selection_evidence_path)


def write_html_report(result, path, *, filtered_events=None, selected_events=None, events=None, selection=None, include_all_events=False,
                      detail_n=INHERIT, plot_top_n=20, **positive_options):
    """Write HTML and linked CSV companions atomically; never rescore."""
    from .result_io import write_report_bundle
    paths = write_report_bundle(result, path, filtered_events=filtered_events, selected_events=selected_events, events=events, selection=selection,
                                include_all_events=include_all_events,
                                detail_n=detail_n, plot_top_n=plot_top_n, **positive_options)
    return paths if any(key.startswith(('event_report:', 'pair_report:')) for key in paths) else paths['report']
