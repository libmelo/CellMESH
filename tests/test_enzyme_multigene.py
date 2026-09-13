import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

from cellmesh import compute_metabolite_availability, load_cell_mesh_database, run_cell_mesh
from cellmesh.core import EVENT_KEY_COLUMNS
from cellmesh.database import normalize_enzyme_database, validate_priors
from cellmesh.permutation import CompiledPermutationScorer


@pytest.fixture
def multigene_case():
    adata = AnnData(
        X=np.array([
            [4, 9, 1, 3, 2, 5, 1, 1], [4, 9, 2, 2, 4, 3, 1, 0],
            [1, 4, 5, 1, 1, 2, 1, 6], [1, 4, 3, 2, 0, 4, 1, 8],
            [8, 3, 2, 1, 3, 2, 2, 0], [6, 5, 1, 1, 5, 1, 2, 2],
            [2, 2, 4, 3, 0, 3, 2, 7], [0, 4, 6, 1, 2, 1, 2, 5],
        ], dtype=float),
        obs=pd.DataFrame(
            {"cell_type": ["A", "A", "B", "B"] * 2,
             "sample": ["D1"] * 4 + ["D2"] * 4},
            index=[f"cell_{i}" for i in range(8)],
        ),
        var=pd.DataFrame(index=["P1", "P2", "C1", "C2", "E1", "E2", "BASE", "S"]),
    )
    enzyme = pd.DataFrame({
        "metabolite": ["代谢物 M"] * 7,
        "hmdb_id": ["HMDB0000001"] * 7,
        "reaction": ["001", "001", "002", "002", "003", "003", "004"],
        "gene": ["P1", "P2", "C1", "C2", "E1", "E2", "BASE"],
        "role": ["production"] * 2 + ["degradation"] * 2 + ["export"] * 2 + ["production"],
        "evidence_level": ["curated"] * 7,
        "source": ["user prior"] * 7,
        "reference": ["00123"] * 7,
    })
    sensor = pd.DataFrame({
        "metabolite": ["代谢物 M"], "hmdb_id": ["HMDB0000001"],
        "sensor_gene": ["S"], "sensor_type": ["Transporter"],
    })
    return adata, enzyme, sensor


def _schema(enzyme, schema):
    out = enzyme.copy(deep=True)
    if schema != "role":
        out["direction"] = out["role"].map({
            "production": "product", "degradation": "substrate", "export": "export",
        })
        if schema != "mixed":
            out = out.drop(columns="role")
        if schema == "legacy":
            out = out.rename(columns={
                "metabolite": "standard_metName", "hmdb_id": "HMDB_ID",
                "reaction": "Reactions", "gene": "Gene_name", "direction": "Direction",
            })
    return out


@pytest.mark.parametrize("schema", ["role", "direction", "mixed", "legacy"])
@pytest.mark.parametrize("gene_field", [
    " P1 ; P2 ; P1 ", " P1, P2, P1 ", " P1| P2| P1 ",
    "P1[reviewed, curated]; P2[database|reviewed] | P1[reviewed, curated]",
])
def test_multigene_fields_normalize_before_expression_matching(multigene_case, schema, gene_field):
    adata, enzyme, sensor = multigene_case
    prior = enzyme.iloc[[0]].copy()
    prior["gene"] = gene_field
    # Expansion must retain typed/custom provenance instead of rebuilding all
    # metadata columns from inferred Python row values.
    prior["confidence"] = pd.Series([0.75], dtype="float32")
    prior = _schema(prior, schema)
    snapshot = prior.copy(deep=True)
    normalized = normalize_enzyme_database(prior)
    expected = enzyme.iloc[:2].copy()
    expected["confidence"] = pd.Series([0.75, 0.75], dtype="float32")
    pd.testing.assert_frame_equal(normalized[expected.columns], expected)
    pd.testing.assert_frame_equal(normalize_enzyme_database(normalized), normalized)
    matched, _ = validate_priors(normalized, sensor, adata.var_names)
    assert matched["gene"].tolist() == ["P1", "P2"]
    pd.testing.assert_frame_equal(prior, snapshot)


def test_multigene_expansion_accepts_categorical_gene_columns(multigene_case):
    _, enzyme, _ = multigene_case
    prior = enzyme.iloc[[0, 2]].copy()
    prior.index = ["duplicate", "duplicate"]
    prior["gene"] = pd.Categorical(["P1; P2", "C1; C2"])
    prior["source"] = pd.Categorical(prior["source"])
    normalized = normalize_enzyme_database(prior)
    assert normalized["gene"].tolist() == ["P1", "P2", "C1", "C2"]
    assert normalized["reaction"].tolist() == ["001", "001", "002", "002"]
    assert normalized["source"].dtype == prior["source"].dtype


@pytest.mark.parametrize("schema", ["role", "direction", "mixed", "legacy"])
@pytest.mark.parametrize("csv_input", [False, True])
@pytest.mark.parametrize("sample_mode", ["pooled_stratified", "sample_aware"])
def test_multigene_analysis_matches_explicit_rows(
    tmp_path, multigene_case, schema, csv_input, sample_mode,
):
    adata, enzyme, sensor = multigene_case
    composite = enzyme.iloc[[0, 2, 4, 6]].copy()
    composite["gene"] = [
        "P1[Enzyme]; P2[Enzyme]; P1[Enzyme]", "C1[Enzyme], C2[Enzyme]",
        "E1[Transporter]|E2[Transporter]", "BASE[Enzyme]",
    ]
    composite = _schema(composite, schema)
    snapshot = composite.copy(deep=True)
    prior_input = composite
    if csv_input:
        prior_input = tmp_path / "enzyme.csv"
        composite.to_csv(prior_input, index=False)
    loaded, _ = load_cell_mesh_database(prior_input, sensor)
    pd.testing.assert_frame_equal(loaded[enzyme.columns], enzyme)
    kwargs = dict(sample_key="sample", sample_mode=sample_mode, min_cells=1,
                  min_expr_frac=0.25, n_perms=4, random_state=11,
                  store_null_scores=True)
    reference = run_cell_mesh(adata, enzyme, sensor, **kwargs)
    actual = run_cell_mesh(adata, prior_input, sensor, **kwargs)
    assert not actual.events.empty
    for name in ("events", "sender_scores", "receiver_scores", "sample_events"):
        expected_table = getattr(reference, name)
        if expected_table is not None:
            pd.testing.assert_frame_equal(getattr(actual, name), expected_table, check_exact=True)

    if sample_mode == "sample_aware":
        pd.testing.assert_frame_equal(
            actual.events.attrs["sample_aware_null_scores"],
            reference.events.attrs["sample_aware_null_scores"], check_exact=True,
        )
        availability = actual.availability_results["availability_by_sample"]
        units = {sample: adata[adata.obs["sample"].eq(sample)].copy() for sample in availability}
    else:
        availability = {"pooled": actual.availability_results}
        units = {"pooled": adata}
    for unit, main in availability.items():
        direct = compute_metabolite_availability(units[unit], composite, min_cells=1)
        for name in ("P", "C", "E", "availability", "metadata"):
            pd.testing.assert_frame_equal(direct[name], main[name], check_exact=True)
        assert main["metadata"]["n_product_reactions"].tolist() == [2]
        assert main["metadata"]["n_substrate_reactions"].tolist() == [1]
        assert main["metadata"]["n_exporter_reactions"].tolist() == [1]

        # Independent numeric check: the production pair is one reaction,
        # duplicates do not add capacity, and BASE is a separate reaction.
        means = pd.DataFrame(units[unit].X, columns=units[unit].var_names).groupby(
            units[unit].obs["cell_type"].to_numpy(), sort=False,
        ).mean()
        expected_p = (np.sqrt((means["P1"] + 1) * (means["P2"] + 1)) - 1 + means["BASE"]) * 0.5
        np.testing.assert_allclose(main["P"].iloc[0], expected_p.reindex(main["P"].columns), rtol=1e-12)
    pd.testing.assert_frame_equal(composite, snapshot)


@pytest.mark.parametrize("schema", ["role", "direction"])
def test_inline_gene_evidence_is_retained_without_inventing_user_source(multigene_case, schema):
    _, enzyme, _ = multigene_case
    prior = enzyme.iloc[[0]].drop(columns=["evidence_level", "source"]).copy()
    prior["gene"] = "P1[reviewed, curated]; P2[database|reviewed]"
    normalized = normalize_enzyme_database(_schema(prior, schema))
    assert normalized["gene"].tolist() == ["P1", "P2"]
    assert normalized["evidence_level"].tolist() == ["reviewed, curated", "database|reviewed"]
    if schema == "role":
        assert "source" not in normalized
    pd.testing.assert_frame_equal(normalize_enzyme_database(normalized), normalized)


@pytest.mark.parametrize("schema", ["role", "direction"])
@pytest.mark.parametrize("gene_field", [None, np.nan, pd.NA, "", " ; , | "])
def test_empty_gene_fields_are_ignored_consistently(multigene_case, schema, gene_field):
    _, enzyme, _ = multigene_case
    prior = enzyme.iloc[[0]].copy()
    prior["gene"] = gene_field
    normalized = normalize_enzyme_database(_schema(prior, schema))
    assert normalized.empty
    assert set(enzyme.columns).issubset(normalized.columns)


@pytest.mark.parametrize("schema", ["role", "direction"])
@pytest.mark.parametrize("available", [["P1", "P2"], ["P1"], []])
def test_partial_gene_matching_keeps_measured_members(multigene_case, schema, available):
    _, enzyme, sensor = multigene_case
    prior = enzyme.iloc[[0]].copy()
    prior["gene"] = "P1[Enzyme]; P2[Enzyme]; NOT_MEASURED[Unknown]"
    # validate_priors is also callable with a standard role table directly.
    if schema == "direction":
        prior = normalize_enzyme_database(_schema(prior, schema))
    validated, _ = validate_priors(prior, sensor, available + ["S"])
    assert validated["gene"].tolist() == available


@pytest.mark.parametrize("sample_mode", ["pooled_stratified", "sample_aware"])
@pytest.mark.parametrize("schema", ["role", "direction"])
def test_compiled_multigene_scores_match_public_recomputation(multigene_case, sample_mode, schema):
    adata, enzyme, sensor = multigene_case
    composite = enzyme.iloc[[0, 2, 4, 6]].copy()
    composite["gene"] = ["P1;P2;P1", "C1,C2", "E1|E2", "BASE[Enzyme]"]
    composite = _schema(composite, schema)
    kwargs = dict(sample_key="sample", sample_mode=sample_mode, min_cells=1,
                  min_expr_frac=0.25, n_perms=0, receiver_reference="mean")
    observed = run_cell_mesh(adata, enzyme, sensor, **kwargs).events
    scorer = CompiledPermutationScorer(
        adata, composite, sensor, observed,
        cell_type_key="cell_type", sample_key="sample", sample_mode=sample_mode,
        layer=None, min_expr_frac=0.25, sender_abundance_exponent=1.0,
        pce_reference="mean", export_weight=0.2, receiver_reference="mean",
        prior_role_coverage=None,
    )
    rng = np.random.default_rng(31)
    observed_keys = pd.MultiIndex.from_frame(observed[EVENT_KEY_COLUMNS])
    for labels in (scorer.original_labels, scorer.permute_codes(rng), scorer.permute_codes(rng)):
        permuted = adata.copy()
        permuted.obs["cell_type"] = labels
        events = run_cell_mesh(permuted, enzyme, sensor, **kwargs).events
        expected = events.set_index(EVENT_KEY_COLUMNS)["cell_mesh_score"].reindex(
            observed_keys, fill_value=0.0,
        )
        np.testing.assert_allclose(scorer.score(labels), expected, rtol=1e-12, atol=1e-14)
