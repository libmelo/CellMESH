"""Exercise actual example cells so input preparation cannot regress independently."""
import json
from pathlib import Path
import subprocess
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

from cellmesh import load_cell_mesh_database, run_cell_mesh
from cellmesh.database import validate_priors
from cellmesh.plotting import plot_communication_network, plot_significant_event_counts

EXAMPLES = Path(__file__).resolve().parents[1] / "examples"
HNSC = ["demo_HNSC_cellmesh_workflow.ipynb", "demo_HNSC_sample_aware_workflow.ipynb"]


def cell(name, index):
    return "".join(json.loads((EXAMPLES / name).read_text())["cells"][index]["source"])


@pytest.fixture
def inputs():
    adata = AnnData(np.array([[4, 9, 1], [1, 4, 7], [8, 3, 2], [2, 5, 7]], float),
                    obs=pd.DataFrame({"cell_type": ["A", "B"] * 2, "sample": ["s1"] * 2 + ["s2"] * 2},
                                     index=["c1", "c2", "c3", "c4"]),
                    var=pd.DataFrame(index=["G1", "G2", "R"]))
    enzyme = pd.DataFrame({"metabolite": ["M"] * 4, "hmdb_id": ["HMDB1"] * 4,
                           "gene": ["G1", "MISSING", "G1", "G2"],
                           "reaction": ["partial", "partial", "pair", "pair"],
                           "role": ["production"] * 4})
    sensor = pd.DataFrame({"metabolite": ["M"], "hmdb_id": ["HMDB1"],
                           "sensor_gene": ["R"], "sensor_type": ["Transporter"]})
    return adata, *load_cell_mesh_database(enzyme, sensor)


@pytest.mark.parametrize("name", HNSC)
@pytest.mark.parametrize("mode", ["pooled_stratified", "sample_aware"])
def test_example_selection_preserves_complete_reactions(inputs, name, mode):
    adata, enzyme, sensor = inputs
    ns = dict(adata=adata, pd=pd, display=lambda *args: None,
              load_cell_mesh_database=lambda: (enzyme.copy(), sensor.copy()))
    exec(cell(name, 6), ns)
    assert "MISSING" in set(ns["enzyme_demo"].gene)
    options = dict(min_cells=1, sample_mode=mode, sample_key="sample", n_perms=5, random_state=7)
    selected = run_cell_mesh(adata, ns["enzyme_demo"], ns["sensor_demo"], **options)
    direct = run_cell_mesh(adata, enzyme, sensor, **options)
    pd.testing.assert_frame_equal(selected.events, direct.events, check_exact=True)
    availability = selected.availability_results
    units = availability["availability_by_sample"].values() if mode == "sample_aware" else [availability]
    for unit in units:
        assert unit["metadata"]["n_product_reactions"].tolist() == [2]


def test_comprehensive_trace_matches_main_analysis(inputs):
    adata, enzyme, sensor = inputs
    ns = dict(adata=adata, enzyme_metabolite=enzyme, metabolite_sensor=sensor,
              validate_priors=validate_priors, pd=pd, np=np)
    exec(cell("cellmesh_comprehensive_walkthrough.ipynb", 7), ns)
    assert "MISSING" in set(ns["enzyme_prior"].gene)
    main = run_cell_mesh(adata, enzyme, sensor, min_cells=1)
    trace = run_cell_mesh(adata, ns["enzyme_prior"], ns["sensor_prior"], min_cells=1)
    pd.testing.assert_frame_equal(main.sender_scores, trace.sender_scores, check_exact=True)


@pytest.mark.parametrize("threshold,empty", [(0.0, True), (1.0, False)])
def test_example_network_empty_and_nonempty_selection(inputs, threshold, empty, capsys):
    adata, enzyme, sensor = inputs
    result = run_cell_mesh(adata, enzyme, sensor, min_cells=1, n_perms=5, random_state=7)
    overview = plot_significant_event_counts(result, max_fdr=threshold)
    ns = dict(res=result, filtered_events=overview["filtered_events"],
              min_cell_mesh_score=None, max_perm_pvalue=None, max_fdr=threshold,
              qc_only=True, plot_communication_network=plot_communication_network,
              plt=plt, display=lambda *args: None)
    try:
        exec(cell(HNSC[0], 15), ns)
        assert ns["max_fdr"] == threshold
        assert ns["filtered_events"].empty == empty
        if empty:
            assert ns["communication_network"] is None
            assert "Thresholds are unchanged" in capsys.readouterr().out
        else:
            assert not ns["edge_table"].empty
    finally:
        plt.close("all")


@pytest.mark.parametrize("path", sorted(EXAMPLES.glob("*.ipynb")), ids=lambda p: p.stem)
def test_notebook_executes_in_fresh_process(path):
    pytest.importorskip("IPython", reason="Install cellmesh[notebook] to execute notebook examples")
    # HNSC data and the default databases intentionally remain local during testing.
    # Only dependent examples are skipped in installations lacking those fixtures.
    data_dir = EXAMPLES.parent / "cellmesh/data"
    if path.name.startswith("demo_HNSC"):
        try:
            load_cell_mesh_database()
        except FileNotFoundError:
            pytest.skip("Local prior files are required for this example")
        if not (data_dir / "demo_HNSC_200cell.h5ad").exists():
            pytest.skip("Local HNSC fixture is required")
    if path.name == "cellmesh_comprehensive_walkthrough.ipynb":
        for name in ["test_single_cell.h5ad", "Enzyme_new.csv", "Interaction1.0.csv"]:
            if not (data_dir / name).exists():
                pytest.skip(f"Local walkthrough fixture is required: {name}")
    script = '''
import json, sys
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from IPython.display import display
ns = {"__name__": "__main__", "display": display}
path = Path(sys.argv[1])
for index, cell in enumerate(json.loads(path.read_text())["cells"]):
    if cell["cell_type"] == "code":
        try:
            exec(compile("".join(cell["source"]), f"{path.name}:cell_{index}", "exec"), ns)
            for number in plt.get_fignums():
                plt.figure(number).canvas.draw()
        finally:
            plt.close("all")
'''
    completed = subprocess.run([sys.executable, "-c", script, str(path)], cwd=EXAMPLES.parent,
                               capture_output=True, text=True, timeout=300)
    assert completed.returncode == 0, completed.stdout + completed.stderr
