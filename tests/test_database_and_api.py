import numpy as np
import pandas as pd

from cell_mesh import load_cell_mesh_database, run_cell_mesh, run_metcomm


class FakeAnnData:
    def __init__(self, X, var_names, obs):
        self.X = X
        self.layers = {}
        self.var_names = pd.Index(var_names)
        self.obs = pd.DataFrame(obs)


def test_load_packaged_database():
    enzyme, sensor = load_cell_mesh_database()
    assert not enzyme.empty
    assert not sensor.empty
    assert {"metabolite", "gene", "role"}.issubset(enzyme.columns)
    assert {"metabolite", "sensor_gene", "sensor_type"}.issubset(sensor.columns)
    assert set(sensor["sensor_type"]).issubset({"surface_receptor", "transporter", "nuclear_receptor", "intracellular_sensor"})


def test_run_cell_mesh_with_packaged_database():
    enzyme, sensor = load_cell_mesh_database()
    # Use a matched metabolite with available enzyme and sensor genes.
    met = sorted(set(enzyme["metabolite"]).intersection(sensor["metabolite"]))[0]
    e_gene = enzyme.loc[enzyme["metabolite"] == met, "gene"].iloc[0]
    s_gene = sensor.loc[sensor["metabolite"] == met, "sensor_gene"].iloc[0]
    genes = [e_gene, s_gene, "BACKGROUND"]
    X = np.array([
        [5, 0, 0], [4, 0, 0], [5, 0, 0],
        [0, 4, 0], [0, 5, 0], [0, 4, 0],
    ], dtype=float)
    adata = FakeAnnData(X, genes, {"cell_type": ["A", "A", "A", "B", "B", "B"]})
    res = run_cell_mesh(adata, cell_type_key="cell_type", min_cells_per_group=2, allow_self=False)
    assert not res.events.empty
    assert "cell_mesh_score" in res.events.columns
    assert res.events.iloc[0]["sender"] == "A"
    assert res.events.iloc[0]["receiver"] == "B"
    res2 = run_metcomm(adata, cell_type_key="cell_type", min_cells_per_group=2, allow_self=False)
    assert not res2.events.empty
