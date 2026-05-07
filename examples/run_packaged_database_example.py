import numpy as np
import pandas as pd

from cell_mesh import load_cell_mesh_database, run_cell_mesh


class FakeAnnData:
    """Small AnnData-like object for a dependency-light example."""

    def __init__(self, X, var_names, obs):
        self.X = X
        self.layers = {}
        self.var_names = pd.Index(var_names)
        self.obs = pd.DataFrame(obs)


enzyme_prior, sensor_prior = load_cell_mesh_database()

# Build a toy expression matrix using genes present in the packaged priors.
genes = sorted(set(enzyme_prior["gene"]).union(sensor_prior["sensor_gene"]))[:40]
if "PTGR1" not in genes:
    genes.append("PTGR1")
if "LTB4R" not in genes:
    genes.append("LTB4R")

rng = np.random.default_rng(7)
cell_types = ["Sender"] * 30 + ["Receiver"] * 30
X = rng.poisson(0.2, size=(60, len(genes))).astype(float)
gene_idx = {g: i for i, g in enumerate(genes)}

if "PTGR1" in gene_idx:
    X[:30, gene_idx["PTGR1"]] += rng.poisson(3, size=30)
if "LTB4R" in gene_idx:
    X[30:, gene_idx["LTB4R"]] += rng.poisson(3, size=30)

adata = FakeAnnData(X, genes, {"cell_type": cell_types})
res = run_cell_mesh(adata, cell_type_key="cell_type", n_perms=0, allow_self=False)
print(res.events.head().to_string(index=False))
