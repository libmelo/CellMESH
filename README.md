# CELL MESH

**CELL MESH** = **Metabolite-mediated Event Scoring with Sensor Hierarchies**.

CELL MESH is a prototype Python package for inferring metabolite-mediated cell-cell communication events from single-cell expression data. It scores directed events of the form:

```text
sender cell type -> receiver cell type : metabolite -> sensor
```

The package includes the two uploaded prior databases as built-in resources:

- `cell_mesh/data/enzyme_test.csv`: enzyme/reaction/metabolite table
- `cell_mesh/data/interaction_test.csv`: metabolite-receptor/sensor interaction table

These files are automatically normalized into the two prior tables used by the algorithm.

## Install

```bash
pip install -e .
```

## Basic usage with packaged database

```python
from cell_mesh import run_cell_mesh

res = run_cell_mesh(
    adata,
    cell_type_key="cell_type",
    sample_key="sample",      # optional
    layer="lognorm",          # optional
    n_perms=1000,
    min_expr_frac=0.10,
)

res.events.head()
```

If `enzyme_metabolite` and `metabolite_sensor` are not supplied, CELL MESH automatically loads the packaged databases.

## Inspect the packaged database

```python
from cell_mesh import load_cell_mesh_database

enzyme_metabolite, metabolite_sensor = load_cell_mesh_database()

print(enzyme_metabolite.head())
print(metabolite_sensor.head())
```

The normalized enzyme table has columns:

```text
metabolite, hmdb_id, gene, role, weight, evidence_level, source, reaction
```

The normalized sensor table has columns:

```text
metabolite, hmdb_id, sensor_gene, sensor_type, weight, evidence_level,
source, protein_name, reference
```

## Supported sensor hierarchy

CELL MESH maps the interaction database annotations to the following sensor types:

| Input annotation | CELL MESH sensor type |
|---|---|
| Cell surface receptor | `surface_receptor` |
| Other receptor | `surface_receptor` |
| Transporter | `transporter` |
| Nuclear receptor | `nuclear_receptor` |
| Intracellular sensor | `intracellular_sensor` |

## Main outputs

`res.events` contains one row per event. Important columns include:

- `sender`, `receiver`
- `metabolite`, `hmdb_id`
- `sensor_gene`, `sensor_type`
- `sender_score`
- `receiver_score`
- `cell_mesh_score`
- `communication_score`, kept as a backward-compatible alias
- `perm_pvalue`, `fdr`
- `confidence_tier`

`res.sender_scores` is a metabolite-by-cell-type matrix.

`res.receiver_scores` is a metabolite-sensor-receiver table.

`res.role_scores` stores intermediate role score matrices for:

- `production`
- `degradation`
- `export`
- `import`
- `usage`

## Backward compatibility

The earlier prototype function name `run_metcomm` is kept as an alias:

```python
from cell_mesh import run_metcomm
res = run_metcomm(adata, cell_type_key="cell_type")
```

New code should use `run_cell_mesh`.

## Notes on the uploaded databases

The uploaded enzyme file has `Direction` values such as `product` and `substrate`. CELL MESH maps these to:

- `product` -> `production`
- `substrate` -> `degradation`

The current uploaded enzyme database does not explicitly encode export/import/usage roles. Those role score matrices are still present but will be zero unless the user supplies additional priors with those roles.

## Caveats

- CELL MESH is currently transcriptomics-only. It estimates metabolite availability using expression proxies and prior knowledge.
- Direct metabolomics, spatial proximity, perturbation data, or downstream response signatures should be used to validate high-confidence predictions.
- Transporter and nuclear receptor events are mechanistically more constrained than surface receptor events, but still rely on indirect expression evidence.
