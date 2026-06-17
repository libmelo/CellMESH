# CellMesh
Metabolite-mediated Event Scoring with Sensor Hierarchies

**CELL MESH** = **Metabolite-mediated Event Scoring with Sensor Hierarchies**.

CELL MESH is a Python package for inferring metabolite-mediated cell-cell communication events from single-cell expression data, **fully based on the metabolite availability algorithm**.

## Algorithm Overview

### 1. Enzyme Prior Normalization
The enzyme-metabolite prior table is the canonical enzyme input. `run_cell_mesh()`
validates this prior once with `validate_priors()` and passes it directly into
the availability calculation. Internally, enzyme roles are normalized to the
availability directions used to build the P/C/E matrices:
| Role in prior | Direction | Matrix | Meaning |
|---|---|---|---|
| `production` | `product` | P | Metabolite production ability |
| `degradation` | `substrate` | C | Metabolite consumption ability |
| `export` | `exporter` | E | Metabolite efflux/transport ability |

### 2. Metabolite Availability Calculation (Unchanged)
The metabolite availability (sender score) is computed as:
```
availability = P_norm * ((1 - C_norm) ** beta) * (0.8 + 0.2 * E_norm)
```
Where:
- `P_norm`: Normalized production score [0, 1]
- `C_norm`: Normalized consumption score [0, 1]
- `E_norm`: Normalized efflux score [0, 1]
- Result range: Strictly between [0, 1], higher value indicates stronger ability to release the metabolite

### 3. Sensor Score Calculation
Sensor score is based on robust min-max normalized sensor gene expression:
- First, compute pseudobulk mean expression of each sensor gene per cell type
- Then apply robust min-max normalization across cell types for each gene
- Genes with `sensor_expr_frac < min_expr_frac` get `sensor_score = 0`

### 4. Communication Score
Communication score is the geometric mean of metabolite availability and sensor score:
```
cell_mesh_score = sqrt(metabolite_availability * sensor_score)
```
Events are matched by `hmdb_id` on both the sender availability side and receiver sensor side. Prior rows without `hmdb_id` are excluded before event construction.

### 5. Sensor Types
Sensors are categorized into three types based on the `Annotation` column in `interaction_test.csv`:
- `Cell surface receptor`
- `Transporter`
- `Other receptor`

Permutation p-values and FDR are computed separately within each sensor type.

## Included Data
The package includes built-in prior databases and a small AnnData example file:

- `cellmesh/data/Enzyme1.0.csv`: versioned enzyme/reaction/metabolite table
- `cellmesh/data/Interaction1.0.csv`: versioned metabolite-sensor interaction table
- `cellmesh/data/enzyme_test.csv`: legacy small enzyme prior used for compatibility tests
- `cellmesh/data/interaction_test.csv`: legacy small sensor prior used for compatibility tests
- `cellmesh/data/Enzyme_new.csv`: walkthrough/test enzyme prior
- `cellmesh/data/test_single_cell.h5ad`: walkthrough/test single-cell data

`load_cell_mesh_database()` automatically loads the highest packaged enzyme file
named `Enzyme<version>.csv` and the highest packaged interaction file named
`Interaction<version>.csv` if file paths are not explicitly provided. The two
version numbers do not need to match: if the highest enzyme prior is
`Enzyme1.2.csv` and the highest interaction prior is `Interaction2.2.csv`, those
two files are selected together. With the current packaged files, the default is
`Enzyme1.0.csv` and `Interaction1.0.csv`. If no versioned files are available,
the loader falls back to the legacy `enzyme_test.csv` / `interaction_test.csv`
pair. The comprehensive walkthrough notebook uses `Enzyme_new.csv`,
`Interaction1.0.csv`, and `test_single_cell.h5ad` so its calculations are fully
reproducible from packaged files.

## Install

From a local checkout:

```bash
pip install -e .
```

For running the walkthrough notebook:

```bash
pip install -e ".[notebook]"
```

For development and tests:

```bash
pip install -e ".[dev]"
pytest -q
```

10X directory loading uses Scanpy and can be enabled with:

```bash
pip install -e ".[scanpy]"
```

## Basic Usage

```python
from cellmesh import run_cell_mesh

res = run_cell_mesh(
    adata,
    cell_type_key="cell_type",
    sample_key="sample",      # optional, for within-sample permutation
    layer="lognorm",          # optional, use specific expression layer
    n_perms=1000,             # optional, number of permutations for p-value calculation
    min_expr_frac=0.10,       # optional, minimum expression fraction for sensor genes
    allow_self=True,          # optional, allow self-communication events
)

# View results
res.events.head()  # All communication events
res.sender_scores  # Metabolite × cell type availability matrix
res.receiver_scores  # Receiver/sensor scores
res.availability_results  # All intermediate calculation results
```

## Full API Reference

### `run_cell_mesh()` Parameters

| Parameter | Default | Description |
|---|---|---|
| `adata` | *required* | AnnData object containing single-cell expression data |
| `enzyme_metabolite` | `None` | Enzyme-metabolite prior table. If None, uses built-in database |
| `metabolite_sensor` | `None` | Metabolite-sensor prior table. If None, uses built-in database |
| `cell_type_key` | `"cell_type"` | Column name in adata.obs containing cell type annotations |
| `sample_key` | `None` | Column name in adata.obs containing sample annotations (for permutation) |
| `layer` | `None` | Name of expression layer to use. If None, uses adata.X |
| `min_expr_frac` | `0.05` | Minimum fraction of cells expressing a sensor gene to be considered |
| `allow_self` | `True` | Whether to allow self-communication events (sender == receiver) |
| `n_perms` | `0` | Number of permutations for empirical p-value calculation. 0 = no permutation |
| `random_state` | `0` | Random seed for reproducibility |
| `lower` | `5` | Lower percentile for robust min-max normalization (for both availability and sensor scoring) |
| `upper` | `95` | Upper percentile for robust min-max normalization (for both availability and sensor scoring) |
| `eps` | `0.05` | Small constant used in robust min-max denominators for P/C/E normalization |
| `beta` | `0.5` | Exponent weight for consumption term in availability formula |
| `missing_C_norm` | `0.2` | Default C_norm value when no consumption evidence exists |
| `missing_E_norm` | `0.5` | Default E_norm value when no efflux evidence exists |
| `min_cells` | `100` | Minimum number of cells per cell type to be included |

### Removed Parameters
The following parameters are no longer available (old scoring mechanism removed):
- `beta_sensor`
- `beta_specificity`

## Inspect the Packaged Database

```python
from cellmesh import load_cell_mesh_database

enzyme_metabolite, metabolite_sensor = load_cell_mesh_database()

print(enzyme_metabolite.head())
print(metabolite_sensor.head())
```

### Enzyme Table Columns
```
metabolite, hmdb_id, gene, role, weight, evidence_level, source, reaction
```

`compute_metabolite_availability()` also accepts this standard
`enzyme_metabolite` schema directly. Legacy direction-style inputs are only
kept as a compatibility path for low-level availability tests.

### Sensor Table Columns (from interaction_test.csv / Interaction1.0.csv)
```
ID, HMDB_ID, standard_metName, Gene_name, Protein_name, Annotation, Database source, Reference
```

## Supported Sensor Types

From the `Annotation` column in the packaged interaction CSV:
- `Cell surface receptor`
- `Transporter`
- `Other receptor` (includes nuclear receptors, intracellular sensors, and other annotations)

## Main Outputs

### `res.events`
Contains one row per communication event. Important columns:
- `sender`, `receiver`: Cell type pair
- `metabolite`, `hmdb_id`: Metabolite information
- `sensor_gene`, `sensor_type`: Sensor information (sensor_type is one of "Cell surface receptor", "Transporter", "Other receptor")
- `metabolite_availability`: Metabolite availability in sender cell type [0, 1]
- `sensor_score`: Robust min-max normalized sensor expression in receiver cell type [0, 1]
- `sensor_expr_frac`: Fraction of cells in receiver cell type expressing the sensor gene
- `cell_mesh_score`: Geometric mean of availability and sensor score [0, 1]
- `perm_pvalue`, `fdr`: Empirical p-value and FDR (computed separately within each sensor type)
- `confidence_tier`: Confidence classification (`Tier1_high`, `Tier2_medium`, `Tier3_exploratory`)

### Other Outputs
- `res.sender_scores`: `(metabolite, hmdb_id)` × cell type matrix of availability scores
- `res.receiver_scores`: Table of sensor scores per metabolite-sensor-cell type combination
- `res.availability_results`: Dictionary containing all intermediate calculation results (P/C/E matrices, pseudobulk, etc.)

## Notes
- **Transcriptomics-only**: CELL MESH estimates metabolite availability using expression proxies and prior knowledge. Direct metabolomics, spatial data, or perturbation experiments should be used to validate predictions.
- **Sensor scoring**: Uses robust min-max normalization of pseudobulk sensor gene expression
- **Communication score**: Geometric mean ensures both sender and receiver have meaningful scores
- **Sensor type stratification**: P-values and FDR are computed separately for each sensor type to avoid confounding
- **Permutation null**: Empirical p-values compare each observed full event key (`sender`, `receiver`, `metabolite`, `hmdb_id`, `sensor_gene`, `sensor_type`) against the same key after cell-type label permutation, with FDR stratified by sensor type.
