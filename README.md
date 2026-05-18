# CellMESH
Metabolite-mediated Event Scoring with Sensor Hierarchies

**CELL MESH** = **Metabolite-mediated Event Scoring with Sensor Hierarchies**.

CELL MESH is a Python package for inferring metabolite-mediated cell-cell communication events from single-cell expression data, **fully based on the metabolite availability algorithm** since v0.4.0. It scores directed events of the form:

```text
sender cell type -> receiver cell type : metabolite -> sensor
```

## Algorithm Overview (v0.4.0+)

CELL MESH now uses a unified metabolite availability calculation to score sender cell types' ability to release each metabolite:

### 1. Enzyme Prior to Reaction Mapping
The enzyme-metabolite prior table is internally converted to reaction matrices using the following mapping:
| Role in prior | Direction | Matrix | Meaning |
|---|---|---|---|
| `production` | `product` | P | Metabolite production ability |
| `degradation` | `substrate` | C | Metabolite consumption ability |
| `usage` | `substrate` | C | Metabolite consumption ability |
| `export` | `exporter` | E | Metabolite efflux/transport ability |
| `import` | `exporter` | E | Metabolite efflux/transport ability |

### 2. Metabolite Availability Calculation
The final sender score (availability) is computed as:
```
availability = (P_norm + eps) * ((1 - C_norm + eps) ** beta) * (E_norm + eps)
```
Where:
- `P_norm`: Normalized production score [0, 1]
- `C_norm`: Normalized consumption score [0, 1]
- `E_norm`: Normalized efflux score [0, 1]
- Result range: Strictly between [0, 1], higher value indicates stronger ability to release the metabolite

### 3. Receiver Score Calculation
Receiver score is based solely on the expression and specificity of the sensor gene for the metabolite.

## Included Databases
The package includes built-in prior databases:
- `cell_mesh/data/enzyme_test.csv`: enzyme/reaction/metabolite table
- `cell_mesh/data/interaction_test.csv`: metabolite-receptor/sensor interaction table

These are automatically loaded if not explicitly provided.

## Install

```bash
pip install -e .
```

## Basic Usage

```python
from cell_mesh import run_cell_mesh

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
res.receiver_scores  # Receiver scores
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
| `beta_sensor` | `1.0` | Weight of sensor expression in receiver score |
| `beta_specificity` | `0.25` | Weight of sensor specificity in receiver score |
| `lower` | `5` | Lower percentile for robust min-max normalization |
| `upper` | `95` | Upper percentile for robust min-max normalization |
| `eps` | `0.05` | Small constant to avoid division by zero in availability calculation |
| `beta` | `0.5` | Exponent weight for consumption term in availability formula |
| `missing_C_norm` | `0.41` | Default C_norm value when no consumption evidence exists |
| `missing_E_norm` | `0.75` | Default E_norm value when no efflux evidence exists |
| `min_cells` | `1` | Minimum number of cells per cell type to be included |

### Removed Parameters (v0.4.0+)
The following parameters are no longer available:
- `reaction_table` (now internally generated from enzyme_metabolite)
- `use_new_availability` (now always enabled)
- `role_agg` (old algorithm removed)
- `min_cells_per_group` (replaced by min_cells)
- `alpha_prod`, `alpha_deg`, `alpha_export`, `alpha_specificity` (old algorithm removed)

## Inspect the Packaged Database

```python
from cell_mesh import load_cell_mesh_database

enzyme_metabolite, metabolite_sensor = load_cell_mesh_database()

print(enzyme_metabolite.head())
print(metabolite_sensor.head())
```

### Enzyme Table Columns
```
metabolite, hmdb_id, gene, role, weight, evidence_level, source, reaction
```

### Sensor Table Columns
```
metabolite, hmdb_id, sensor_gene, sensor_type, weight, evidence_level,
source, protein_name, reference
```

## Supported Sensor Types

| Input annotation | CELL MESH sensor type |
|---|---|
| Cell surface receptor | `surface_receptor` |
| Other receptor | `surface_receptor` |
| Transporter | `transporter` |
| Nuclear receptor | `nuclear_receptor` |
| Intracellular sensor | `intracellular_sensor` |

## Main Outputs

### `res.events`
Contains one row per communication event. Important columns:
- `sender`, `receiver`: Cell type pair
- `metabolite`, `hmdb_id`: Metabolite information
- `sensor_gene`, `sensor_type`: Sensor information
- `sender_score`: Metabolite availability in sender cell type [0, 1]
- `receiver_score`: Sensor activation score in receiver cell type [0, 1]
- `cell_mesh_score`: Combined event score = sender_score * receiver_score * prior_weight [0, 1]
- `perm_pvalue`, `fdr`: Empirical p-value and FDR (if n_perms > 0)
- `confidence_tier`: Confidence classification (`Tier1_high`, `Tier2_medium`, `Tier3_exploratory`)

### Other Outputs
- `res.sender_scores`: Metabolite × cell type matrix of availability scores
- `res.receiver_scores`: Table of receiver scores per metabolite-sensor-cell type combination
- `res.availability_results`: Dictionary containing all intermediate calculation results (P/C/E matrices, pseudobulk, etc.)
- `res.role_scores`: Empty dict (kept for backward compatibility, intermediate results now in availability_results)

## Notes
- **Transcriptomics-only**: CELL MESH estimates metabolite availability using expression proxies and prior knowledge. Direct metabolomics, spatial data, or perturbation experiments should be used to validate predictions.
- **Backward compatibility**: `run_cell_mesh` is the main entry point, the old `run_metcomm` alias is still available but deprecated.
- **Weight support**: The algorithm respects the `weight` column in prior tables, using weighted geometric mean for reaction scores.
- **Multi-gene support**: Genes in the `gene` column separated by `;`, `,` or `|` are parsed as a gene set for the same reaction.
