from __future__ import annotations

import re
from importlib.resources import files
from typing import Iterable, Tuple

import numpy as np
import pandas as pd

from .config import VALID_ROLES, VALID_SENSOR_TYPES


def _data_path(filename: str):
    return files("cellmesh.data").joinpath(filename)


def _version_tuple(version: str) -> tuple[int, ...]:
    return tuple(int(part) for part in version.split("."))


def _find_versioned_database_files(entries: Iterable[object], prefix: str) -> dict[tuple[int, ...], object]:
    pattern = re.compile(rf"^{re.escape(prefix)}(\d+(?:\.\d+)*)\.csv$")
    matches = {}
    for entry in entries:
        name = getattr(entry, "name", str(entry))
        match = pattern.match(name)
        if match:
            matches[_version_tuple(match.group(1))] = entry
    return matches


def _select_default_database_paths(entries: Iterable[object], enzyme_fallback, interaction_fallback):
    entries = list(entries)
    enzyme_files = _find_versioned_database_files(entries, "Enzyme")
    interaction_files = _find_versioned_database_files(entries, "Interaction")

    enzyme_path = enzyme_files[max(enzyme_files)] if enzyme_files else enzyme_fallback
    interaction_path = interaction_files[max(interaction_files)] if interaction_files else interaction_fallback
    return enzyme_path, interaction_path


def _default_database_paths():
    data_dir = files("cellmesh.data")
    return _select_default_database_paths(
        data_dir.iterdir(),
        _data_path("enzyme_test.csv"),
        _data_path("interaction_test.csv"),
    )


def _first_present(row: pd.Series, *columns: str):
    for column in columns:
        if column in row.index:
            return row.get(column)
    return None


def _split_gene_field(value: object) -> list[tuple[str, str | None]]:
    """Parse fields such as 'PTGR2[Enzyme]; PTGR1[Enzyme]'."""
    if pd.isna(value):
        return []
    genes: list[tuple[str, str | None]] = []
    for part in str(value).split(";"):
        part = part.strip()
        if not part:
            continue
        m = re.match(r"^([^\[]+)(?:\[([^\]]+)\])?$", part)
        if m:
            genes.append((m.group(1).strip(), m.group(2).strip() if m.group(2) else None))
        else:
            genes.append((part, None))
    return genes


def _valid_hmdb_mask(values: pd.Series) -> pd.Series:
    text = values.astype(str).str.strip().str.lower()
    return pd.notna(values) & ~text.isin({"", "nan", "none", "null"})


def normalize_enzyme_database(enzyme_df: pd.DataFrame) -> pd.DataFrame:
    """Convert the packaged enzyme table into CELL MESH enzyme prior format.

    Input columns expected from the uploaded file:
    standard_metName, HMDB_ID, Reactions, Gene_name, Direction.
    Also accepts the newer lower-case schema:
    metabolite, HMDB_ID/hmdb_id, reaction, gene, direction.

    Output columns:
    metabolite, hmdb_id, gene, role, weight, evidence_level, source, reaction.
    """
    role_map = {
        "product": "production",
        "substrate": "degradation",
        "exporter": "export",
    }
    rows = []
    for _, row in enzyme_df.iterrows():
        direction = str(row.get("Direction", row.get("direction", ""))).strip().lower()
        role = role_map.get(direction)
        if role is None:
            continue
        for gene, evidence in _split_gene_field(row.get("Gene_name", row.get("gene"))):
            if not gene:
                continue
            rows.append(
                {
                    "metabolite": row.get("standard_metName", row.get("metabolite")),
                    "hmdb_id": row.get("HMDB_ID", row.get("hmdb_id")),
                    "gene": gene,
                    "role": role,
                    "weight": 1.0,
                    "evidence_level": evidence or "database",
                    "source": "packaged_enzyme_test",
                    "reaction": row.get("Reactions", row.get("reaction")),
                }
            )
    out = pd.DataFrame(rows)
    if out.empty:
        return pd.DataFrame(columns=["metabolite", "hmdb_id", "gene", "role", "weight", "evidence_level", "source", "reaction"])
    return out.drop_duplicates().reset_index(drop=True)


def _normalize_sensor_type(annotation: object) -> str:
    """Normalize sensor type to the three required categories.
    
    Maps to: "Cell surface receptor", "Transporter", "Other receptor"
    """
    if pd.isna(annotation):
        return "Other receptor"
    text = str(annotation).strip()
    text_lower = text.lower()
    
    # Exact matches first
    if text in ["Cell surface receptor", "Transporter", "Other receptor"]:
        return text
    
    # Case-insensitive matching
    if "cell surface" in text_lower or "surface receptor" in text_lower:
        return "Cell surface receptor"
    if "transport" in text_lower:
        return "Transporter"
    
    # Everything else goes to Other receptor
    return "Other receptor"


def normalize_interaction_database(interaction_df: pd.DataFrame) -> pd.DataFrame:
    """Convert the packaged metabolite-receptor table into CELL MESH sensor prior format.

    Input columns expected from the uploaded file:
    ID, HMDB_ID, standard_metName, Gene_name, Protein_name, Annotation,
    Database source, Reference.
    Also accepts normalized/internal or lower-case schema:
    metabolite, hmdb_id, sensor_gene, sensor_type, weight, evidence_level,
    source, protein_name, reference.

    Output columns:
    metabolite, hmdb_id, sensor_gene, sensor_type, weight, evidence_level,
    source, protein_name, reference.
    """
    rows = []
    for _, row in interaction_df.iterrows():
        gene = str(_first_present(row, "Gene_name", "gene_name", "sensor_gene", "gene") or "").strip()
        if not gene or gene.lower() == "nan":
            continue
        sensor_type = _first_present(row, "Annotation", "annotation", "sensor_type")
        rows.append(
            {
                "metabolite": _first_present(row, "standard_metName", "standard_metname", "metabolite"),
                "hmdb_id": _first_present(row, "HMDB_ID", "hmdb_id"),
                "sensor_gene": gene,
                "sensor_type": _normalize_sensor_type(sensor_type),
                "weight": pd.to_numeric(_first_present(row, "weight"), errors="coerce") if "weight" in row.index else 1.0,
                "evidence_level": _first_present(row, "evidence_level", "Annotation", "annotation"),
                "source": _first_present(row, "source", "Database source", "database_source"),
                "protein_name": _first_present(row, "protein_name", "Protein_name"),
                "reference": _first_present(row, "reference", "Reference"),
            }
        )
    out = pd.DataFrame(rows)
    if out.empty:
        return pd.DataFrame(columns=["metabolite", "hmdb_id", "sensor_gene", "sensor_type", "weight", "evidence_level", "source", "protein_name", "reference"])
    out["weight"] = pd.to_numeric(out["weight"], errors="coerce").fillna(1.0)
    return (
        out.sort_values("weight", ascending=False)
        .drop_duplicates(subset=["hmdb_id", "sensor_gene"])
        .reset_index(drop=True)
    )


def load_cell_mesh_database(
    enzyme_file: str | None = None,
    interaction_file: str | None = None,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Load packaged or user-provided CELL MESH prior databases.

    Returns
    -------
    enzyme_metabolite, metabolite_sensor
        Two normalized prior tables ready for :func:`cellmesh.run_cell_mesh`.
    """
    default_enzyme_path, default_interaction_path = _default_database_paths()
    enzyme_path = enzyme_file if enzyme_file is not None else default_enzyme_path
    interaction_path = interaction_file if interaction_file is not None else default_interaction_path

    enzyme_raw = pd.read_csv(enzyme_path, encoding="utf-8-sig")
    interaction_raw = pd.read_csv(interaction_path, encoding="utf-8-sig")
    return normalize_enzyme_database(enzyme_raw), normalize_interaction_database(interaction_raw)


def load_default_priors() -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Alias for :func:`load_cell_mesh_database`."""
    return load_cell_mesh_database()


def validate_priors(
    enzyme_metabolite: pd.DataFrame,
    metabolite_sensor: pd.DataFrame,
    var_names,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Runtime prior cleaning for run_cell_mesh().

    This validates the normalized prior schemas, filters invalid HMDB IDs,
    restricts roles/sensor types to supported values, keeps only genes present
    in the expression matrix, and normalizes weights.
    """
    genes = set(pd.Index(var_names).astype(str))

    enz = enzyme_metabolite.copy()
    required_enz = {"metabolite", "hmdb_id", "gene", "role"}
    missing = required_enz - set(enz.columns)
    if missing:
        raise ValueError(f"enzyme_metabolite is missing columns: {sorted(missing)}")

    enz["gene"] = enz["gene"].astype(str)
    enz["hmdb_id"] = enz["hmdb_id"].astype(object).where(pd.notna(enz["hmdb_id"]), np.nan)
    enz["role"] = enz["role"].astype(str).str.lower()
    enz = enz[_valid_hmdb_mask(enz["hmdb_id"])]
    enz = enz[enz["role"].isin(VALID_ROLES)]
    enz = enz[enz["gene"].isin(genes)]

    if "weight" not in enz:
        enz["weight"] = 1.0
    enz["weight"] = pd.to_numeric(enz["weight"], errors="coerce").fillna(1.0)

    sen = metabolite_sensor.copy()
    required_sen = {"metabolite", "hmdb_id", "sensor_gene", "sensor_type"}
    missing = required_sen - set(sen.columns)
    if missing:
        raise ValueError(f"metabolite_sensor is missing columns: {sorted(missing)}")

    sen["sensor_gene"] = sen["sensor_gene"].astype(str)
    sen["hmdb_id"] = sen["hmdb_id"].astype(object).where(pd.notna(sen["hmdb_id"]), np.nan)
    sen["sensor_type"] = sen["sensor_type"].astype(str)
    sen = sen[_valid_hmdb_mask(sen["hmdb_id"])]
    sen = sen[sen["sensor_type"].isin(VALID_SENSOR_TYPES)]
    sen = sen[sen["sensor_gene"].isin(genes)]

    if "weight" not in sen:
        sen["weight"] = 1.0
    sen["weight"] = pd.to_numeric(sen["weight"], errors="coerce").fillna(1.0)

    return enz.reset_index(drop=True), sen.reset_index(drop=True)
