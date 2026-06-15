from __future__ import annotations

import re
from importlib.resources import files
from typing import Tuple

import pandas as pd


def _data_path(filename: str):
    return files("cellmesh.data").joinpath(filename)


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

    Output columns:
    metabolite, hmdb_id, sensor_gene, sensor_type, weight, evidence_level,
    source, protein_name, reference.
    """
    rows = []
    for _, row in interaction_df.iterrows():
        gene = str(row.get("Gene_name", "")).strip()
        if not gene or gene.lower() == "nan":
            continue
        rows.append(
            {
                "metabolite": row.get("standard_metName"),
                "hmdb_id": row.get("HMDB_ID"),
                "sensor_gene": gene,
                "sensor_type": _normalize_sensor_type(row.get("Annotation")),
                "weight": 1.0,
                "evidence_level": row.get("Annotation"),
                "source": row.get("Database source"),
                "protein_name": row.get("Protein_name"),
                "reference": row.get("Reference"),
            }
        )
    out = pd.DataFrame(rows)
    if out.empty:
        return pd.DataFrame(columns=["metabolite", "hmdb_id", "sensor_gene", "sensor_type", "weight", "evidence_level", "source", "protein_name", "reference"])
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
        Two normalized prior tables ready for :func:`cell_mesh.run_cell_mesh`.
    """
    enzyme_path = enzyme_file if enzyme_file is not None else _data_path("enzyme_test.csv")
    interaction_path = interaction_file if interaction_file is not None else _data_path("interaction_test.csv")

    enzyme_raw = pd.read_csv(enzyme_path, encoding="utf-8-sig")
    interaction_raw = pd.read_csv(interaction_path, encoding="utf-8-sig")
    return normalize_enzyme_database(enzyme_raw), normalize_interaction_database(interaction_raw)


def load_default_priors() -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Alias for :func:`load_cell_mesh_database`."""
    return load_cell_mesh_database()
