from __future__ import annotations

import re
from os import PathLike
from importlib.resources import files
from typing import Iterable, Tuple

import numpy as np
import pandas as pd

from .config import VALID_ROLES, VALID_SENSOR_TYPES
from .preprocess import _normalized_gene_names
from ._table_io import _delimited_records, _record_location


_ENZYME_ALIASES = {
    "metabolite": ("standard_metName",),
    "hmdb_id": ("HMDB_ID",),
    "gene": ("Gene_name",),
    "reaction": ("Reactions",),
}
_SENSOR_ALIASES = {
    "metabolite": ("standard_metName", "standard_metname"),
    "hmdb_id": ("HMDB_ID",),
    "sensor_gene": ("Gene_name", "gene_name", "gene"),
    "sensor_type": ("Annotation", "annotation"),
    "source": ("Database source", "database_source"),
    "protein_name": ("Protein_name",),
    "reference": ("Reference",),
}
_DIRECTION_ROLES = {
    "product": "production",
    "substrate": "degradation",
    "exporter": "export",
    # 原始酶库也使用 export 表示外排，不能因只识别 exporter 而丢弃证据。
    # 回归测试：tests/test_prior_input_consistency.py。
    "export": "export",
    # 保留旧评分接口的外排别名，并由两入口共用，避免再次出现解析差异。
    "transporter": "export",
}


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


def _read_prior_table(value, *, table_name):
    if isinstance(value, pd.DataFrame):
        return value
    # 先检查原始逻辑记录：pandas 可把超宽行的首项当作隐式索引，也会
    # 给短行补空值。不能用 index_col=False 或跳过坏行来掩盖错位/丢字段。
    # 引号内逗号、换行和合法空字段必须保留；不要为提速仅检查首行。
    # 同名表头检查仍不能替代后续同义列检查。回归：test_prior_csv_structure.py。
    with _delimited_records(value) as records:
        try:
            start, end, header = next(records)
        except StopIteration:
            raise ValueError(f"{table_name} file {value}: CSV header is missing") from None
        header_index = pd.Index(header)
        duplicates = header_index[header_index.duplicated()].tolist()
        if duplicates:
            raise ValueError(
                f"{table_name} has duplicate column names: {duplicates} "
                f"({value}, {_record_location(start, end)})"
            )
        for start, end, row in records:
            if len(row) != len(header):
                raise ValueError(
                    f"{table_name} file {value}, {_record_location(start, end)}: "
                    f"expected {len(header)} fields from header, got {len(row)}"
                )
    # CSV type inference must not turn reaction/reference IDs such as 001 into 1,
    # or interpret a literal name such as NA as missing. Blank fields stay missing.
    text_columns = {"role", "Direction", "direction", "evidence_level"}
    for aliases in (_ENZYME_ALIASES, _SENSOR_ALIASES):
        text_columns.update(aliases)
        for alternatives in aliases.values():
            text_columns.update(alternatives)
    result = pd.read_csv(
        value, encoding="utf-8-sig", dtype={column: str for column in text_columns},
        keep_default_na=False, na_values=[""],
    )
    if not isinstance(result.index, pd.RangeIndex):
        raise ValueError(f"{table_name} file {value}: unexpected implicit CSV index")
    return result


def _split_gene_field(value: object) -> list[tuple[str, str | None]]:
    """Split gene symbols on ; , or | outside optional [evidence] annotations."""
    if pd.isna(value):
        return []
    genes: list[tuple[str, str | None]] = []
    # Evidence text may itself contain separators, e.g. G1[reviewed, curated].
    # It must stay attached to its gene instead of becoming another symbol.
    for part in re.split(r"[;,|](?![^\[]*\])", str(value)):
        part = part.strip()
        if not part:
            continue
        m = re.fullmatch(r"([^\[\]]+)(?:\[([^\[\]]*)\])?", part)
        if m:
            gene = m.group(1).strip()
            if gene:
                genes.append((gene, m.group(2).strip() if m.group(2) else None))
        else:
            genes.append((part, None))
    return genes


def _valid_hmdb_mask(values: pd.Series) -> pd.Series:
    return _normalize_hmdb_series(values).notna()


def _normalize_hmdb_id(value: object):
    """Return one canonical HMDB identifier or ``np.nan`` when it is missing."""
    if pd.isna(value):
        return np.nan
    text = str(value).strip()
    if text.lower() in {"", "nan", "none", "null"}:
        return np.nan
    return text.upper()


def _normalize_hmdb_series(values: pd.Series) -> pd.Series:
    """Strip and case-normalize HMDB identifiers without stringifying missing data."""
    return values.map(_normalize_hmdb_id).astype(object)


def _prior_metabolite_names(*priors: pd.DataFrame) -> dict[str, str]:
    """Choose display names from normalized priors in order, then fall back to ID.

    Call with Enzyme before Interaction to prefer the first non-empty Enzyme
    name. Standalone scoring can supply just its own prior. This mapping is
    annotation only: never use its values as reaction, event or group keys.
    """
    identifiers = {}
    names = {}
    for prior in priors:
        for identifier, name in prior[["hmdb_id", "metabolite"]].itertuples(index=False, name=None):
            identifier = _normalize_hmdb_id(identifier)
            if pd.isna(identifier):
                continue
            identifiers[identifier] = None
            if identifier not in names and pd.notna(name) and str(name).strip():
                names[identifier] = str(name).strip()
    return {identifier: names.get(identifier, identifier) for identifier in identifiers}


def _enzyme_alias_key(target, column, value):
    """Return comparable field meaning, using None only for an empty value."""
    if target == "gene":
        return frozenset(gene for gene, _ in _split_gene_field(value)) or None
    if target == "hmdb_id":
        value = _normalize_hmdb_id(value)
    if pd.isna(value) or not str(value).strip():
        return None
    text = str(value).strip()
    if target == "role":
        text = text.lower()
        if column == "role":
            return ("role", text)
        if text in _DIRECTION_ROLES:
            return ("role", _DIRECTION_ROLES[text])
        # An unknown non-empty direction is not a blank eligible for fallback.
        return ("unknown_direction", text)
    return text


def _merge_gene_alias_evidence(values):
    """Combine equivalent gene lists without discarding inline provenance."""
    evidence_by_gene = {}
    for value in values:
        for gene, evidence in _split_gene_field(value):
            evidence_list = evidence_by_gene.setdefault(gene, [])
            if evidence and evidence not in evidence_list:
                evidence_list.append(evidence)
    return "; ".join(
        f"{gene}[{evidence}]" if evidence else gene
        for gene, annotations in evidence_by_gene.items()
        for evidence in (annotations or [None])
    )


def _canonical_enzyme_columns(frame):
    """Resolve and validate all enzyme aliases before expanding or filtering rows."""
    if frame.columns.has_duplicates:
        duplicates = frame.columns[frame.columns.duplicated()].tolist()
        raise ValueError(f"enzyme_metabolite has duplicate column names: {duplicates}")
    out = frame.copy()
    fields = [*_ENZYME_ALIASES.items(), ("role", ("Direction", "direction"))]
    for target, aliases in fields:
        columns = [name for name in (target, *aliases) if name in out]
        if not columns:
            description = "role or Direction/direction" if target == "role" else target
            raise ValueError(f"enzyme_metabolite is missing {description} columns")
        if len(columns) == 1 and target != "role":
            out = out.rename(columns={columns[0]: target})
            continue

        # Compare semantic values, not whichever spelling happens to come first.
        # In particular, export/exporter/transporter agree, and gene order and
        # separators do not change a reaction's gene set. Blanks can be filled;
        # conflicting non-empty aliases must fail before any records are lost.
        values = [out[column].astype(object).tolist() for column in columns]
        resolved = []
        conflicts = []
        for position, row_values in enumerate(zip(*values)):
            candidates = [
                (value, _enzyme_alias_key(target, column, value))
                for column, value in zip(columns, row_values)
            ]
            provided = [(value, key) for value, key in candidates if key is not None]
            if not provided:
                resolved.append(None)
                continue
            key = provided[0][1]
            if any(other_key != key for _, other_key in provided[1:]):
                conflicts.append(out.index[position])
            if target == "role":
                resolved.append(key[1] if key[0] == "role" else None)
            elif target == "gene":
                resolved.append(_merge_gene_alias_evidence([value for value, _ in provided]))
            else:
                resolved.append(provided[0][0])
        if conflicts:
            field = "role and direction" if target == "role" else f"{target} alias"
            raise ValueError(
                f"enzyme_metabolite has conflicting {field} values "
                f"in columns {columns} at rows {conflicts}"
            )
        # Drop consumed aliases now. After expansion, a retained Gene_name could
        # still contain G1;G2 beside gene=G1 and falsely conflict on the next load.
        if target == "role" and target not in out:
            # Raw direction tables historically append role after the existing
            # metadata columns. Preserve that output order for CSV round trips.
            out[target] = np.asarray(resolved, dtype=object)
            out = out.drop(columns=columns)
        else:
            out = out.rename(columns={columns[0]: target})
            out[target] = np.asarray(resolved, dtype=object)
            out = out.drop(columns=columns[1:])
    return out


def normalize_enzyme_database(enzyme_df: pd.DataFrame) -> pd.DataFrame:
    """Normalize raw direction tables and standard role tables alike.

    Input columns expected from the uploaded file:
    standard_metName, HMDB_ID, Reactions, Gene_name, Direction.
    Also accepts lower-case fields and the normalized ``role`` schema.
    Simultaneous aliases must agree after normalization. Empty values are filled
    from non-empty aliases; conflicting fields raise with row/column details.

    Output columns:
    metabolite, hmdb_id, gene, role, reaction, plus available metadata.
    Multi-gene fields in either schema are expanded before expression-gene
    matching. Reaction identity and existing metadata are preserved.
    """
    has_role = "role" in enzyme_df
    # 同义列先统一校验、补空值并移除旧列，不能静默择一或在丢行后校验。
    # 回归测试：tests/test_enzyme_aliases.py。
    out = _canonical_enzyme_columns(enzyme_df)
    out = out.loc[out["role"].isin(VALID_ROLES)].copy()
    out["hmdb_id"] = _normalize_hmdb_series(out["hmdb_id"])
    reaction = out["reaction"]
    invalid_reaction = reaction.isna() | reaction.astype(str).str.strip().eq("")
    if invalid_reaction.any():
        raise ValueError("enzyme_metabolite reaction values must be non-empty")
    out["reaction"] = reaction.astype(str).str.strip()

    # role 和 direction 表都必须先拆分基因，再与表达矩阵匹配。
    # 不能为提速跳过 role 表的拆分："G1;G2" 不是单个基因名，会被
    # gene.isin(var_names) 整行删除。保留 reaction，后续仍归并为一个反应。
    # 评分及置换只使用这里的单基因结果，不另设解析规则。
    # 回归测试：tests/test_enzyme_multigene.py。
    gene_lists = out["gene"].astype(object).map(_split_gene_field).tolist()
    positions = [i for i, entries in enumerate(gene_lists) for _ in entries]
    gene_evidence = [entry for entries in gene_lists for entry in entries]
    # Repeat existing rows so custom metadata keeps its dtype as well as value.
    out = out.iloc[positions].copy()
    out["gene"] = np.asarray([gene for gene, _ in gene_evidence], dtype=object)
    evidence = [value for _, value in gene_evidence]
    if "evidence_level" not in out and (not has_role or any(evidence)):
        out["evidence_level"] = np.asarray(
            [value or ("database" if not has_role else None) for value in evidence],
            dtype=object,
        )
    if not has_role and "source" not in out:
        out["source"] = "packaged_enzyme_test"
    # 标准化后只保留 role 表示反应角色，避免后续评分重新采用原始方向列，
    # 使大小写、空白或空值绕过这里的标准化，导致错误评分。
    # 回归测试：test_raw_direction_cannot_override_normalized_role。
    return out.drop(columns=["weight"], errors="ignore").drop_duplicates().reset_index(drop=True)


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


def _sensor_alias_key(target, value):
    """Compare field meaning, distinguishing missing types from Other receptor."""
    if target == "hmdb_id":
        value = _normalize_hmdb_id(value)
    if pd.isna(value) or not str(value).strip():
        return None
    text = str(value).strip()
    if target == "sensor_gene" and text.lower() == "nan":
        return None
    if target == "sensor_type":
        return _normalize_sensor_type(text)
    return text


def _canonical_sensor_columns(frame):
    """Validate and merge all sensor aliases before filtering or deduplication."""
    if frame.columns.has_duplicates:
        duplicates = frame.columns[frame.columns.duplicated()].tolist()
        raise ValueError(f"metabolite_sensor has duplicate column names: {duplicates}")
    out = frame.copy()
    required = {"metabolite", "hmdb_id", "sensor_gene", "sensor_type"}
    missing = [target for target in required
               if not any(name in out for name in (target, *_SENSOR_ALIASES[target]))]
    if missing:
        raise ValueError(f"metabolite_sensor is missing columns: {sorted(missing)}")
    for target, aliases in _SENSOR_ALIASES.items():
        columns = [name for name in (target, *aliases) if name in out]
        if not columns:
            continue
        if len(columns) == 1:
            out = out.rename(columns={columns[0]: target})
            continue

        # 不能为提速改回“标准列存在就忽略别名”：标准列可能为空，或与别名
        # 冲突。先逐行比较、补空，再过滤基因和去重，否则有效证据会静默丢失。
        # 类型先判断空值再归类，避免空值默认的 Other receptor 阻止正确补齐。
        # Interaction 基因按单个符号比较，不复用 Enzyme 的多基因拆分规则。
        # 回归测试：tests/test_interaction_aliases.py。
        values = [out[column].astype(object).tolist() for column in columns]
        resolved = []
        conflicts = []
        for position, row_values in enumerate(zip(*values)):
            provided = [(value, _sensor_alias_key(target, value)) for value in row_values]
            provided = [(value, key) for value, key in provided if key is not None]
            if not provided:
                resolved.append(None)
                continue
            key = provided[0][1]
            if any(other_key != key for _, other_key in provided[1:]):
                details = ", ".join(f"{column}={value!r}"
                                    for column, value in zip(columns, row_values))
                conflicts.append(f"{position + 1} (index={out.index[position]!r}): {details}")
            resolved.append(provided[0][0])
        if conflicts:
            raise ValueError(
                f"metabolite_sensor has conflicting {target} alias values in columns "
                f"{columns} at data rows (1-based): " + "; ".join(conflicts)
            )
        out = out.rename(columns={columns[0]: target})
        out[target] = np.asarray(resolved, dtype=object)
        out = out.drop(columns=columns[1:])
    return out


def _retain_sensor_annotation_evidence(out, original):
    """Keep raw type annotations as evidence without overwriting supplied evidence."""
    columns = [name for name in _SENSOR_ALIASES["sensor_type"] if name in original]
    if not columns:
        return
    evidence = []
    for row_values in zip(*(original[column].tolist() for column in columns)):
        annotations = {}
        for value in row_values:
            key = _sensor_alias_key("evidence_level", value)
            if key is not None:
                annotations.setdefault(key, value)
        values = list(annotations.values())
        if not values:
            evidence.append(None)
        elif len(values) == 1:
            evidence.append(values[0])
        else:
            evidence.append("; ".join(str(value) for value in values))
    if "evidence_level" not in out:
        out["evidence_level"] = np.asarray(evidence, dtype=object)
    else:
        missing = out["evidence_level"].map(
            lambda value: _sensor_alias_key("evidence_level", value) is None,
        )
        if missing.any():
            out["evidence_level"] = out["evidence_level"].astype(object).where(~missing, evidence)


def normalize_interaction_database(interaction_df: pd.DataFrame) -> pd.DataFrame:
    """Convert the packaged metabolite-receptor table into CELL MESH sensor prior format.

    Input columns expected from the uploaded file:
    ID, HMDB_ID, standard_metName, Gene_name, Protein_name, Annotation,
    Database source, Reference.
    Also accepts normalized/internal or lower-case schema:
    metabolite, hmdb_id, sensor_gene, sensor_type, evidence_level,
    source, protein_name, reference.

    Output columns:
    metabolite, hmdb_id, sensor_gene, sensor_type, evidence_level,
    source, protein_name, reference.

    Simultaneous aliases are compared by field meaning, merged when equivalent,
    and used to fill missing values. Conflicts and duplicate column names raise
    before any rows are filtered. Consumed aliases are removed. Raw Annotation/
    annotation text fills missing evidence_level values; explicit evidence wins.
    Same-type duplicate HMDB/sensor pairs keep their first occurrence; type
    conflicts are rejected before deduplication. Other metadata is preserved. Legacy/custom
    ``weight`` columns are intentionally ignored because the interaction
    database does not define a quantitative weight for sensor scoring.
    """
    out = _canonical_sensor_columns(interaction_df)
    _retain_sensor_annotation_evidence(out, interaction_df)
    if "sensor_type" not in interaction_df:
        for column in ("source", "protein_name", "reference"):
            if column not in out:
                out[column] = None
    genes = out["sensor_gene"].astype(str).str.strip()
    valid_genes = out["sensor_gene"].notna() & genes.ne("") & genes.str.lower().ne("nan")
    out["sensor_gene"] = genes
    out = out.loc[valid_genes].copy()
    out["sensor_type"] = out["sensor_type"].map(_normalize_sensor_type)
    out["hmdb_id"] = _normalize_hmdb_series(out["hmdb_id"])
    out = out.drop(columns=["weight"], errors="ignore")
    return _deduplicate_sensor_pairs(out).reset_index(drop=True)


def _deduplicate_sensor_pairs(sensor: pd.DataFrame) -> pd.DataFrame:
    # 必须在按 HMDB + 基因去重之前检查类型冲突，否则 CSV 会先丢掉
    # 冲突行，使后续校验失效，并由文件行顺序决定类型及分类型 FDR。
    # 标准化与运行时校验共用此规则；同类型重复只保留首条元数据。
    # 回归测试：tests/test_prior_input_consistency.py。
    sensor_key = ["hmdb_id", "sensor_gene"]
    # Missing identifiers are filtered by validate_priors, not valid sensor keys.
    type_counts = sensor.groupby(sensor_key)["sensor_type"].nunique()
    conflicting_keys = type_counts[type_counts > 1]
    if not conflicting_keys.empty:
        conflict_details = []
        for hmdb_id, sensor_gene in conflicting_keys.index:
            mask = sensor["hmdb_id"].eq(hmdb_id) & sensor["sensor_gene"].eq(sensor_gene)
            sensor_types = sorted(sensor.loc[mask, "sensor_type"].unique())
            conflict_details.append(f"({hmdb_id}, {sensor_gene}): {sensor_types}")
        raise ValueError(
            "metabolite_sensor has conflicting sensor_type values for the same "
            "canonical hmdb_id and sensor_gene: " + "; ".join(conflict_details)
        )
    return sensor.drop_duplicates(subset=sensor_key, keep="first")


def load_cell_mesh_database(
    enzyme_file: str | PathLike[str] | pd.DataFrame | None = None,
    interaction_file: str | PathLike[str] | pd.DataFrame | None = None,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Load packaged or user-provided CELL MESH prior databases.

    Omitted paths use the independently highest numeric versions of
    ``Enzyme<version>.csv`` and ``Interaction<version>.csv`` in ``cellmesh/data``.
    Selection is by version number, not file modification time. Each database
    falls back to its legacy test CSV when no versioned file is available.
    Paths may be strings or path-like objects. CSVs and DataFrames use the same
    schema normalization and sensor conflict checks. DataFrames are copied;
    each prior can be supplied independently in raw or normalized form.

    Returns
    -------
    enzyme_metabolite, metabolite_sensor
        Two normalized prior tables ready for :func:`cellmesh.run_cell_mesh`.
    """
    for name, value in (("enzyme_file", enzyme_file), ("interaction_file", interaction_file)):
        if value is not None and not isinstance(value, (str, PathLike, pd.DataFrame)):
            raise TypeError(f"{name} must be None, a CSV path, or a pandas DataFrame")

    if enzyme_file is None or interaction_file is None:
        default_enzyme_path, default_interaction_path = _default_database_paths()
        enzyme_file = default_enzyme_path if enzyme_file is None else enzyme_file
        interaction_file = default_interaction_path if interaction_file is None else interaction_file

    enzyme = normalize_enzyme_database(_read_prior_table(enzyme_file, table_name="enzyme_metabolite"))
    sensor = normalize_interaction_database(_read_prior_table(interaction_file, table_name="metabolite_sensor"))
    return enzyme, sensor


def load_default_priors() -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Alias for :func:`load_cell_mesh_database`."""
    return load_cell_mesh_database()


def validate_priors(
    enzyme_metabolite: pd.DataFrame,
    metabolite_sensor: pd.DataFrame,
    var_names,
    *,
    filter_enzyme_genes: bool = True,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Runtime prior cleaning for run_cell_mesh().

    This validates the normalized prior schemas, filters invalid HMDB IDs,
    restricts roles/sensor types to supported values, removes unsupported legacy
    weight columns, and guarantees one sensor record per canonical HMDB ID and
    sensor gene. Sensors are restricted to measured genes. The default also
    filters enzyme genes for compatibility with standalone callers;
    run_cell_mesh uses filter_enzyme_genes=False to preserve complete reaction
    identities until reaction grouping and deduplication have finished.
    """
    genes = set(_normalized_gene_names(var_names))

    # Direct callers must get the same gene expansion as database loading;
    # filtering composite strings first silently discards valid reaction genes.
    enz = normalize_enzyme_database(enzyme_metabolite)
    enz = enz[_valid_hmdb_mask(enz["hmdb_id"])]
    if filter_enzyme_genes:
        enz = enz[enz["gene"].isin(genes)]

    # Enzyme priors are unweighted. Drop legacy/custom weight columns so they
    # cannot imply a scoring effect that the packaged database does not define.
    enz = enz.drop(columns=["weight"], errors="ignore")

    sen = metabolite_sensor.copy()
    required_sen = {"metabolite", "hmdb_id", "sensor_gene", "sensor_type"}
    missing = required_sen - set(sen.columns)
    if missing:
        raise ValueError(f"metabolite_sensor is missing columns: {sorted(missing)}")

    sen["sensor_gene"] = sen["sensor_gene"].astype(str)
    sen["hmdb_id"] = _normalize_hmdb_series(sen["hmdb_id"])
    sen["sensor_type"] = sen["sensor_type"].astype(str)
    sen = sen[_valid_hmdb_mask(sen["hmdb_id"])]
    sen = sen[sen["sensor_type"].isin(VALID_SENSOR_TYPES)]
    sen = sen[sen["sensor_gene"].isin(genes)]

    # Sensor priors are unweighted.  In particular, do not silently turn a
    # user-supplied evidence field into a numeric receiver-score multiplier.
    sen = sen.drop(columns=["weight"], errors="ignore")

    sen = _deduplicate_sensor_pairs(sen)

    return enz.reset_index(drop=True), sen.reset_index(drop=True)
