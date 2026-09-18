"""Strict logical-record inspection before table parsers can infer structure."""
from __future__ import annotations

from contextlib import contextmanager
import csv
import sys
from threading import RLock

from pandas.io.common import get_handle


_CSV_LOCK = RLock()


def _record_location(start: int, end: int) -> str:
    return f"line {start}" if start == end else f"lines {start}-{end}"


def _reject_nul_text(path, *, encoding=None, compression="infer",
                     encoding_errors="strict", storage_options=None):
    """Reject NUL in decoded text, including records excluded by parser options."""
    # Scan BEFORE pandas: the C parser can silently turn G\0suffix into G, even
    # in a header probe. Scan decoded characters, not bytes (UTF-16/32 normally
    # contain zero bytes). Do not remove NUL or change the caller's parser.
    # Regression: tests/test_text_nul_inputs.py.
    with get_handle(path, "r", encoding=encoding, compression=compression,
                    errors=encoding_errors, storage_options=storage_options) as handles:
        for line_number, line in enumerate(handles.handle, start=1):
            if "\0" in line:
                raise ValueError(f"{path}, line {line_number}: NUL characters are not valid CSV/TSV text")


@contextmanager
def _delimited_records(path, *, delimiter=",", skip_blank_lines=True):
    """Yield (first physical line, last physical line, raw string fields).

    Use pandas' compression-aware handle so CSV validation supports the same
    encodings/archives as its final read. csv.reader preserves quoted commas,
    newlines and empty fields; splitting physical lines cannot validate CSV.
    """
    with get_handle(path, "r", encoding="utf-8-sig", compression="infer") as handles, _CSV_LOCK:
        # The csv module's small default field limit must not reject metadata
        # accepted by pandas. Restore it on all exits; serialize our readers
        # because this setting belongs to the csv module, not to one reader.
        limit = sys.maxsize
        while True:
            try:
                previous_limit = csv.field_size_limit(limit)
                break
            except OverflowError:
                # Some Python/platform combinations use a narrower C long.
                limit //= 10
        try:
            physical_lines = []

            def lines():
                for line in handles.handle:
                    physical_lines.append(line)
                    yield line

            reader = csv.reader(lines(), delimiter=delimiter, strict=True)

            def records():
                while True:
                    physical_lines.clear()
                    start = reader.line_num + 1
                    try:
                        row = next(reader)
                    except StopIteration:
                        return
                    except csv.Error as error:
                        location = _record_location(start, max(start, reader.line_num))
                        raise ValueError(f"{path}, {location}: invalid CSV/TSV syntax: {error}") from error
                    if skip_blank_lines and all(not line.strip() for line in physical_lines):
                        continue
                    if any("\0" in field for field in row):
                        raise ValueError(
                            f"{path}, {_record_location(start, reader.line_num)}: "
                            "NUL characters are not valid CSV/TSV text"
                        )
                    yield start, reader.line_num, row

            yield records()
        finally:
            csv.field_size_limit(previous_limit)
