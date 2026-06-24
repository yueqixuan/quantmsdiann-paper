"""Regression tests for the paper_numbers collectors in scripts/rebuild.py.

Locks the producer/collector schema contract that previously drifted: the
single-cell producer writes mv_per_cell.tsv with a 'pg_count' column, but the
collector read a non-existent 'n_proteins' column and silently dropped the
per-cell median protein-group numbers from the manuscript audit.
"""
import pathlib
import sys

import pandas as pd
import pytest

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
import scripts.rebuild as R  # noqa: E402


def _notes(numbers):
    """map note -> value for the rows a collector added."""
    return {note: value for _key, value, _source, note in numbers.rows}


def test_collect_single_cell_reads_pg_count(monkeypatch):
    """mv_per_cell.tsv carries 'pg_count'; the collector must emit the per-cell
    median from it (the bug read 'n_proteins' and emitted nothing)."""
    def fake_read(path):
        if path.name == 'sc_totals.tsv':
            return pd.DataFrame({'dataset': ['HeLa'], 'version': ['v1'],
                                 'precursors': [10], 'proteins': [5]})
        if path.name == 'mv_per_cell.tsv':
            return pd.DataFrame({'dataset': ['HeLa', 'HeLa'], 'version': ['v1', 'v1'],
                                 'pg_count': [100, 200]})
        return None

    monkeypatch.setattr(R, '_read', fake_read)
    numbers = R.Numbers()
    R.collect_single_cell(numbers)
    notes = _notes(numbers)
    assert notes.get('HeLa v1 median protein groups/cell') == '150'


def test_collect_single_cell_fails_closed_on_schema_drift(monkeypatch):
    """A present-but-wrong-schema mv_per_cell.tsv must raise, not silently skip,
    so the audit can never claim success while dropping a user-visible metric."""
    def fake_read(path):
        if path.name == 'mv_per_cell.tsv':
            return pd.DataFrame({'dataset': ['HeLa'], 'version': ['v1'],
                                 'n_proteins': [100]})  # wrong column
        return None

    monkeypatch.setattr(R, '_read', fake_read)
    with pytest.raises(ValueError, match='schema drift'):
        R.collect_single_cell(R.Numbers())


def test_single_cell_table_builder_declares_pg_count():
    """Guard the other side of the contract: the producer's declared columns for
    mv_per_cell.tsv include 'pg_count' (matches what the collector reads)."""
    import inspect
    src = inspect.getsource(R)
    assert "'mv_per_cell.tsv': pd.DataFrame(per_cell, columns=['dataset', 'version', 'pg_count']" in src
