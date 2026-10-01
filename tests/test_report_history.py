"""Automatic migration/transaction boundaries; not a real-data acceptance case."""
import copy
import json
from pathlib import Path
import zipfile

import h5py
import pandas as pd
import pytest

import cellmesh.result_io as io
from test_reporting_and_persistence import make_result


def old_files(prefix):
    files = {prefix.with_name(prefix.name + '.selected.report.html'): b'<html>Old Dot area</html>',
             prefix.with_name(prefix.name + '.selected_events.csv'): b'event,score\ne1,0.5\n'}
    for path, data in files.items():
        path.write_bytes(data)
    return files


def assert_archive(root, originals):
    archives = list(root.rglob('legacy-selection-*.zip'))
    assert len(archives) == 1
    with zipfile.ZipFile(archives[0]) as z:
        record = json.loads(z.read('MANIFEST.json'))
        for path, data in originals.items():
            assert not path.exists()
            name = path.relative_to(root).as_posix()
            assert z.read(name) == data
            assert record['files'][name]['size'] == len(data)
    return archives[0]


@pytest.mark.parametrize('operation,layout', [('save', 3), ('filtered', 3), ('pairs', 4)])
def test_save_append_archive_old_names(operation, layout, tmp_path):
    result = make_result(n_perms=0)
    filtered = result.filter_events(source='all')
    files = result.save(tmp_path/'analysis', filtered_events=filtered)
    science = result.events.copy(deep=True)
    parameters = copy.deepcopy(result.parameters)
    general = files['report'].read_bytes()
    originals = old_files(tmp_path/'analysis')
    if layout == 3:
        # Simulate the old naming/metadata contract, not a scientific data change.
        files['filtered_report'].unlink()
        files['filtered_csv'].unlink()
        result.delivery['report_layout'] = 3
        with h5py.File(files['result'], 'r+') as h:
            del h['metadata/delivery']
            io._write(h['metadata'], 'delivery', result.delivery, include_null_scores=False)
    unrelated = tmp_path/'my.selected.report.html'
    unrelated.write_text('user document')
    if operation == 'save':
        files = result.save(tmp_path/'analysis', filtered_events=filtered)
    elif operation == 'filtered':
        files = result.add_filtered_report(filtered, replace=True)
    else:
        selected = result.select_events(metabolite_sensor_pairs=[('HMDB0000001', 'R')], source='all')
        files = result.add_event_reports(selected)
    assert_archive(tmp_path, originals)
    assert files['filtered_report'].is_file() and files['filtered_csv'].is_file()
    assert 'analysis.selected.report.html' not in files['navigation'].read_text()
    assert unrelated.read_text() == 'user document'
    if operation != 'save':
        assert files['report'].read_bytes() == general
    restored = io.load_result(files['result'])
    pd.testing.assert_frame_equal(restored.events, science, check_exact=True)
    assert io._pack(restored.parameters) == io._pack(parameters)


@pytest.fixture
def migrated(tmp_path):
    """Small file-management fixture; no scientific calculation or visual review."""
    prefix = tmp_path/'analysis'
    with h5py.File(tmp_path/'analysis.cellmesh.h5', 'w') as h:
        h.attrs.update(format='CELL MESH result', schema_version=io.FORMAT_VERSION)
        io._write(h.create_group('metadata'), 'delivery',
                  dict(report_layout=4, selection={'keys': []}, report_enabled=True), include_null_scores=False)
    (tmp_path/'analysis.filtered.report.html').write_text('current HTML')
    (tmp_path/'analysis.filtered_events.csv').write_text('current CSV')
    originals = old_files(prefix)
    return prefix, originals


def test_repair_is_lossless_idempotent_and_preserves_current_files(migrated):
    prefix, originals = migrated
    before = {p: p.read_bytes() for p in prefix.parent.iterdir() if p not in originals}
    files = io.archive_legacy_reports(prefix)
    archive = assert_archive(prefix.parent, originals)
    assert files == {'history_archive': archive}
    assert all(p.read_bytes() == data for p, data in before.items())
    archived = archive.read_bytes()
    assert io.archive_legacy_reports(prefix) == {}
    # A copied legacy file reappearing uses the same verified archive.
    for p, data in originals.items():
        p.write_bytes(data)
    io.archive_legacy_reports(prefix)
    assert_archive(prefix.parent, originals)
    assert archive.read_bytes() == archived


def test_incomplete_replacement_and_legacy_layout_remain_untouched(migrated):
    prefix, originals = migrated
    current = prefix.with_name('analysis.filtered.report.html')
    current.unlink()
    with pytest.raises(ValueError, match='replacement.*missing'):
        io.archive_legacy_reports(prefix)
    assert all(p.read_bytes() == data for p, data in originals.items())
    with h5py.File(prefix.with_name('analysis.cellmesh.h5'), 'r+') as h:
        del h['metadata/delivery']
        io._write(h['metadata'], 'delivery', dict(report_layout=3), include_null_scores=False)
    with pytest.raises(ValueError, match='Migrate reports'):
        io.archive_legacy_reports(prefix)
    assert not list(prefix.parent.rglob('*.zip'))


@pytest.mark.parametrize('failure', ['archive', 'remove_second_original'])
def test_archive_transaction_rolls_back(migrated, monkeypatch, failure):
    prefix, originals = migrated
    before = {p: p.read_bytes() for p in prefix.parent.iterdir()}
    original_replace, original_unlink = io.os.replace, Path.unlink
    failed = False
    def replace(src, dst):
        nonlocal failed
        if failure == 'archive' and str(dst).endswith('.zip') and not failed:
            failed = True
            raise OSError('injected archive failure')
        return original_replace(src, dst)
    def unlink(path, *args, **kwargs):
        nonlocal failed
        if failure == 'remove_second_original' and path.name == 'analysis.selected_events.csv' and not failed:
            failed = True
            raise OSError('injected archive failure')
        return original_unlink(path, *args, **kwargs)
    with monkeypatch.context() as m:
        m.setattr(io.os, 'replace', replace)
        m.setattr(Path, 'unlink', unlink)
        with pytest.raises(OSError, match='injected'):
            io.archive_legacy_reports(prefix)
    assert all(p.read_bytes() == data for p, data in before.items())
    assert not list(prefix.parent.rglob('*.zip'))
    assert not list(prefix.parent.glob('.cellmesh-*'))


def test_archive_conflict_preserves_originals(migrated):
    prefix, originals = migrated
    io.archive_legacy_reports(prefix)
    archive = assert_archive(prefix.parent, originals)
    archive.write_bytes(b'not the original archive')
    for p, data in originals.items():
        p.write_bytes(data)
    with pytest.raises(ValueError, match='history archive is invalid'):
        io.archive_legacy_reports(prefix)
    assert all(p.read_bytes() == data for p, data in originals.items())
    assert archive.read_bytes() == b'not the original archive'


@pytest.mark.parametrize('operation', ['save', 'append'])
def test_bundle_updates_roll_back_when_archiving_fails(tmp_path, monkeypatch, operation):
    result = make_result(n_perms=0)
    filtered = result.filter_events(source='all')
    result.save(tmp_path/'analysis', filtered_events=filtered)
    old_files(tmp_path/'analysis')
    before = {p: p.read_bytes() for p in tmp_path.rglob('*') if p.is_file()}
    delivery = copy.deepcopy(result.delivery)
    replace = io.os.replace
    def fail(src, dst):
        if str(dst).endswith('.zip'):
            raise OSError('archive publication failed')
        return replace(src, dst)
    with monkeypatch.context() as m:
        m.setattr(io.os, 'replace', fail)
        with pytest.raises(OSError, match='archive publication failed'):
            if operation == 'save':
                result.save(tmp_path/'analysis', filtered_events=filtered)
            else:
                result.add_filtered_report(filtered, replace=True)
    assert result.delivery == delivery
    assert {p: p.read_bytes() for p in tmp_path.rglob('*') if p.is_file()} == before


def test_symlink_is_not_archived(migrated):
    prefix, originals = migrated
    old = next(iter(originals))
    old.unlink()
    other = prefix.with_name('user.html')
    other.write_text('user content')
    old.symlink_to(other)
    with pytest.raises(ValueError, match='symlink'):
        io.archive_legacy_reports(prefix)
    assert old.is_symlink() and other.read_text() == 'user content'
    assert not list(prefix.parent.rglob('*.zip'))


def test_only_superseded_checkpoints_are_archived(migrated):
    prefix, originals = migrated
    checkpoints = prefix.parent/'.ipynb_checkpoints'
    checkpoints.mkdir()
    old = checkpoints/'analysis.selected.report-checkpoint.html'
    old.write_bytes(b'old notebook checkpoint')
    originals[old] = old.read_bytes()
    current = checkpoints/'analysis.filtered.report-checkpoint.html'
    current.write_bytes(b'current checkpoint')
    io.archive_legacy_reports(prefix)
    assert_archive(prefix.parent, originals)
    assert current.read_bytes() == b'current checkpoint'
