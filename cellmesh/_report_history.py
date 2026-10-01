"""Archive superseded filenames without changing their original contents."""
from __future__ import annotations

import hashlib
import json
import shutil
import zipfile


def stage_legacy_archive(staging, path, delivery, published):
    """Stage a lossless archive; the bundle publisher removes originals atomically."""
    from .result_io import _staged, _validate_targets
    if delivery.get('report_layout', 0) < 4 or delivery.get('selection') is None:
        return {}, []
    stem = path.name[:-len('.cellmesh.h5')]
    old_names = [('.selected.report.html', '.filtered.report.html'),
                 ('.selected_events.csv', '.filtered_events.csv')]
    originals, replacements = {}, {}
    for old_suffix, new_suffix in old_names:
        old = path.with_name(stem + old_suffix)
        checkpoint = old.parent / '.ipynb_checkpoints' / (old.stem + '-checkpoint' + old.suffix)
        candidates = [old, checkpoint]
        _validate_targets(candidates)
        candidates = [candidate for candidate in candidates if candidate.exists()]
        if not candidates:
            continue
        new = path.with_name(stem + new_suffix)
        _validate_targets([new])
        report_disabled = old_suffix.endswith('.html') and not delivery.get('report_enabled', True)
        if not report_disabled and new not in published.values() and not new.is_file():
            raise ValueError(f'Cannot archive {old.name}: current replacement {new.name} is missing')
        for candidate in candidates:
            originals[candidate] = candidate.read_bytes()
            replacements[candidate.relative_to(path.parent).as_posix()] = None if report_disabled else new.name
    if not originals:
        return {}, []
    records = {p.relative_to(path.parent).as_posix():
               dict(sha256=hashlib.sha256(data).hexdigest(), size=len(data),
                    replacement=replacements[p.relative_to(path.parent).as_posix()])
               for p, data in originals.items()}
    manifest = json.dumps(dict(reason='superseded selection filenames', files=records),
                          ensure_ascii=False, sort_keys=True, indent=2).encode('utf-8')
    identifier = hashlib.sha256(manifest).hexdigest()[:20]
    target = path.parent / (stem + '.reports') / 'history' / f'legacy-selection-{identifier}.zip'
    _validate_targets([target])
    members = {p.relative_to(path.parent).as_posix(): data for p, data in originals.items()}
    members['MANIFEST.json'] = manifest
    members['README.txt'] = (
        '历史归档：这些旧文件已被当前交付替代，不作为当前分析结果。\n'
        '原文件内容完整保留，摘要和当前文件名见 MANIFEST.json。\n'
        '旧 HTML 的相对链接基于原交付目录；本 ZIP 不包含整个历史交付及其依赖。\n'
        f'查看当前报告请使用 {stem}.reports/index.html。\n'
    ).encode('utf-8')
    staged = _staged(staging, target)
    staged.parent.mkdir(parents=True, exist_ok=True)
    if target.exists():
        try:
            with zipfile.ZipFile(target) as archive:
                matches = (sorted(archive.namelist()) == sorted(members) and
                           all(archive.read(name) == data for name, data in members.items()))
        except (OSError, zipfile.BadZipFile, KeyError) as error:
            raise ValueError(f'Existing history archive is invalid: {target}') from error
        if not matches:
            raise ValueError(f'Existing history archive differs: {target}')
        shutil.copy2(target, staged)
    else:
        with zipfile.ZipFile(staged, 'w', compression=zipfile.ZIP_DEFLATED) as archive:
            for name, data in members.items():
                archive.writestr(name, data)
    return {'history_archive': target}, list(originals)
