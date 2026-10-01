"""Independent condition filters and explicit pair reports; synthetic boundary tests only."""
import copy
import json
import re
import h5py
import pandas as pd
import pytest
from cellmesh import CellMeshResult
from cellmesh._pair_selection import pair_id
from cellmesh._report_pages import full_table
from cellmesh.result_io import _read
from test_reporting_and_persistence import make_result, SELECTOR

PAIR = ('HMDB0000001', 'R')

def payload(text, name):
    return json.loads(re.search(r'<script type="application/json" id="' + name + r'">(.*?)</script>', text, re.S)[1])

@pytest.fixture
def result():
    return make_result(n_perms=0)

def choose(result, **kwargs):
    return result.select_events(metabolite_sensor_pairs=[PAIR], source='all', **kwargs)

def test_lists_are_many_to_many_pairs_are_exact_and_original_values_unchanged():
    rows = [dict(sender='A', receiver='B', hmdb_id=m, sensor_gene=s, sensor_type='Receptor',
                 cell_mesh_score=.5, fdr_sensor_type=.03, perm_pvalue=.02)
            for m in ['HMDB0000001','HMDB0000002'] for s in ['R','S']]
    r = CellMeshResult(pd.DataFrame(rows), pd.DataFrame(), pd.DataFrame(), {})
    before = r.events.copy(deep=True)
    f = r.filter_events(hmdb_ids=['HMDB0000001','HMDB0000002'], sensor_genes=['R','S'])
    selected = r.select_events(metabolite_sensor_pairs=[('hmdb0000001','R'),('HMDB0000002','S')])
    assert len(f.events)==4 and len(selected.events)==2
    assert set(zip(selected.events.hmdb_id,selected.events.sensor_gene))==set(selected.pairs)
    assert selected.events.fdr_sensor_type.eq(.03).all()
    pd.testing.assert_frame_equal(r.events,before,check_exact=True)

def test_four_scopes_and_no_implicit_restore(result, tmp_path):
    f = result.filter_events(source='all', hmdb_ids=[PAIR[0]], sensor_genes=[PAIR[1]])
    selected = choose(result)
    before = copy.deepcopy(result.events)
    general = None
    for filtered, pairs in [(False,False),(True,False),(False,True),(True,True)]:
        files=result.save(tmp_path/f'{filtered}-{pairs}'/'analysis',
                          filtered_events=f if filtered else None,selected_events=selected if pairs else None)
        assert ('filtered_report' in files)==filtered
        assert ('filtered_csv' in files)==filtered
        assert ('selected_csv' in files)==pairs
        pages=[p for k,p in files.items() if k.startswith('pair_report:')]
        assert len(pages)==int(pairs)
        if filtered:
            data=payload(files['filtered_report'].read_text(),'interest-data')
            assert len(data['rows'])==len(f.events)
            assert all(bool(row['details'])==pairs for row in data['rows'])
        if pairs:
            text=pages[0].read_text()
            assert 'content="metabolite_sensor"' in text and 'id="pair-details"' in text
            assert len(payload(text,'pair-report-data')['event_keys'])==len(selected.events)>1
            assert len(re.findall(r'class="event-explanation"',text))==len(selected.events)
            assert len(pd.read_csv(files['selected_csv']))==len(selected.events)
        pd.testing.assert_frame_equal(result.events,before,check_exact=True)
        # Input-file locations differ by prefix; report science and global scope never change.
        if general is None: general=payload(files['report'].read_text(),'event-data')
        assert payload(files['report'].read_text(),'event-data')==general
    restored=CellMeshResult.load(files['result'])
    assert len(restored.get_selected_events())==1
    pd.testing.assert_frame_equal(restored.get_selected_events()[0].events,selected.events)
    fresh=restored.save(tmp_path/'fresh')
    assert 'filtered_csv' not in fresh and 'selected_csv' not in fresh

def test_append_preserves_general_science_other_pair_and_filtered_csv(result,tmp_path):
    files=result.save(tmp_path/'a')
    original={k:files[k].read_bytes() for k in ['report','positive_csv']}
    files=result.add_event_reports(choose(result,cell_pairs=[('B','A')]))
    pair_page=files['pair_report:'+pair_id(PAIR)]
    old_pair=pair_page.read_bytes()
    assert 'filtered_report' not in files
    assert '一级统计报告：沿用现有交付' in files['navigation'].read_text()
    assert '本次调用生成或更新 1 份' in files['navigation'].read_text()
    files=result.add_filtered_report(result.filter_events(source='all',sender_labels=['A']))
    assert pair_page.read_bytes()==old_pair
    csv=files['filtered_csv'].read_bytes()
    second=result.select_events(metabolite_sensor_pairs=[('HMDB0000002','R')],source='all')
    files=result.add_event_reports(second)
    assert pair_page.read_bytes()==old_pair and files['filtered_csv'].read_bytes()==csv
    assert all(files[k].read_bytes()==v for k,v in original.items())
    assert '原结果未生成' in files['pair_report:'+pair_id(second.pairs[0])].read_text()
    with pytest.raises(ValueError,match='replace=True'):result.add_event_reports(choose(result))
    result.add_event_reports(choose(result),replace=True)
    assert len(result.get_selected_events())==2
    with h5py.File(files['result']) as h:
        pd.testing.assert_frame_equal(_read(h['results/events']),result.events,check_exact=True)
        assert _read(h['metadata/delivery'])['report_layout']==4

@pytest.mark.parametrize('pairs',[['bad'],[PAIR,PAIR],[('HMDB0000001',)],'bad',[('HMDB9999999','R')]])
def test_bad_pair_requests_fail_before_writes(result,tmp_path,pairs):
    with pytest.raises((ValueError,TypeError)):
        result.save(tmp_path/'bad',selected_events=result.select_events(metabolite_sensor_pairs=pairs))
    assert not list(tmp_path.iterdir())

def test_mixed_scopes_and_legacy_calls_give_migration_errors(result,tmp_path):
    with pytest.raises(TypeError,match='filter_events'):
        result.select_events(hmdb_ids=[PAIR[0]],sensor_genes=[PAIR[1]])
    with pytest.raises(TypeError,match='filter_events'):
        result.save(tmp_path/'bad',filtered_events=choose(result))
    with pytest.raises(TypeError,match='select_events'):
        result.save(tmp_path/'bad',selected_events=result.filter_events(source='all'))
    with pytest.raises(TypeError,match='retired'):result.save(tmp_path/'old',events=[SELECTOR])
    assert not list(tmp_path.iterdir())

def test_empty_and_single_filter_remain_tier_two_and_pair_pages_not_capped(result,tmp_path):
    f=result.filter_events(source='all',hmdb_ids=[PAIR[0]],sensor_genes=[PAIR[1]],cell_pairs=[('A','B')])
    out=result.save(tmp_path/'filtered',filtered_events=f)
    assert 'filtered_report' in out and not any(k.startswith('pair_report:') for k in out)
    selected=result.select_events(metabolite_sensor_pairs=[PAIR,('HMDB0000002','R')],source='all')
    out=result.save(tmp_path/'selected',selected_events=selected,plot_top_n=1)
    assert len([k for k in out if k.startswith('pair_report:')])==2
    assert len(payload(out['pair_report:'+pair_id(PAIR)].read_text(),'pair-report-data')['event_keys'])==4
    empty=choose(result,sender_labels=[])
    out=result.save(tmp_path/'empty',selected_events=empty)
    assert '均未满足' in out['pair_report:'+pair_id(PAIR)].read_text()
    assert pd.read_csv(out['selected_csv']).empty

def test_standalone_csv_and_all_detail_tables_are_complete(result,tmp_path):
    selected=choose(result)
    files=selected.to_html(tmp_path/'pairs.html')
    assert files['report'].is_file() and 'positive_csv' not in files
    selected.to_csv(tmp_path/'manual.csv')
    assert len(pd.read_csv(tmp_path/'manual.csv'))==len(selected.events)
    frame=pd.DataFrame({'gene':[f'g{i}' for i in range(153)],'score':range(153)})
    text=full_table(frame,'genes','基因')
    assert len(payload(text,'detail-genes-data')['rows'])==153
    assert len(re.findall('data-event-index=',text))==50

def test_render_failure_rolls_back(result,tmp_path,monkeypatch):
    files=result.save(tmp_path/'a',filtered_events=result.filter_events(source='all'))
    before={p:p.read_bytes() for p in files.values()}; delivery=copy.deepcopy(result.delivery)
    import cellmesh._pair_reports as pages
    def fail(*args,**kwargs):raise RuntimeError('render failure')
    monkeypatch.setattr(pages,'render_pair',fail)
    with pytest.raises(RuntimeError,match='render failure'):result.add_event_reports(choose(result))
    assert result.delivery==delivery and all(p.read_bytes()==v for p,v in before.items())
    assert not list(tmp_path.rglob('pair_*.html'))

def test_report_false_still_exports_requested_csv_and_validates_pair_criterion(result,tmp_path):
    selected=choose(result)
    files=result.save(tmp_path/'a',report=False,selected_events=selected)
    assert set(files)=={'result','positive_csv','selected_csv'}
    with pytest.raises(ValueError,match='disabled'):result.add_event_reports(selected,replace=True)
    with pytest.raises(ValueError,match='Existing independent'):result.save(tmp_path/'a')
    scoped=result.select_events(metabolite_sensor_pairs=[PAIR],positive_threshold=.2)
    with pytest.raises(ValueError,match='different positive criterion'):
        result.save(tmp_path/'different',selected_events=scoped,positive_threshold=.05)
