"""Automated boundary tests; real-data human acceptance lives in the sibling notebook."""
import copy
import json
import re

import h5py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest

from cellmesh import CellMeshResult, plot_communication_dotmap
from cellmesh._communication_dotmap import communication_scope, matrix_axes, dotmap_encoding, dot_color, dot_area
from cellmesh._report_dotmap import communication_dotmap_html
from cellmesh._report_plots import selection_plots_html
from cellmesh._selection_evidence import START, END
from cellmesh.result_io import _read, _write
from test_input_archive import inputs, run
from test_reporting_and_persistence import make_result, SELECTOR


def payload(text, key):
    return json.loads(re.search(r'id="' + key + r'">(.*?)</script>', text, re.S)[1])


def block(path):
    return path.read_text().split(START)[1].split(END)[0]


def path_frame():
    return pd.DataFrame([dict(sender=f'S{i % 3}', receiver=f'R{i}', hmdb_id=f'HMDB{1+i%2:07d}',
        sensor_gene='R<script>unsafe()</script>', metabolite='Same name', cell_mesh_score=.4,
        fdr_sensor_type=.9, passes_min_cells=False) for i in range(105)])


def test_one_matrix_preserves_complete_events_and_composite_axis_identities():
    frame = path_frame()
    frame.loc[0, 'cell_mesh_score'] = np.nan
    before = frame.copy(deep=True)
    scope = communication_scope(frame)
    view = plot_communication_dotmap(frame)
    try:
        keys = ['sender','receiver','hmdb_id','sensor_gene']
        plotted = view['plot_events']
        assert len(plotted) == len(frame) and not plotted.duplicated(keys).any()
        assert set(map(tuple,plotted[keys].values)) == set(map(tuple,frame[keys].values))
        assert plotted.cell_mesh_score.isna().sum() == 1
        assert view['selection']['n_in_scope'] == 105
        dots = view['dot_table']
        assert dots.id.is_unique and not dots.duplicated(['x','y']).any()
        assert len(view['x_keys']) == 105 and len(view['y_keys']) == 6
        for dot in dots.itertuples():
            assert view['x_keys'][dot.x] == (dot.receiver, dot.sensor_gene)
            assert view['y_keys'][dot.y] == (dot.sender, dot.hmdb_id)
        pd.testing.assert_frame_equal(frame,before)
    finally: plt.close(view['fig'])
    text, record = communication_dotmap_html(view)
    assert '<script>unsafe()' not in text
    data = payload(text,'selection-dotmap-data')
    assert len(data['records']) == record['n_plotted'] == 105
    assert text.count('<svg ') == record['n_figures'] == 1
    assert data['x_fields'] == ['receiver','sensor_gene']
    assert data['y_fields'] == ['sender','hmdb_id']
    assert '<select' not in text and '<button' not in text and 'pages' not in data


@pytest.mark.parametrize('top_n',[True,0,-1,1.5,'3'])
def test_invalid_matrix_limits(top_n):
    with pytest.raises(ValueError):
        plot_communication_dotmap(path_frame(),top_n=top_n)


def test_pair_scope_retains_unscored_and_does_not_apply_hidden_qc_or_fdr():
    frame=path_frame()
    frame.loc[frame.hmdb_id.eq('HMDB0000002'),'cell_mesh_score']=np.nan
    all_scope=communication_scope(frame)
    assert len(all_scope['records']) == 105
    limited=communication_scope(frame,top_n=1)
    assert limited['pairs'] == [('HMDB0000001','R<script>unsafe()</script>')]
    result=CellMeshResult(frame,pd.DataFrame(),pd.DataFrame(),{})
    selection=result.filter_events(source='all')
    view=plot_communication_dotmap(selection)
    try:
        text=selection_plots_html(selection, dotmap=view, secondary=True)
    finally:
        plt.close(view['fig'])
    plots={r['kind']:r for r in payload(text,'selection-plot-data')}
    assert plots['dotmap']['n_plotted'] == 105
    assert plots['dotmap']['n_missing_score'] == 52
    assert plots['network']['edge_event_total'] == 53
    assert 'counts' not in plots
    with pytest.raises(ValueError, match='unique'):
        communication_scope(pd.concat([frame,frame.iloc[:1]]))
    with pytest.raises(ValueError, match='aggregate'):
        communication_scope(frame.assign(sample='x'))


def test_count_heatmap_requires_multiple_relations_and_unequal_observed_counts():
    frame=path_frame().iloc[:4].copy()
    frame.loc[1,['sender','receiver']]=frame.loc[0,['sender','receiver']].values
    result=CellMeshResult(frame,pd.DataFrame(),pd.DataFrame(),{})
    text=selection_plots_html(result.filter_events(source='all'),secondary=True)
    plots={r['kind']:r for r in payload(text,'selection-plot-data')}
    assert plots['counts']['count_total'] == 4
    assert plots['network']['edge_event_total'] == 4


@pytest.mark.parametrize('mode',['pooled_stratified','sample_aware'])
def test_violin_grouping_and_append_cache_without_implicit_input(inputs,tmp_path,monkeypatch,mode):
    result=run(inputs,sample_mode=mode,n_perms=0,retain_input=True)
    focus=result.filter_events(source='all')
    paths=result.save(tmp_path/'a',selection=focus,include_input=True)
    expected=block(paths['filtered_report'])
    assert expected.count('data-expression-index=') == 2
    assert expected.count('<svg') == 2
    assert ('汇集' in expected) == (mode=='sample_aware')
    loaded=CellMeshResult.load(paths['result'])
    def forbidden(*args,**kwargs):
        raise AssertionError('must not implicitly load expression archive')
    monkeypatch.setattr(loaded,'load_input',forbidden)
    event=loaded.events[['sender','receiver','hmdb_id','sensor_gene']].iloc[0].to_dict()
    paths=loaded.add_event_reports(loaded.select_events(metabolite_sensor_pairs=[(event["hmdb_id"], event["sensor_gene"])], source="all"))
    assert block(paths['filtered_report']) == expected
    # Replacement must discard figures with a stale selection scope.
    paths=loaded.add_selection(loaded.filter_events(source='all',sender_labels=[]),replace=True)
    assert START not in paths['filtered_report'].read_text()
    assert not plt.get_fignums()


def test_missing_input_notice_can_be_filled_on_later_append(inputs,tmp_path):
    result=run(inputs,n_perms=0)
    paths=result.save(tmp_path/'a',selection=result.filter_events(source='all'))
    assert '未提供细胞级表达输入' in block(paths['filtered_report'])
    event=result.events[['sender','receiver','hmdb_id','sensor_gene']].iloc[0].to_dict()
    paths=result.add_event_reports(result.select_events(metabolite_sensor_pairs=[(event["hmdb_id"], event["sensor_gene"])], source="all"), adata=inputs[0])
    assert block(paths['filtered_report']).count('<svg') == 2


@pytest.mark.parametrize('report',[True,False])
def test_full_csv_opt_in_and_fresh_save_does_not_inherit(report,tmp_path):
    result=make_result(n_perms=0)
    paths=result.save(tmp_path/'a',report=report,include_all_events=True)
    original=paths['events_csv'].read_bytes()
    restored=CellMeshResult.load(paths['result'])
    if report:
        paths['events_csv'].unlink()
        paths=restored.add_selection(restored.filter_events(source='all'))
        assert paths['events_csv'].read_bytes()==original
    fresh=restored.save(tmp_path/'fresh',report=report)
    assert 'events_csv' not in fresh
    assert not (tmp_path/'fresh.all_events.csv').exists()
    with h5py.File(fresh['result']) as h:
        assert _read(h['metadata/delivery'])['include_all_events'] is False
        pd.testing.assert_frame_equal(_read(h['results/events']),result.events)


def test_disable_full_csv_leaves_existing_file_unlinked_and_legacy_append_does_not_restore(tmp_path):
    result=make_result(n_perms=0)
    paths=result.save(tmp_path/'a',include_all_events=True)
    full=paths['events_csv']; content=full.read_bytes()
    paths=result.save(tmp_path/'a')
    assert full.read_bytes()==content and 'events_csv' not in paths
    assert full.name not in paths['report'].read_text()
    # Legacy delivery had no flag: an append migrates links, never requests full CSV.
    result.delivery.pop('include_all_events')
    with h5py.File(paths['result'],'r+') as h:
        del h['metadata/delivery']
        _write(h['metadata'],'delivery',result.delivery,include_null_scores=False)
    paths=result.add_selection(result.filter_events(source='all'))
    assert 'events_csv' not in paths and full.read_bytes()==content
    assert full.name not in paths['report'].read_text()


def test_standalone_full_csv_only_for_general_report(tmp_path):
    result=make_result(n_perms=0)
    for options in [dict(selection=result.filter_events(source='all')),dict(selected_events=result.select_events(metabolite_sensor_pairs=[(SELECTOR["hmdb_id"], SELECTOR["sensor_gene"])], source="all"))]:
        with pytest.raises(ValueError,match='general'):
            result.to_html(tmp_path/'bad.html',include_all_events=True,**options)
        assert not list(tmp_path.iterdir())
    general=result.to_html(tmp_path/'default.html')
    assert not (tmp_path/'default.all_events.csv').exists()
    assert (tmp_path/'default.positive_events.csv').exists()
    explicit=result.to_html(tmp_path/'explicit.html',include_all_events=True)
    assert 'explicit.all_events.csv' in explicit.read_text()
    pd.testing.assert_frame_equal(pd.read_csv(tmp_path/'explicit.all_events.csv'),result.events,check_dtype=False)


def test_dot_encoding_matches_score_and_chosen_statistic_without_mutating_data():
    frame=path_frame().iloc[:5].copy()
    frame['cell_mesh_score']=[0,.5,.9,np.nan,.9]
    frame['perm_pvalue']=[.1,.01,0,np.nan,.01]
    before=frame.copy(deep=True)
    view=plot_communication_dotmap(frame,fdr_col='perm_pvalue')
    try:
        encoding=view['encoding']
        assert encoding['statistic']=='perm_pvalue'
        assert dot_area(0,encoding) < dot_area(.5,encoding) < dot_area(.9,encoding)
        assert len({dot_color(v,encoding) for v in [.1,.01,0,None]})==4
        dots=view['dot_table']
        assert (dots.marker=='x').sum()==1
        assert encoding['has_zero'] and encoding['missing_statistic']
        assert view['colorbar'] is not None  # A real rectangular continuous colorbar.
        assert view['layout']['legend_kind']=='continuous'
        pd.testing.assert_frame_equal(frame,before)
    finally: plt.close(view['fig'])


@pytest.mark.parametrize('stat',[.05,0,np.nan])
def test_constant_zero_and_unavailable_statistics_have_honest_static_keys(stat):
    frame=path_frame().iloc[:3].copy()
    frame['fdr_sensor_type']=stat
    view=plot_communication_dotmap(frame)
    try:
        assert view['dot_table'].color.nunique()==1
        assert view['colorbar'] is None  # Never invent a gradient for a constant.
        labels=' '.join(text.get_text() for text in view['legend_ax'].texts)
        if pd.isna(stat): assert 'unavailable' in labels
        elif stat==0: assert '= 0' in labels
        else: assert 'All available values equal' in labels
    finally: plt.close(view['fig'])


def test_readable_full_labels_and_legend_numbers():
    from cellmesh._communication_dotmap import legend_numbers, node_label_lines
    frame=path_frame().iloc[:5].copy()
    frame['metabolite']='A long metabolite name with several meaningful words'
    view=plot_communication_dotmap(frame,metabolite_label='name')
    try:
        labels=[tick.get_text() for tick in view['ax'].get_yticklabels()]
        assert all('HMDB' in text for text in labels)  # Disambiguate identical names with different IDs.
        assert view['labels']['context']==''
        assert all('R<script>unsafe()</script>' in tick.get_text() for tick in view['ax'].get_xticklabels())
        assert all('several meaningful words' in ' '.join(text.splitlines()) for text in labels)
        assert all('\n' not in text for text in labels)
        assert all(t.get_fontsize()>=12 for t in view['ax'].get_xticklabels())
    finally: plt.close(view['fig'])
    assert legend_numbers([.141918,.310998,.491789])==['0.142','0.311','0.492']
    assert len(set(legend_numbers([.001000001,.001000002])))==2
    assert float(legend_numbers([1e-30])[0])>0
    label='Long metabolite name with words and punctuation'
    assert ' '.join(node_label_lines(label))==label


def test_multiple_metabolites_and_sensors_share_one_matrix_with_correct_coordinates():
    # Explicitly specified coordinates catch an accidental return to cell-pair × metabolite-sensor axes.
    identities=[('S1','R1','HMDB1','G1'),('S1','R2','HMDB1','G2'),
                ('S2','R1','HMDB2','G1'),('S1','R1','HMDB2','G1'),('S2','R2','HMDB1','G2')]
    frame=pd.DataFrame(identities,columns=['sender','receiver','hmdb_id','sensor_gene'])
    frame=frame.assign(metabolite='Shared display name',cell_mesh_score=.5,fdr_sensor_type=.03)
    view=plot_communication_dotmap(frame)
    try:
        assert view['x_keys']==[('R1','G1'),('R2','G2')]
        assert view['y_keys']==[('S1','HMDB1'),('S1','HMDB2'),('S2','HMDB1'),('S2','HMDB2')]
        expected=dict(zip(identities,[(0,0),(1,0),(0,3),(0,1),(1,2)]))
        for dot in view['dot_table'].itertuples():
            assert (dot.x,dot.y)==expected[(dot.sender,dot.receiver,dot.hmdb_id,dot.sensor_gene)]
        assert view['colorbar'] is None
        assert len(view['ax'].collections)==len(frame)  # All event markers share one data axis.
        assert len(view['plot_events'])==5
    finally: plt.close(view['fig'])


def test_dotmap_full_composite_ticks_even_for_shared_molecules_and_readable_colorbar():
    frame=path_frame().iloc[:5].copy()
    frame['hmdb_id']='HMDB0000961'
    frame['metabolite']='Farnesyl pyrophosphate'
    frame['sensor_gene']='CYSLTR1'
    frame['fdr_sensor_type']=[.05,.01,.001,0,np.nan]
    view=plot_communication_dotmap(frame)
    try:
        for (sender, hmdb), tick in zip(view['y_keys'], view['ax'].get_yticklabels()):
            assert tick.get_text() == f'{sender} – {hmdb}'
            assert tick.get_rotation()==0
        for (receiver, sensor), tick in zip(view['x_keys'], view['ax'].get_xticklabels()):
            assert tick.get_text() == f'{receiver} – {sensor}'
        assert view['context_artist'] is None
        assert all(t.get_rotation()==45 for t in view['ax'].get_xticklabels())
        assert view['title_artist'].get_fontsize()==14
        assert view['title_artist'].get_text()=='Selected cell–cell communication events (n = 5)'
        assert view['title_artist'].get_ha()=='center'
        assert view['colorbar'] is not None
        assert len(view['layout']['colorbar_ticks'])==5
        fig=view['fig'];fig.canvas.draw()
        renderer=fig.canvas.get_renderer()
        bar=view['colorbar'].ax.get_window_extent(renderer)
        assert bar.width/fig.dpi==pytest.approx(.28)
        assert bar.height/fig.dpi==pytest.approx(1.8)
        assert not bar.overlaps(view['ax'].get_tightbbox(renderer))
        ticks=view['ax'].get_yticklabels()
        for first, second in zip(ticks, ticks[1:]):
            assert not first.get_window_extent(renderer).overlaps(second.get_window_extent(renderer))
        assert view['ax'].get_window_extent(renderer).height < view['legend_ax'].get_window_extent(renderer).height
        legend_text=' '.join(t.get_text() for t in view['legend_ax'].texts)
        assert '−log10(BH-adjusted p)' in legend_text
        assert 'BH correction: within each sensor type' in legend_text
        assert 'Dot area' not in legend_text
        assert 'BH-adjusted p = 0' in legend_text and 'Statistic unavailable' in legend_text
        from matplotlib.colors import to_hex
        for row in frame.loc[frame.fdr_sensor_type.gt(0)].itertuples():
            expected=to_hex(view['colorbar'].cmap(view['colorbar'].norm(-np.log10(row.fdr_sensor_type))))
            assert dot_color(row.fdr_sensor_type,view['encoding'])==expected
    finally: plt.close(view['fig'])


def test_metabolite_label_name_changes_only_display_and_falls_back_to_id():
    frame=path_frame().iloc[:5].copy()
    frame.loc[frame.hmdb_id.eq('HMDB0000001'),'metabolite']='Methionine'
    frame.loc[frame.hmdb_id.eq('HMDB0000002'),'metabolite']=None
    ids=plot_communication_dotmap(frame)
    names=plot_communication_dotmap(frame,metabolite_label='name')
    try:
        assert ids['x_keys']==names['x_keys'] and ids['y_keys']==names['y_keys']
        pd.testing.assert_frame_equal(ids['dot_table'],names['dot_table'])
        for (sender,hmdb),tick in zip(names['y_keys'],names['ax'].get_yticklabels()):
            label='Methionine' if hmdb=='HMDB0000001' else hmdb
            assert tick.get_text()==f'{sender} – {label}'
    finally:
        plt.close(ids['fig']);plt.close(names['fig'])


@pytest.mark.parametrize('label',['short','HMDB',None,[]])
def test_invalid_metabolite_label(label):
    with pytest.raises(ValueError,match='metabolite_label'):
        plot_communication_dotmap(path_frame(),metabolite_label=label)


def test_constant_tiny_probability_swatch_keeps_original_value():
    frame=path_frame().iloc[:2].copy().assign(fdr_sensor_type=1e-320)
    view=plot_communication_dotmap(frame)
    try:
        assert view['colorbar'] is None
        assert any('1e-320' in t.get_text() for t in view['legend_ax'].texts)
        assert view['dot_table'].color.nunique()==1
    finally: plt.close(view['fig'])
