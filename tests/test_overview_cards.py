"""Scope-aware counts and displayed precision must preserve the decision rule."""
import re

import numpy as np
import pytest

from cellmesh.reporting import render_html_report
from cellmesh._report_layout import selection_overview_html
from cellmesh.positivity import positive_calls, resolve_positive, selection_table
from test_positive_delivery import result
from test_reporting_and_persistence import SELECTOR


def card_texts(text):
    return re.findall(r'<div class="overview-card[^\"]*">(.*?)</div>',text,re.S)


def test_general_cards_show_positive_and_full_denominator(result):
    text=render_html_report(result)
    cards=card_texts(text)
    assert len(cards)==2
    assert '达标事件' in cards[0] and '<strong>2</strong>' in cards[0]
    assert '50%' in cards[0]
    assert '全部已生成事件' in cards[1] and '<strong>4</strong>' in cards[1]
    assert '无法判定的原因' in text and '未达标' in text


def test_selected_all_positive_and_mixed_counts_do_not_duplicate_cards(result):
    for source in ['positive','all']:
        focus=result.filter_events(source=source)
        spec=resolve_positive(result)
        text=selection_overview_html(result,focus,selection_table(focus,spec),spec)
        assert len(card_texts(text))==1
        assert ('全部达标' in text)==(source=='positive')
        if source=='all':
            assert '未达标' in text and '无法判定' in text
    focus=result.filter_events(source='all',sender_labels=[])
    text=selection_overview_html(result,focus,selection_table(focus,spec),spec)
    assert '<strong>0</strong>' in text and '全部达标' not in text


@pytest.mark.parametrize('options,reason',[
    ({},''),
    ({'positive_qc_only':True},'未通过细胞数 QC'),
    ({'positive_min_score':1.},'通讯评分低于门槛'),
])
def test_event_card_status_uses_enabled_criterion(result,options,reason):
    # The chosen row has p/BH below threshold but fails QC and has score < 1.
    selector=result.events[['sender','receiver','hmdb_id','sensor_gene']].iloc[0].to_dict()
    result.events.loc[0,'cell_mesh_score']=.1234567890123
    text=render_html_report(result,events=[selector],**options)
    cards=card_texts(text)
    assert len(cards)==3  # score, chosen statistic, available sample coverage
    assert '<strong title="原值：0.1234567890123">0.123457</strong>' in cards[0]
    if reason:
        assert reason in text and '未达标原因' in text
    else:
        assert '未达标原因' not in text


def test_missing_statistics_remain_undetermined_not_zero(result):
    result.events['fdr_sensor_type']=np.nan
    text=render_html_report(result)
    assert '尚未判定' in card_texts(text)[0]
    selector=result.events[['sender','receiver','hmdb_id','sensor_gene']].iloc[0].to_dict()
    text=render_html_report(result,events=[selector])
    assert '无法判定' in text and '>—</strong>' in text
    assert result.events.fdr_sensor_type.isna().all()
