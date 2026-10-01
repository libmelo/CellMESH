"""Offline report contracts: all rows, bounded DOM, stable identities and scope."""
import copy
from html.parser import HTMLParser
import json
import os
from pathlib import Path
import re
import subprocess
import time

import numpy as np
import pandas as pd
import pytest

from cellmesh._report_browser import event_browser, explanation_entries
from cellmesh.reporting import render_html_report
from test_reporting_and_persistence import make_result, SELECTOR, assert_nested_equal


@pytest.mark.skipif(os.environ.get('CELLMESH_BROWSER_TESTS') != '1', reason='optional system WebKit browser')
def test_shared_event_browser_empty_and_large_table_interactions(tmp_path):
    """Exercise event-table JS separately from tier-3 record-table acceptance."""
    from cellmesh._report_pages import document
    server = subprocess.Popen(['/usr/bin/broadwayd', ':28'], stdout=subprocess.DEVNULL,
                              stderr=subprocess.DEVNULL)
    env = dict(os.environ, GDK_BACKEND='broadway', BROADWAY_DISPLAY=':28', LIBGL_ALWAYS_SOFTWARE='1')
    tests = Path(__file__).parent
    try:
        time.sleep(.5)
        for n in [0, 10001]:
            frame = pd.DataFrame([dict(sender='A', receiver='B', hmdb_id=f'HMDB{i+1:07d}',
                                       sensor_gene='R', cell_mesh_score=.2) for i in range(n)],
                                 columns=['sender', 'receiver', 'hmdb_id', 'sensor_gene', 'cell_mesh_score'])
            path = tmp_path / f'events-{n}.html'
            path.write_text(document('Browser test', [event_browser(frame, [], show_details=False)], scope='general'))
            done = subprocess.run(['/usr/bin/python3', str(tests/'report_browser_webkit.py'), str(path),
                                   str(tests/'report_browser_boundaries.js')],
                                  env=env, capture_output=True, text=True, timeout=60)
            assert done.returncode == 0, done.stdout + done.stderr
            assert json.loads(done.stdout.strip().splitlines()[-1])['passed'] > 0
    finally:
        server.terminate()
        server.wait(timeout=5)


class ReportParser(HTMLParser):
    def __init__(self, text):
        super().__init__()
        self.tags = []
        self.feed(text)

    def handle_starttag(self, tag, attrs):
        self.tags.append((tag, dict(attrs)))


def payload(text, key="event-data"):
    match = re.search(r'<script type="application/json" id="' + key + r'">(.*?)</script>', text, re.S)
    assert match
    return json.loads(match.group(1))


@pytest.fixture(scope="module")
def captured():
    return make_result(n_perms=0)


def many_events(result, n=125):
    """Synthetic table-only variants; retain the four real events for explanations."""
    result = copy.deepcopy(result)
    original = pd.DataFrame(result.events)
    extra = [original.iloc[[0]].assign(sensor_gene=f"SENSOR_{i:03d}") for i in range(n - len(original))]
    result.events = pd.concat([original, *extra], ignore_index=True)
    # Includes small numeric strings, NA, non-ASCII, case-insensitive search and
    # literal markup. None of these browser fixtures claim new biological scores.
    result.events["fdr_sensor_type"] = [.02, .1, .003, np.nan] + [.5] * (n - 4)
    result.events["fdr_global"] = result.events.fdr_sensor_type.map(lambda x: None if pd.isna(x) else str(x))
    result.events.loc[n - 1, "metabolite"] = '末页代谢物 <b>literal</b> </script><script>window.PWNED=1</script>'
    result.events["extra_annotation"] = None
    return result


@pytest.mark.parametrize("n", [0, 1, 49, 50, 51, 100, 125, 1001])
def test_browser_embeds_every_row_but_only_first_page_html(n):
    frame = pd.DataFrame([dict(sender="A", receiver="B", hmdb_id=f"HMDB{i:07d}",
                               metabolite="同名", sensor_gene="R", cell_mesh_score=i / max(n, 1))
                          for i in range(n)], columns=["sender", "receiver", "hmdb_id", "metabolite",
                                                        "sensor_gene", "cell_mesh_score"])
    text = event_browser(frame, [])
    data = payload(text)
    assert len(data["rows"]) == n
    assert data["columns"] == list(frame.columns)
    assert data["page_size"] == 50
    assert len([attrs for tag, attrs in ReportParser(text).tags if "data-event-index" in attrs]) == min(n, 50)
    if n:
        assert data["rows"][-1]["values"][data["columns"].index("hmdb_id")] == f"HMDB{n - 1:07d}"
    assert 'id="event-browser"' in text


def test_payload_preserves_precision_and_missingness():
    frame = pd.DataFrame([dict(sender="001", receiver="B", hmdb_id="HMDB0000001", sensor_gene="R",
                               cell_mesh_score=0., fdr_sensor_type=pd.NA, perm_pvalue="0.00000001",
                               exact_id=2**60 + 1, note="", passes_min_cells=False)])
    data = payload(event_browser(frame, []))
    values = dict(zip(data["columns"], data["rows"][0]["values"]))
    assert values["sender"] == "001"
    assert values["exact_id"] == str(2**60 + 1)
    assert values["fdr_sensor_type"] is None
    assert values["cell_mesh_score"] == "0.0"
    assert values["passes_min_cells"] == "False"
    assert values["note"] == ""
    assert values["perm_pvalue"] == "0.00000001"
    assert data["columns"].index("perm_pvalue") in data["numeric"]


def test_inert_json_and_html_escape_names_headers_and_unicode(captured):
    result = many_events(captured)
    bad = '</ScRiPt><img src=x onerror="window.PWNED=1">\u2028\u2029&'
    result.events[bad] = bad
    # Exercise standalone rendering without filesystem side effects.
    text = render_html_report(result, selection=result.filter_events(source="all"), detail_n=0, positive_threshold=1)
    data = payload(text, "interest-data")
    assert bad in data["columns"]
    assert data["rows"][0]["values"][data["columns"].index(bad)] == bad
    assert bad not in text
    tags = ReportParser(text).tags
    scripts = [attrs for tag, attrs in tags if tag == "script"]
    assert {attrs.get("id") for attrs in scripts} == {"interest-plot-data", "sample-heatmap-data", "interest-data", "event-browser-script", "event-sample-script", "interest-dotmap-data", "communication-features-script"}
    assert all("src" not in attrs for attrs in scripts)
    assert not [attrs for tag, attrs in tags if tag == "img"]
    assert not any(key.lower().startswith("on") for _, attrs in tags for key in attrs)
    assert any("window.PWNED=1" in row["values"][data["columns"].index("metabolite")] for row in data["rows"])


@pytest.mark.parametrize("detail_n,expected", [(0, 0), (1, 0), (10, 0), (None, 0)])
def test_explanation_directory_counts_links_and_default_open(captured, detail_n, expected):
    focus = captured.filter_events(source="all")
    before = copy.deepcopy(captured)
    text = render_html_report(captured, selection=focus, detail_n=detail_n)
    data = payload(text, "interest-data")
    assert len(data["details"]) == expected
    assert sum(bool(row["details"]) for row in data["rows"]) == expected
    assert f"共 4 条事件，其中 {expected} 条可查看对应通讯报告。" in text
    tags = ReportParser(text).tags
    sections = [attrs for tag, attrs in tags if tag == "details" and attrs.get("class") == "event-explanation"]
    assert len(sections) == expected
    assert ["open" in attrs for attrs in sections] == [False] * expected
    assert len(re.findall(r'<div class="selection-plot"><svg ', text)) == 2  # network and sample heatmap; matrix uses static page containers
    assert 'id="interest-dotmap"' in text and 'class="flow-legends"' not in text
    ids = [attrs["id"] for attrs in sections]
    for row in data["rows"]:
        for index in row["details"]:
            assert data["details"][index]["id"] in ids
    assert_nested_equal(captured.events, before.events)
    assert captured.parameters == before.parameters


def test_manual_sample_scopes_and_unscored_candidates_remain_in_directory(captured, tmp_path):
    second = dict(SELECTOR, sample="S2")
    unscored = dict(SELECTOR, sensor_gene="ABSENT_R")
    with pytest.raises(TypeError, match='retired'):
        captured.to_html(tmp_path/'duplicate.html', events=[SELECTOR, second])
    text = render_html_report(captured, events=[second])
    assert '<option value="1" selected>S2</option>' in text
    assert len(payload(text, 'detail-sample-events-data')['rows']) == 3
    missing = render_html_report(captured, events=[unscored])
    assert '未生成该事件' in missing and 'ABSENT_R' in missing
    assert 'id="interest-data"' not in missing


def test_normalized_identity_links_do_not_merge_same_names():
    selectors = [dict(sender=" A ", receiver="B", hmdb_id="hmdb0000001", sensor_gene=" R ")]
    frame = pd.DataFrame([dict(sender="A", receiver="B", hmdb_id=identifier, sensor_gene="R",
                               metabolite="same name", cell_mesh_score=.5)
                          for identifier in ["HMDB0000001", "HMDB0000002"]])
    data = payload(event_browser(frame, explanation_entries(selectors)))
    assert [row["details"] for row in data["rows"]] == [[0], []]


def test_full_and_selected_delivery_scopes_unchanged(captured, tmp_path):
    result = many_events(captured)
    before = copy.deepcopy(result.events)
    focus = result.filter_events(source="all", sender_labels=["A"])
    files = result.save(tmp_path / "analysis", report=True, selection=focus, detail_n=0, include_all_events=True)
    selected_path = files["filtered_csv"]
    data = payload(files["filtered_report"].read_text(), "interest-data")
    assert len(data["rows"]) == len(focus.events) > 100
    assert len(pd.read_csv(files["events_csv"])) == len(before)
    assert len(pd.read_csv(selected_path)) == len(data["rows"])
    pd.testing.assert_frame_equal(result.events, before, check_exact=True)
    assert data["rows"][-1]["values"][data["columns"].index("sensor_gene")] == focus.events.iloc[-1].sensor_gene
    assert set(path.name for path in tmp_path.iterdir()) == {
        "analysis.cellmesh.h5", "analysis.all_events.csv", "analysis.positive_events.csv", "analysis.report.html", "analysis.filtered_events.csv", "analysis.filtered.report.html", "analysis.reports"}


def test_unfiltered_report_shows_positive_table_and_no_explanations(captured):
    text = render_html_report(captured)
    data = payload(text)
    assert len(data["rows"]) == 0
    assert data["details"] == []
    assert "筛选条件与范围" not in text
    assert "共 0 条事件。" in text
    assert 'id="report-explanations"' not in text


def test_bad_header_aborts_before_overwriting_html(captured, tmp_path):
    result = copy.deepcopy(captured)
    result.events.columns = ["sender"] * len(result.events.columns)
    path = tmp_path / "report.html"
    path.write_text("previous report")
    with pytest.raises((ValueError, KeyError)):
        result.to_html(path)
    assert path.read_text() == "previous report"
