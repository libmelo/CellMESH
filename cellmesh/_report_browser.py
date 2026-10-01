"""Offline event-table browsing and links to precomputed explanations.

All rows are embedded as inert JSON; only one page becomes table DOM. This
module changes presentation only, never scientific values or selection scope.
"""
from __future__ import annotations

import html
import json
import re

import pandas as pd

from .selection import event_id, event_identity as _identity


_KEYS = ("sender", "receiver", "hmdb_id", "sensor_gene")
_SEARCH_FIELDS = ("sender", "receiver", "hmdb_id", "metabolite", "sensor_gene")
_PAGE_SIZE = 50
_COLUMN_LABELS = {
    "sender": "发送细胞", "receiver": "接收细胞", "hmdb_id": "HMDB 编号",
    "metabolite": "代谢物", "sensor_gene": "受体基因", "sensor_type": "受体类型",
    "cell_mesh_score": "通讯评分", "perm_pvalue": "置换 P 值",
    "fdr_sensor_type": "BH 校正 p 值（按受体类型）", "fdr_global": "BH 校正 p 值（全局）",
    "positive_status": "阳性判定状态",
    "passes_min_cells": "细胞数质控", "metabolite_availability": "代谢物可用性",
    "sensor_score": "受体评分", "sensor_expr_frac": "受体表达比例",
    "n_samples_coobserved": "可计算样本数", "n_samples_positive": "正评分样本数",
    "n_samples_passing_min_cells": "通过质控样本数", "event_prevalence": "正评分样本比例",
    "min_cells_pass_prevalence": "质控通过样本比例",
}


def page_entries(selectors, hrefs):
    entries = []
    for selector in selectors:
        sender, receiver, hmdb, sensor = _identity(selector)
        key = event_id(selector)
        entries.append(dict(id=key, label=f"{sender} → {hmdb} / {sensor} → {receiver}",
                            scope="单事件报告", identity=(sender, receiver, hmdb, sensor),
                            href=str(hrefs[key])))
    return entries


def explanation_entries(selectors):
    entries = []
    for i, selector in enumerate(selectors):
        sender, receiver, hmdb, sensor = _identity(selector)
        sample = selector.get("sample")
        scope = "全部样本" if sample is None else f"样本 {str(sample).strip()}"
        label = f"{sender} → {hmdb} / {sensor} → {receiver}"
        if sample is not None:
            label += f" · {scope}"
        entries.append(dict(id=f"event-explanation-{i + 1}", label=label, scope=scope,
                            identity=(sender, receiver, hmdb, sensor)))
    return entries


def _cell_text(value):
    # Strings preserve IDs and large integers exactly in the browser. Do not
    # serialize DataFrame attrs (which may hold permutation matrices).
    if pd.api.types.is_scalar(value) and pd.isna(value):
        return None
    return str(value)


def _script_json(payload):
    # A JSON script element is still parsed by HTML first: escape every '<',
    # including mixed-case closing script tags in user-supplied annotations.
    encoded = json.dumps(payload, ensure_ascii=False, allow_nan=False, separators=(",", ":"))
    for char, replacement in (("&", "\\u0026"), ("<", "\\u003c"), (">", "\\u003e"),
                              ("\u2028", "\\u2028"), ("\u2029", "\\u2029")):
        encoded = encoded.replace(char, replacement)
    return encoded


def _links(details, entries, code=None):
    if not details:
        text = '<span class="event-unexplained">尚未生成对应通讯报告</span>'
        if code:
            text += '<details class="event-request"><summary>追加代码</summary><pre>' + html.escape(code) + '</pre></details>'
        return text
    return " · ".join('<a data-detail-link href="' + html.escape(entries[i].get("href", "#" + entries[i]["id"]), quote=True) +
                      '">查看' + (entries[i].get('scope', '通讯报告') if "href" in entries[i] else '解释（' + html.escape(entries[i]["scope"]) + '）') + '</a>' for i in details)


def event_browser(frame, entries, *, prefix="event", title="阳性事件列表", show_details=True,
                  event_table=True, compact=False):
    """Return the complete event browser; original column keys stay in its data."""
    columns = [str(name) for name in frame.columns]
    if len(set(columns)) != len(columns):
        raise ValueError("Report event table requires unique column names")
    by_key = {}
    for i, entry in enumerate(entries):
        by_key.setdefault(entry["identity"], []).append(i)
    positions = {key: columns.index(key) for key in _KEYS} if event_table else {}
    rows = []
    for values in frame.itertuples(index=False, name=None):
        selector = {key: values[index] for key, index in positions.items()}
        identity = _identity(selector) if event_table else None
        code = ("selected = result.select_events(metabolite_sensor_pairs=" + repr([(identity[2], identity[3])]) + ")\nresult.add_event_reports(selected)") if show_details and event_table else None
        rows.append(dict(values=[_cell_text(value) for value in values],
                         details=by_key.get(identity, []), code=code,
                         event_id=event_id(selector) if show_details and event_table else None))
    numeric = [i for i, name in enumerate(frame.columns)
               if ((pd.api.types.is_numeric_dtype(frame[name].dtype)
                    and not pd.api.types.is_bool_dtype(frame[name].dtype)
                    and not pd.api.types.is_complex_dtype(frame[name].dtype))
                   or str(name) in {"cell_mesh_score", "perm_pvalue", "fdr_sensor_type", "fdr_global"})]
    explained = sum(bool(row["details"]) for row in rows)
    unit = '事件' if event_table else '记录'
    payload = dict(columns=columns, labels=[_COLUMN_LABELS.get(name, name) for name in columns],
                   row_unit=unit, compact=compact,
                   numeric=numeric, rows=rows, page_size=_PAGE_SIZE, show_details=show_details,
                   search_columns=([columns.index(name) for name in _SEARCH_FIELDS if name in columns]
                                   if event_table else list(range(len(columns)))),
                   details=[{k: entry[k] for k in ("id", "label", "scope", "href") if k in entry} for entry in entries])
    parts = [f'<section class="event-browser' + (' event-browser-compact' if compact else '') +
             f'" data-prefix="{prefix}" id="event-browser" aria-label="{html.escape(title, quote=True)}">']
    if not compact:
        parts += [f'<p id="event-scope">共 {len(rows):,} 条{unit}' +
                  (f'，其中 {explained:,} 条可查看对应通讯报告。' if show_details else '。') + '</p>',
                  '<p class="report-note">搜索、排序和翻页只调整表格视图，不改变图形范围、原始评分、CSV 或解释对象。</p>']
    parts += [
             '<div id="event-controls" class="event-controls" hidden>',
             '<label>搜索范围 <select id="event-search-field"><option value="">' +
             ('全部名称与编号' if event_table else '全部字段') + '</option>']
    for index in payload["search_columns"]:
        parts.append(f'<option value="{index}">{html.escape(payload["labels"][index])}</option>')
    parts += ['</select></label>',
              '<label>搜索 <input id="event-search" type="search" placeholder="' +
              ('细胞、HMDB、代谢物名称或受体' if event_table else '输入关键词') + '" '
              'aria-controls="event-rows"></label>',
              '<label>每页 <select id="event-page-size"><option>25</option><option selected>50</option>'
              '<option>100</option></select> 条</label>',
              '<button type="button" id="event-reset">恢复初始视图</button></div>',
              f'<p id="event-status" role="status" aria-live="polite">共 {len(rows)} 条{unit}；'
              f'静态预览前 {min(len(rows), _PAGE_SIZE)} 条。</p>',
              '<div class="table event-table"><table id="event-table"><caption' +
              (' class="visually-hidden"' if compact else '') + '>' + html.escape(title) + '</caption><thead><tr>',
              '<th scope="col">' + (('事件标识' if show_details else '报告编号') if event_table else '序号') + '</th>' +
              ('<th scope="col">详细解释</th>' if show_details else '')]
    for i, name in enumerate(columns):
        parts.append(f'<th scope="col" aria-sort="none"><button type="button" data-sort-column="{i}" '
                     f'title="{html.escape(name)}" disabled>{html.escape(payload["labels"][i])}</button></th>')
    parts.append('</tr></thead><tbody id="event-rows">')
    for i, row in enumerate(rows[:_PAGE_SIZE]):
        parts += [f'<tr data-event-index="{i}"><td>{row["event_id"] or i + 1}</td>' +
                  (f'<td>{_links(row["details"], entries, row["code"])}</td>' if show_details else ''),
                  ''.join('<td' + (' class="event-identity"' if columns[j] in _KEYS else '') + '>' +
                          ('<span class="event-na">NA</span>' if value is None else html.escape(value)) +
                          '</td>' for j, value in enumerate(row["values"])), '</tr>']
    if not rows:
        parts.append(f'<tr><td colspan="{len(columns) + 1 + int(show_details)}">无对应{unit}。</td></tr>')
    parts += ['</tbody></table></div>',
              '<div id="event-pager" class="event-controls" hidden>',
              '<button type="button" id="event-prev">上一页</button>',
              '<label>跳至第 <input id="event-page" type="number" min="1" value="1" step="1" '
              'aria-controls="event-rows"> 页</label>',
              '<button type="button" id="event-next">下一页</button></div>']
    if not compact:
        parts.append('<noscript><p>启用 JavaScript 后可离线搜索、排序和翻页，查看全部内嵌' + unit + '。</p></noscript>')
    parts += ['<script type="application/json" id="event-data">' + _script_json(payload) + '</script>', '</section>']
    content = ''.join(parts)
    if prefix != "event":
        # Namespace markup only; never rewrite event labels or JSON values.
        content = re.sub(r'<[^>]+>', lambda m: re.sub(r'(\b(?:id|aria-controls)="|href="#)event-(?!explanation-)',
                         lambda a: a[1] + prefix + '-', m[0]), content)
    return content


BROWSER_CSS = """
.visually-hidden{position:absolute!important;width:1px;height:1px;padding:0;overflow:hidden;clip:rect(0,0,0,0);white-space:nowrap;border:0}
.event-controls{display:flex;flex-wrap:wrap;align-items:center;gap:12px;margin:12px 0}
.event-controls[hidden]{display:none}
.event-controls label{display:flex;align-items:center;gap:6px}
.event-controls input,.event-controls select,.event-controls button,.event-table button{
 font:inherit;padding:6px 8px;border:1px solid #b6c8ce;border-radius:4px;background:white;color:#22313b}
.event-controls input[type=search]{width:clamp(180px,30vw,330px)}
.event-controls input[type=number]{width:6em}
button:not(:disabled),summary{cursor:pointer}button:disabled{opacity:.5}
button:focus-visible,input:focus-visible,select:focus-visible,a:focus-visible,summary:focus-visible{
 outline:3px solid #287f90;outline-offset:2px}
.event-table{max-height:65vh;border:1px solid #dde4e8}
.event-table table{width:100%}.event-table caption{text-align:left;padding:8px;font-weight:600}
.event-table thead{position:sticky;top:0;z-index:1}
.event-table th{white-space:nowrap}.event-table td{min-width:80px;max-width:420px;overflow-wrap:anywhere}
.event-table td.event-identity,.event-table td:nth-child(2){white-space:nowrap}
.event-table th[aria-sort=ascending] button::after{content:" ↑"}
.event-table th[aria-sort=descending] button::after{content:" ↓"}
.event-unexplained,.event-na{color:#667780}
.event-explanation{border:1px solid #cbdcde;border-radius:6px;padding:12px;scroll-margin-top:16px}
.event-explanation>summary{font-weight:600;color:#165d69}
.event-explanation>summary span{font-weight:400;color:#667780;margin-left:8px}
"""


BROWSER_JS = r"""
(() => {
  "use strict";
  for (const root of document.querySelectorAll('.event-browser')) {
  const prefix = root.dataset.prefix;
  const byId = id => document.getElementById(id.startsWith('event-') && !id.startsWith('event-explanation-')
    ? prefix + id.slice(5) : id);
  const data = JSON.parse(byId("event-data").textContent);
  const unit = data.row_unit || "事件";
  const body = byId("event-rows"), status = byId("event-status");
  const search = byId("event-search"), field = byId("event-search-field");
  const size = byId("event-page-size"), pageInput = byId("event-page");
  const previous = byId("event-prev"), next = byId("event-next");
  const headers = Array.from(root.querySelectorAll("[data-sort-column]"));
  const numeric = new Set(data.numeric);
  const state = {page: 0, size: data.page_size, sort: null, direction: 1};
  let visible = [];

  function compareRows(a, b) {
    const column = state.sort;
    let left = data.rows[a].values[column], right = data.rows[b].values[column];
    if (numeric.has(column)) {
      const number = value => value === null || value.trim() === "" || !Number.isFinite(Number(value))
        ? null : Number(value);
      left = number(left); right = number(right);
    }
    // NA is last in both directions. Ties retain original report order.
    if (left === null || right === null) {
      if (left === right) return a - b;
      return left === null ? 1 : -1;
    }
    return (left < right ? -state.direction : left > right ? state.direction : a - b);
  }

  function render() {
    const pages = Math.ceil(visible.length / state.size);
    state.page = Math.max(0, Math.min(state.page, pages - 1));
    const start = state.page * state.size, end = Math.min(start + state.size, visible.length);
    const fragment = document.createDocumentFragment();
    for (const index of visible.slice(start, end)) {
      const row = data.rows[index], tr = document.createElement("tr");
      tr.dataset.eventIndex = String(index);
      const ordinal = document.createElement("td"); ordinal.textContent = row.event_id || String(index + 1); tr.append(ordinal);
      const detail = document.createElement("td");
      if (row.details.length) {
        row.details.forEach((entryIndex, i) => {
          if (i) detail.append(document.createTextNode(" · "));
          const entry = data.details[entryIndex], link = document.createElement("a");
          link.href = entry.href || "#" + entry.id; link.dataset.detailLink = "";
          link.textContent = entry.href ? "查看" + (entry.scope || "通讯报告") : "查看解释（" + entry.scope + "）"; detail.append(link);
        });
      } else {
        detail.className = "event-unexplained"; detail.textContent = "尚未生成对应通讯报告";
        if (row.code) {
          const request = document.createElement("details"), summary = document.createElement("summary");
          summary.textContent = "追加代码";
          const code = document.createElement("pre"); code.textContent = row.code;
          request.append(summary, code); detail.append(request);
        }
      }
      if (data.show_details !== false) tr.append(detail);
      for (const [column, value] of row.values.entries()) {
        const cell = document.createElement("td");
        cell.textContent = value === null ? "NA" : value;
        if (["sender", "receiver", "hmdb_id", "sensor_gene"].includes(data.columns[column]))
          cell.classList.add("event-identity");
        if (value === null) cell.classList.add("event-na");
        tr.append(cell);
      }
      fragment.append(tr);
    }
    if (!visible.length) {
      const tr = document.createElement("tr"), cell = document.createElement("td");
      cell.colSpan = data.columns.length + (data.show_details === false ? 1 : 2);
      cell.textContent = data.rows.length ? "没有匹配的" + unit + "；可恢复初始视图。" : "无对应" + unit + "。";
      tr.append(cell); fragment.append(tr);
    }
    body.replaceChildren(fragment);
    root.querySelector(".event-table").scrollTop = 0;
    const order = state.sort === null ? "初始排序" : (data.labels || data.columns)[state.sort] +
      (state.direction === 1 ? " 升序" : " 降序") + "（NA 末尾）";
    status.textContent = "共 " + data.rows.length + " 条" + unit + "；匹配 " + visible.length +
      "；显示 " + (visible.length ? start + 1 : 0) + "–" + end + "；第 " +
      (pages ? state.page + 1 : 0) + " / " + pages + " 页；" + order + "。";
    previous.disabled = !pages || state.page === 0;
    next.disabled = !pages || state.page === pages - 1;
    pageInput.disabled = !pages; pageInput.max = String(Math.max(1, pages));
    pageInput.value = String(pages ? state.page + 1 : 1);
    for (const button of headers) {
      const active = Number(button.dataset.sortColumn) === state.sort;
      button.parentElement.setAttribute("aria-sort", active ?
        (state.direction === 1 ? "ascending" : "descending") : "none");
    }
  }

  function update() {
    const query = search.value.trim().toLowerCase();
    const columns = field.value === "" ? data.search_columns : [Number(field.value)];
    visible = [];
    data.rows.forEach((row, index) => {
      if (!query || columns.some(column => row.values[column] !== null &&
          row.values[column].toLowerCase().includes(query))) visible.push(index);
    });
    if (state.sort !== null) visible.sort(compareRows);
    state.page = 0; render();
  }
  search.addEventListener("input", update);
  field.addEventListener("change", update);
  size.addEventListener("change", () => {state.size = Number(size.value); state.page = 0; render();});
  previous.addEventListener("click", () => {state.page -= 1; render();});
  next.addEventListener("click", () => {state.page += 1; render();});
  pageInput.addEventListener("change", () => {
    const page = Number(pageInput.value);
    if (Number.isFinite(page) && Number.isInteger(page)) state.page = page - 1;
    render();
  });
  byId("event-reset").addEventListener("click", () => {
    search.value = ""; field.value = ""; size.value = String(data.page_size);
    state.size = data.page_size; state.sort = null; state.direction = 1; update();
  });
  document.addEventListener('cellmesh:focus-event', event => {
    if (!event.detail || event.detail.browser !== prefix) return;
    const index=data.rows.findIndex(row=>row.event_id===event.detail.eventId);
    if (index<0) return;
    search.value=''; field.value=''; state.sort=null; state.direction=1;
    update(); state.page=Math.floor(visible.indexOf(index)/state.size); render();
    for (let parent=root.parentElement;parent;parent=parent.parentElement) {
      if (parent.tagName==='DETAILS') parent.open=true;
    }
    const row=body.querySelector(`[data-event-index="${index}"]`);
    if (row) {row.tabIndex=-1; row.focus(); row.scrollIntoView({block:'center'});}
  });
  for (const button of headers) {
    button.disabled = false;
    button.addEventListener("click", () => {
      const column = Number(button.dataset.sortColumn);
      state.direction = state.sort === column ? -state.direction : 1;
      state.sort = column; update();
    });
  }
  const detailIds = new Set(data.details.map(entry => entry.id));
  function openDetail(id) {
    if (!detailIds.has(id)) return;
    const section = byId(id);
    if (!section) return;
    section.open = true;
    section.querySelector("summary").focus();
    section.scrollIntoView({block: "start"});
  }
  document.addEventListener("click", event => {
    const link = event.target.closest("a[data-detail-link]");
    if (link && link.getAttribute("href").startsWith("#")) openDetail(link.getAttribute("href").slice(1));
  });
  window.addEventListener("hashchange", () => openDetail(window.location.hash.slice(1)));
  update();
  byId("event-controls").hidden = false; byId("event-pager").hidden = false;
  openDetail(window.location.hash.slice(1));
  }
  function openEvidence(id) {
    const section = document.getElementById(id);
    if (!section || !section.classList.contains('expression-evidence')) return;
    section.open = true;
    section.scrollIntoView({block:'start'});
  }
  document.addEventListener('click', event => {
    const link = event.target.closest('a[data-evidence-link]');
    if (link) openEvidence(link.getAttribute('href').slice(1));
  });
  window.addEventListener('hashchange', () => openEvidence(window.location.hash.slice(1)));
  openEvidence(window.location.hash.slice(1));
})();
"""
