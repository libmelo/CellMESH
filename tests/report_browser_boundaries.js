(() => {
  const get = id => document.getElementById(id);
  const check = (ok, message) => {if (!ok) throw new Error(message);};
  if (get('event-controls').hidden) {
    check(document.querySelectorAll('#event-rows tr[data-event-index]').length === 50,
          'without JavaScript, first 50 rows remain readable');
    check(get('event-status').textContent.includes('静态预览前 50 条'), 'static fallback is explicitly labeled');
    check(document.querySelectorAll('.event-explanation').length === 2, 'explanations remain present without JavaScript');
    return JSON.stringify({passed: 3, mode: 'no-javascript fallback'});
  }
  if (window.location.hash) {
    const id = window.location.hash.slice(1);
    check(get(id).open, 'fragment on initial load opens the corresponding explanation');
    check(document.querySelectorAll('#event-rows tr[data-event-index]').length === 50, 'fragment does not change event scope');
    return JSON.stringify({passed: 2, mode: 'initial fragment'});
  }
  const data = JSON.parse(get('event-data').textContent);
  if (data.rows.length === 10001) {
    const started = performance.now();
    const rows = () => Array.from(document.querySelectorAll('#event-rows tr[data-event-index]'));
    check(rows().length === 50, 'large table renders only 50 rows');
    get('event-page').value = '201'; get('event-page').dispatchEvent(new Event('change'));
    check(rows().length === 1 && rows()[0].dataset.eventIndex === '10000', 'large table last row is reachable');
    get('event-search').value = 'HMDB0010001'; get('event-search').dispatchEvent(new Event('input'));
    check(rows().length === 1 && rows()[0].dataset.eventIndex === '10000', 'large table search covers last row');
    get('event-reset').click();
    check(rows().length === 50, 'large-table reset restores bounded DOM');
    return JSON.stringify({passed: 4, mode: '10001-row table', interaction_ms: performance.now() - started});
  }
  check(data.rows.length === 0 && data.details.length === 0, 'empty selection remains empty');
  check(get('event-next').disabled && get('event-prev').disabled && get('event-page').disabled,
        'empty report disables pagination');
  check(get('event-status').textContent.includes('第 0 / 0 页'), 'empty-page indicator');
  get('event-search').value = 'anything'; get('event-search').dispatchEvent(new Event('input'));
  document.querySelector('[data-sort-column="0"]').click();
  get('event-reset').click();
  check(get('event-status').textContent.includes('共 0 条事件'), 'empty search/sort/reset succeeds');
  check(!document.querySelector('.event-explanation'), 'empty report does not generate explanation sections');
  return JSON.stringify({passed: 5, mode: 'empty report'});
})();
