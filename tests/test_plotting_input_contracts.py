"""A06: input types must not change selection, missingness or reaction activity."""
from copy import deepcopy

import numpy as np
import pandas as pd
import pytest

matplotlib = pytest.importorskip('matplotlib')
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from cellmesh.plotting import (
    plot_significant_event_counts as counts,
    plot_communication_network as network,
    plot_event_dotplot as dots,
    plot_sample_event_scores as samples,
    plot_metabolite_secretion_violin as secretion,
    plot_receptor_expression_violin as receptor,
)
from test_plotting import single_cell_violin_inputs


@pytest.fixture(autouse=True)
def close_figures():
    yield
    plt.close('all')


@pytest.fixture
def events():
    return pd.DataFrame(dict(sender=['A'] * 3, receiver=['B'] * 3,
                             hmdb_id=['HMDB1', 'HMDB2', 'HMDB3'], sensor_gene=['R'] * 3,
                             cell_mesh_score=[0.4, 0.6, 0.8],
                             perm_pvalue=[0.01, 0.5, 0.9], fdr_sensor_type=[0.01, 0.5, 0.9]))


INVALID = [-0.1, np.inf, -np.inf, 1 + 2j, True, 'bad', 'NaN', '', pd.Timestamp('2020-01-01')]


@pytest.mark.parametrize('plot', [counts, network, dots])
@pytest.mark.parametrize('column', ['cell_mesh_score', 'perm_pvalue', 'fdr_sensor_type'])
@pytest.mark.parametrize('value', INVALID + [1.01])
def test_invalid_values_rejected_before_selection(events, plot, column, value):
    events[column] = events[column].astype(object)
    events.at[2, column] = value  # even when the offending row would lose top_n/thresholds
    kwargs = dict(max_perm_pvalue=0.1, max_fdr=0.1, min_cell_mesh_score=0.2)
    if plot is dots:
        kwargs['top_n'] = 1
    with pytest.raises(ValueError, match=column + '.*invalid records.*position 2'):
        plot(events, **kwargs)


@pytest.mark.parametrize('column', ['perm_pvalue', 'fdr_sensor_type'])
@pytest.mark.parametrize('plot', [network, dots])
def test_unused_threshold_statistics_still_validated_when_used_for_ranking(events, column, plot):
    events.loc[2, column] = -np.inf
    with pytest.raises(ValueError, match=column):
        plot(events, max_fdr=None, max_perm_pvalue=None)


def test_numeric_scientific_notation_rank_and_input_immutability(events):
    events['fdr_sensor_type'] = ['1e-2', '0.5', '0.9']
    events['cell_mesh_score'] = ['0.4', '0.6', '0.8']
    original = events.copy(deep=True)
    result = dots(events, top_n=1)
    assert result['selected_event_keys'] == [('HMDB1', 'R')]
    assert result['plot_events']['fdr_sensor_type'].tolist() == [0.01]
    pd.testing.assert_frame_equal(events, original)


@pytest.mark.parametrize('plot', [counts, network, dots])
def test_numeric_strings_match_numeric_filtering(events, plot):
    numeric = plot(events, max_fdr=0.1, max_perm_pvalue=0.1)
    strings = events.copy()
    for column in ['cell_mesh_score', 'perm_pvalue', 'fdr_sensor_type']:
        strings[column] = strings[column].astype(str)
    actual = plot(strings, max_fdr='1e-1', max_perm_pvalue='0.1')
    key = 'plot_events' if plot is dots else 'filtered_events'
    pd.testing.assert_frame_equal(actual[key], numeric[key])


@pytest.mark.parametrize('plot', [counts, network, dots])
def test_missing_score_excluded_with_reason_and_visible_count(events, plot):
    events.loc[1, 'cell_mesh_score'] = np.nan
    result = plot(events, max_fdr=None)
    assert result['missing_score_events']['hmdb_id'].tolist() == ['HMDB2']
    assert result['missing_score_events']['plot_exclusion_reason'].tolist() == ['missing_score']
    key = 'plot_events' if plot is dots else 'filtered_events'
    assert len(result[key]) == 2
    if plot is dots:
        result['fig'].canvas.draw()
        visible = sum(len(artist.get_offsets()) for axis in result['axes'] for artist in axis.collections)
        assert visible == result['thresholds']['events_plotted'] == 2
    elif plot is counts:
        assert result['counts'].to_numpy().sum() == 2
    else:
        assert result['edge_table']['event_count'].sum() == 2


@pytest.mark.parametrize('plot', [network, dots])
def test_normalized_selectors_and_context_scope(events, plot):
    events['sender'] = [' A ', ' A ', 'unused']
    events['receiver'] = ' B '
    events['cell_mesh_score'] = ['0.4', '0.6', 'bad']
    result = plot(events, sender_labels=[' A '], receiver_labels=[' B '], max_fdr=None)
    key = 'plot_events' if plot is dots else 'filtered_events'
    assert result[key]['sender'].astype(str).tolist() == ['A', 'A']
    assert result[key]['receiver'].astype(str).tolist() == ['B', 'B']


@pytest.mark.parametrize('plot', [network, dots])
@pytest.mark.parametrize('selector', ['sender_labels', 'receiver_labels'])
def test_selector_collisions_are_explicit(events, plot, selector):
    with pytest.raises(ValueError, match='duplicates after normalization'):
        plot(events, **{selector: ['A', ' A ']})


def test_explicit_dot_context_does_not_validate_unselected_event(events):
    events['cell_mesh_score'] = ['0.4', 'bad', '0.8']
    result = dots(events, event_keys=[('HMDB1', ' R ')])
    assert result['selected_event_keys'] == [('HMDB1', 'R')]


@pytest.mark.parametrize('plot', [counts, network, dots])
@pytest.mark.parametrize('parameter', ['max_fdr', 'max_perm_pvalue', 'min_cell_mesh_score'])
@pytest.mark.parametrize('value', [np.nan, np.inf, -0.1, 1.1, True, 1 + 0j])
def test_threshold_contracts(events, plot, parameter, value):
    with pytest.raises(ValueError, match=parameter):
        plot(events, **{parameter: value})


@pytest.mark.parametrize('bounds', [(30, 20), (-1, 20), (0, 20), (np.nan, 20), (20, np.inf), (True, 20)])
def test_dot_size_contract(events, bounds):
    with pytest.raises(ValueError, match='dot size range'):
        dots(events, min_dot_size=bounds[0], max_dot_size=bounds[1])


def test_fixed_positive_dot_size_allowed(events):
    result = dots(events, min_dot_size=25, max_dot_size=25)
    for artist in result['ax'].collections:
        assert np.all(artist.get_sizes() == 25)


@pytest.mark.parametrize('parameter', ['node_size_range', 'edge_width_range', 'node_value_range',
                                      'edge_width_value_range', 'edge_color_value_range'])
@pytest.mark.parametrize('bounds', [(np.nan, 2), (1, np.inf), (2, 1), (True, 2), (1 + 0j, 2)])
def test_network_ranges(events, parameter, bounds):
    with pytest.raises(ValueError, match=parameter):
        network(events, **{parameter: bounds})


@pytest.mark.parametrize('parameter', ['curve', 'edge_alpha', 'node_label_size'])
@pytest.mark.parametrize('value', [np.nan, np.inf, -1, True])
def test_network_numeric_parameters(events, parameter, value):
    with pytest.raises(ValueError, match=parameter):
        network(events, **{parameter: value})


@pytest.fixture
def sample_events(events):
    out = pd.concat([events.iloc[[0]]] * 3, ignore_index=True)
    out['sample'] = [' S1 ', ' S2 ', ' S3 ']
    out['cell_mesh_score'] = [0.0, 0.8, np.nan]
    return out


def sample_plot(events, **kwargs):
    return samples(events, sender=' A ', receiver=' B ', hmdb_id='HMDB1', sensor_gene=' R ', **kwargs)


@pytest.mark.parametrize('value', INVALID + [1.01])
def test_sample_bad_value_cannot_become_missing(sample_events, value):
    sample_events['cell_mesh_score'] = sample_events['cell_mesh_score'].astype(object)
    sample_events.at[1, 'cell_mesh_score'] = value
    with pytest.raises(ValueError, match='cell_mesh_score.*invalid records'):
        sample_plot(sample_events)


def test_sample_true_missing_zero_order_and_typed_aliases(sample_events):
    alias = sample_events.iloc[[1]].copy()
    alias['cell_mesh_score'] = '8e-1'
    result = sample_plot(pd.concat([sample_events, alias]), sample_order=[' S3 ', ' S1 ', ' S2 '])
    assert result['sample_scores']['sample'].astype(str).tolist() == ['S3', 'S1', 'S2']
    assert result['missing_samples'] == ['S3']
    assert result['median_score'] == 0.4
    assert result['n_samples_coobserved'] == 2
    assert [bar.get_height() for bar in result['ax'].patches] == [0.0, 0.8]


def test_sample_order_duplicates_after_normalization(sample_events):
    with pytest.raises(ValueError, match='duplicates after normalization'):
        sample_plot(sample_events, sample_order=['S1', ' S1 '])


@pytest.mark.parametrize('plot', [network, dots, sample_plot])
def test_custom_scores_can_exceed_one(events, sample_events, plot):
    data = sample_events if plot is sample_plot else events
    data = data.assign(custom_score=['2e0', '4', '6'])
    result = plot(data, score_col='custom_score')
    key = 'sample_scores' if plot is sample_plot else 'plot_events' if plot is dots else 'filtered_events'
    assert result[key]['custom_score'].min() >= 2


@pytest.mark.parametrize('plot,column', [(secretion, 'sender'), (receptor, 'sensor_score'),
                                       (receptor, 'sensor_expr_frac')])
@pytest.mark.parametrize('value', INVALID + [1.01])
def test_violin_metadata_contract(single_cell_violin_inputs, plot, column, value):
    result, adata = single_cell_violin_inputs
    if column == 'sender':
        result.sender_scores['A'] = pd.Series([value], index=result.sender_scores.index, dtype=object)
    else:
        result.receiver_scores[column] = result.receiver_scores[column].astype(object)
        result.receiver_scores.loc[0, column] = value
    kwargs = {'receptor_gene': 'G1'} if plot is receptor else {}
    with pytest.raises(ValueError, match='invalid records'):
        plot(result, adata, hmdb_id='HMDB1', **kwargs)


@pytest.mark.parametrize('plot', [secretion, receptor])
@pytest.mark.parametrize('kwargs', [dict(vmin=np.nan), dict(vmax=np.inf), dict(vmin=2, vmax=1),
                                   dict(vmin=True), dict(vmax=1 + 2j)])
def test_violin_color_ranges(single_cell_violin_inputs, plot, kwargs):
    result, adata = single_cell_violin_inputs
    if plot is receptor:
        kwargs = dict(kwargs, receptor_gene='G1')
    with pytest.raises(ValueError, match='vmin|vmax'):
        plot(result, adata, hmdb_id='HMDB1', **kwargs)


@pytest.mark.parametrize('variant', [[' P1 ', 'P2'], ['P1', 'P1', 'P2'], [' P2 ', 'P1', ' P1 ']])
def test_reaction_canonical_list_used_for_calculation(single_cell_violin_inputs, variant):
    result, adata = single_cell_violin_inputs
    adata.X[:, 1] = [9, 7, 5, 3, 2, 1]  # unequal genes reveal unintended duplicate weights
    expected = secretion(result, adata, hmdb_id='HMDB1')
    assert expected['plot_data']['production_score'].iloc[0] == pytest.approx(np.sqrt(50) - 1)
    changed = deepcopy(result)
    changed.availability_results['reaction_genes'].at[0, 'genes'] = variant
    original = deepcopy(changed.availability_results['reaction_genes'])
    actual = secretion(changed, adata, hmdb_id='HMDB1', sender_labels=[' A ', ' B '])
    pd.testing.assert_frame_equal(actual['plot_data'], expected['plot_data'])
    assert actual['reaction_definitions'].iloc[0]['genes'] == ['P1', 'P2']
    assert actual['reaction_metadata'].iloc[0]['genes'] == ['P1', 'P2']
    pd.testing.assert_frame_equal(changed.availability_results['reaction_genes'], original)


@pytest.mark.parametrize('gene', [' ', None, np.nan])
def test_missing_reaction_gene_identifier_rejected(single_cell_violin_inputs, gene):
    result, adata = single_cell_violin_inputs
    result.availability_results['reaction_genes'].at[0, 'genes'] = [gene, 'P2']
    with pytest.raises(ValueError, match='reaction_genes.genes'):
        secretion(result, adata, hmdb_id='HMDB1')


def test_receptor_numeric_strings_and_selector_equivalence(single_cell_violin_inputs):
    result, adata = single_cell_violin_inputs
    expected = receptor(result, adata, hmdb_id='HMDB1', receptor_gene='G1')
    for column in ['sensor_score', 'sensor_expr_frac']:
        result.receiver_scores[column] = result.receiver_scores[column].astype(str)
    actual = receptor(result, adata, hmdb_id='HMDB1', receptor_gene=' G1 ', receiver_labels=[' A ', ' B '])
    pd.testing.assert_frame_equal(actual['summary'], expected['summary'])
    assert actual['summary']['mean_receptor_expression'].max() > 1  # raw expression is unbounded


@pytest.mark.parametrize('plot', [secretion, receptor])
def test_violin_missing_scores_reported_and_not_filled(single_cell_violin_inputs, plot):
    result, adata = single_cell_violin_inputs
    kwargs = dict(hmdb_id='HMDB1')
    if plot is secretion:
        result.sender_scores.loc[:, 'A'] = np.nan
    else:
        result.receiver_scores.loc[0, 'sensor_score'] = np.nan
        kwargs['receptor_gene'] = 'G1'
    plotted = plot(result, adata, **kwargs)
    assert plotted['missing_score_cell_types'] == ['A']
    assert plotted['cell_order'] == ['B']
    assert len(plotted['plot_data']) == 3


@pytest.mark.parametrize('plot', [secretion, receptor])
def test_violin_unselected_context_not_validated(single_cell_violin_inputs, plot):
    result, adata = single_cell_violin_inputs
    kwargs = dict(hmdb_id='HMDB1')
    if plot is secretion:
        result.sender_scores['A'] = 'bad'
        for key in ('P', 'C', 'E'):
            result.availability_results[key]['A'] = 'bad'
        kwargs['sender_labels'] = [' B ']
    else:
        result.receiver_scores['sensor_score'] = ['bad', '0.75']
        result.receiver_scores['sensor_expr_frac'] = ['bad', '1.0']
        kwargs.update(receiver_labels=[' B '], receptor_gene='G1')
    plotted = plot(result, adata, **kwargs)
    assert plotted['cell_order'] == ['B']


@pytest.mark.parametrize('first_variant', [True, False])
def test_equivalent_reaction_duplicates_count_once(single_cell_violin_inputs, first_variant):
    result, adata = single_cell_violin_inputs
    adata.X[:, 1] = [9, 7, 5, 3, 2, 1]
    expected = secretion(result, adata, hmdb_id='HMDB1')
    reactions = result.availability_results['reaction_genes']
    alias = reactions.iloc[[0]].copy(deep=True)
    alias.at[0, 'genes'] = [' P2 ', 'P1', 'P1']
    result.availability_results['reaction_genes'] = pd.concat(
        [alias, reactions] if first_variant else [reactions, alias], ignore_index=True,
    )
    actual = secretion(result, adata, hmdb_id='HMDB1')
    pd.testing.assert_frame_equal(actual['plot_data'], expected['plot_data'])
    assert len(actual['reaction_definitions']) == 3


def test_sender_numeric_alias_rows_and_intermediate_labels(single_cell_violin_inputs):
    result, adata = single_cell_violin_inputs
    expected = secretion(result, adata, hmdb_id='HMDB1')
    result.sender_scores = pd.concat([result.sender_scores, result.sender_scores.astype(str)])
    result.sender_scores.columns = [' A ', ' B ']
    for key in ('P', 'C', 'E'):
        result.availability_results[key] = result.availability_results[key].astype(str)
        result.availability_results[key].columns = [' A ', ' B ']
    actual = secretion(result, adata, hmdb_id='HMDB1')
    pd.testing.assert_frame_equal(actual['summary'], expected['summary'])


@pytest.mark.parametrize('plot', [counts, network, dots])
def test_qc_excluded_bad_values_do_not_block_selected_plot(events, plot):
    events['passes_min_cells'] = [True, False, True]
    events['cell_mesh_score'] = ['0.4', 'bad', '0.8']
    result = plot(events, max_fdr=None)
    assert len(result['qc_excluded_events']) == 1
    with pytest.raises(ValueError, match='cell_mesh_score'):
        plot(events, max_fdr=None, qc_only=False)


@pytest.mark.parametrize('plot', [counts, network, dots, sample_plot])
@pytest.mark.parametrize('native_column', [np.array([1+2j, 0+0j, 1+0j]), np.array([True, False, True])])
def test_invalid_native_dtype_raises_clear_error_without_cast_warnings(events, sample_events, plot, native_column):
    import warnings
    data = (sample_events if plot is sample_plot else events).copy()
    data['cell_mesh_score'] = native_column
    with warnings.catch_warnings():
        warnings.simplefilter('error')
        with pytest.raises(ValueError, match='cell_mesh_score.*invalid records'):
            plot(data)


@pytest.mark.parametrize('plot', [secretion, receptor])
def test_violin_complex_metadata_no_silent_cast(single_cell_violin_inputs, plot):
    import warnings
    result, adata = single_cell_violin_inputs
    kwargs = dict(hmdb_id='HMDB1')
    if plot is secretion:
        result.sender_scores['A'] = 1 + 2j
    else:
        result.receiver_scores['sensor_score'] = [1 + 2j, 0 + 0j]
        kwargs['receptor_gene'] = 'G1'
    with warnings.catch_warnings():
        warnings.simplefilter('error')
        with pytest.raises(ValueError, match='invalid records'):
            plot(result, adata, **kwargs)
