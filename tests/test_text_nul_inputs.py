"""B04: validate decoded text before parsers can truncate NUL fields."""
import gzip

import numpy as np
import pandas as pd
import pytest
from scipy.io import mmwrite
from scipy import sparse

from cellmesh import read_anndata


@pytest.mark.parametrize('mode,sep', [('csv', ','), ('tsv', '\t')])
@pytest.mark.parametrize('engine', ['c', 'python'])
@pytest.mark.parametrize('encoding', ['utf-8', 'utf-16'])
@pytest.mark.parametrize('compressed', [False, True])
@pytest.mark.parametrize('field,line', [('header', 1), ('cell', 2), ('value', 2)])
def test_nul_expression_rejected_before_parsing(tmp_path, mode, sep, engine, encoding, compressed, field, line):
    text = 'cell,G,R\na,1,2\nb,3,4\n'
    old, new = {'header': ('G,R', 'G\0suffix,R'), 'cell': ('a,1', 'a\0suffix,1'),
                'value': ('a,1', 'a,1\0junk')}[field]
    payload = text.replace(old, new).replace(',', sep).encode(encoding)
    path = tmp_path / ('input.' + mode + ('.gz' if compressed else ''))
    path.write_bytes(gzip.compress(payload) if compressed else payload)
    with pytest.raises(ValueError, match=rf'input.*line {line}.*NUL'):
        read_anndata(path, mode=mode, engine=engine, encoding=encoding)


@pytest.mark.parametrize('kind', ['cell_meta', 'gene_meta', 'genes', 'barcodes'])
@pytest.mark.parametrize('compressed', [False, True])
def test_nul_auxiliary_files_rejected_by_public_reader(tmp_path, kind, compressed):
    expression = tmp_path / 'expression.csv'
    expression.write_text('cell,G,R\na,1,2\nb,3,4\n')
    matrix = tmp_path / 'matrix.mtx'
    mmwrite(matrix, sparse.coo_matrix(np.array([[1., 3.], [2., 4.]])))
    content = {'cell_meta': 'cell,cell_type\na,A\0suffix\nb,B\n',
               'gene_meta': 'gene,note\nG,bad\0suffix\nR,ok\n',
               'genes': 'G\0suffix\nR\n', 'barcodes': 'a\0suffix\nb\n'}[kind]
    path = tmp_path / ('bad.csv' if 'meta' in kind else 'bad.tsv')
    if compressed:
        path = path.with_suffix(path.suffix + '.gz')
    path.write_bytes(gzip.compress(content.encode()) if compressed else content.encode())
    source, mode = (expression, 'csv') if 'meta' in kind else (matrix, 'mtx')
    with pytest.raises(ValueError, match=r'bad.*line [12].*NUL'):
        read_anndata(source, mode=mode, **{kind + '_path': path})


@pytest.mark.parametrize('encoding', ['utf-8-sig', 'utf-16', 'utf-32'])
@pytest.mark.parametrize('compressed', [False, True])
def test_decoded_scan_preserves_quotes_newlines_and_literal_escape(tmp_path, encoding, compressed):
    frame = pd.DataFrame([[1, 2], [3, 4]], index=['a\nb', r'c\0d'], columns=['G,1', '基因'])
    path = tmp_path / ('input.csv.gz' if compressed else 'input.csv')
    frame.to_csv(path, encoding=encoding)
    actual = read_anndata(path, mode='csv', encoding=encoding)
    assert actual.obs_names.tolist() == frame.index.tolist()
    assert actual.var_names.tolist() == frame.columns.tolist()
    np.testing.assert_array_equal(actual.X, frame.to_numpy())


def test_nul_physical_line_inside_quoted_record(tmp_path):
    path = tmp_path / 'input.csv'
    path.write_text('cell,G\n"a\nb\0c",1\n')
    with pytest.raises(ValueError, match=r'line 3.*NUL'):
        read_anndata(path, mode='csv')


@pytest.mark.parametrize('options', [{'nrows': 1}, {'usecols': ['cell', 'G']}, {'skiprows': [2]}])
def test_parser_selection_does_not_hide_corrupt_text(tmp_path, options):
    path = tmp_path / 'input.csv'
    path.write_text('cell,G,R\na,1,2\nb,3,4\0junk\n')
    with pytest.raises(ValueError, match=r'line 3.*NUL'):
        read_anndata(path, mode='csv', **options)
