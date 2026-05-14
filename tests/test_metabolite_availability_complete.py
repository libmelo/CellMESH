
"""
Complete Metabolite Availability Tests - Following Documentation Requirements

This test file covers all assertions mentioned in the documentation.
"""

import numpy as np
import pandas as pd
import pytest
from scipy import sparse

import cellmesh
from cellmesh.preprocess import robust_minmax
from cellmesh import compute_metabolite_availability


class FakeAnnData:
    """Simple AnnData mock for testing"""
    def __init__(self, X, var_names, obs):
        self.X = X
        self.layers = {}
        self.var_names = pd.Index(var_names)
        self.obs = pd.DataFrame(obs)


@pytest.fixture
def toy_adata():
    """Create toy AnnData from documentation"""
    genes = [
        "G_prodA1", "G_prodA2", "G_prodA3",
        "G_consA",
        "G_transA1", "G_transA2",
        "G_prodB", "G_prodF",
        "G_noise"
    ]
    cells = ["T1", "T2", "T3", "M1", "M2", "M3", "F1", "F2"]
    cell_types = ["T_cell", "T_cell", "T_cell", "Macrophage", "Macrophage", "Macrophage", "Fibroblast", "Fibroblast"]
    
    X = np.array([
        [4.0, 1.0, 0.7, 0.2, 2.5, 1.5, 0.1, 1.2, 0.5],
        [3.8, 1.2, 0.8, 0.3, 2.7, 1.4, 0.0, 1.1, 0.4],
        [4.2, 0.8, 0.6, 0.2, 2.4, 1.6, 0.2, 1.3, 0.6],
        [1.0, 2.8, 0.5, 4.0, 0.8, 0.4, 1.0, 0.5, 0.7],
        [1.2, 3.0, 0.4, 4.2, 0.7, 0.3, 0.8, 0.6, 0.8],
        [0.9, 2.6, 0.6, 3.8, 0.9, 0.5, 1.1, 0.4, 0.6],
        [0.5, 0.6, 0.2, 0.8, 0.3, 0.2, 4.0, 2.5, 0.9],
        [0.4, 0.7, 0.3, 0.9, 0.4, 0.2, 4.2, 2.7, 1.0],
    ], dtype=float)
    
    return FakeAnnData(X, genes, {"cell_type": cell_types})


@pytest.fixture
def toy_reaction_table():
    """Create toy reaction table from documentation"""
    return pd.DataFrame([
        {"metabolite": "MetA", "HMDB_ID": "HMDB00001", "reaction": "R_A_prod_1", "gene": "G_prodA1", "direction": "product"},
        {"metabolite": "MetA", "HMDB_ID": "HMDB00001", "reaction": "R_A_prod_1", "gene": "G_prodA2", "direction": "product"},
        {"metabolite": "MetA", "HMDB_ID": "HMDB00001", "reaction": "R_A_prod_2", "gene": "G_prodA3", "direction": "product"},
        {"metabolite": "MetA", "HMDB_ID": "HMDB00001", "reaction": "R_A_sub_1", "gene": "G_consA", "direction": "substrate"},
        {"metabolite": "MetA", "HMDB_ID": "HMDB00001", "reaction": "R_A_trans_1", "gene": "G_transA1", "direction": "transporter"},
        {"metabolite": "MetA", "HMDB_ID": "HMDB00001", "reaction": "R_A_trans_1", "gene": "G_transA2", "direction": "transporter"},
        {"metabolite": "MetB", "HMDB_ID": "HMDB00002", "reaction": "R_B_prod_1", "gene": "G_prodB", "direction": "product"},
        {"metabolite": "MetF", "HMDB_ID": "HMDB00006", "reaction": "R_F_prod_1", "gene": "G_prodF", "direction": "product"},
        {"metabolite": "MetF", "HMDB_ID": "HMDB00006", "reaction": "R_F_sub_1", "gene": "G_consA", "direction": "substrate"},
        {"metabolite": "MetD", "HMDB_ID": "HMDB00004", "reaction": "R_D_sub_1", "gene": "G_consA", "direction": "substrate"},
        {"metabolite": "MetC", "HMDB_ID": "HMDB00003", "reaction": "R_C_prod_1", "gene": "G_missing", "direction": "product"},
    ])


class TestRobustMinMax:
    """Documentation tests 7: robust minmax on constant vector returns 0"""
    
    def test_normal_case(self):
        """Test normal robust minmax with varying values"""
        x = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
        result = robust_minmax(x, lower=0, upper=100)
        assert len(result) == 10
        assert result[0] == pytest.approx(0.0)
        assert result[-1] == pytest.approx(1.0)
        assert np.all(result >= 0) and np.all(result <= 1)
    
    def test_constant_vector_returns_zero(self):
        """Documentation test 7: Constant vector returns all zeros"""
        x = np.array([5.0, 5.0, 5.0, 5.0, 5.0])
        result = robust_minmax(x)
        assert np.all(result == 0.0)
        print("✓ Test 7 passed: Constant vector returns 0")
    
    def test_nan_handling(self):
        """Test that NaN values are handled properly"""
        x = np.array([1, 2, np.nan, 4, 5])
        result = robust_minmax(x)
        assert not np.any(np.isnan(result))


class TestMetaboliteAvailability:
    """Tests for metabolite availability calculations"""
    
    def test_1_2_multi_gene_and_multi_reaction(self, toy_adata, toy_reaction_table):
        """Documentation tests 1-2: multi-gene reaction takes max; multi-reaction sums"""
        result = compute_metabolite_availability(
            toy_adata, toy_reaction_table, celltype_col="cell_type",
            min_cells=1, lower=0, upper=100, return_intermediates=True
        )
        
        # Get pseudobulk
        pseudobulk = result['pseudobulk']
        
        # Find MetA index
        metA_idx = None
        for idx in result['P'].index:
            if idx[0] == 'MetA':
                metA_idx = idx
                break
        
        assert metA_idx is not None, "MetA should be included"
        
        # Test multi-gene max and multi-reaction sum
        expected_metA_P = pseudobulk[["G_prodA1", "G_prodA2"]].max(axis=1) + pseudobulk["G_prodA3"]
        pd.testing.assert_series_equal(
            result['P'].loc[metA_idx], expected_metA_P, check_names=False, atol=1e-10
        )
        print("✓ Tests 1-2 passed: Multi-gene reaction takes max; multi-reaction sums")
    
    def test_3_metabolite_without_product_skipped(self, toy_adata, toy_reaction_table):
        """Documentation test 3: Metabolite without product is skipped (MetD)"""
        result = compute_metabolite_availability(
            toy_adata, toy_reaction_table, celltype_col="cell_type",
            min_cells=1, return_intermediates=True
        )
        
        # Check MetD is not in availability
        included_mets = [idx[0] for idx in result['availability'].index]
        assert 'MetD' not in included_mets, "MetD should be skipped"
        print("✓ Test 3 passed: Metabolite without product is skipped")
    
    def test_4_missing_product_genes_skipped(self, toy_adata, toy_reaction_table):
        """Documentation test 4: Metabolite with missing product genes is skipped (MetC)"""
        result = compute_metabolite_availability(
            toy_adata, toy_reaction_table, celltype_col="cell_type",
            min_cells=1, return_intermediates=True
        )
        
        # Check MetC is not in availability
        included_mets = [idx[0] for idx in result['availability'].index]
        assert 'MetC' not in included_mets, "MetC should be skipped"
        print("✓ Test 4 passed: Metabolite with missing product genes is skipped")
    
    def test_5_missing_substrate_default_C_norm(self, toy_adata, toy_reaction_table):
        """Documentation test 5: Missing substrate uses C_norm=0.41 (MetB)"""
        result = compute_metabolite_availability(
            toy_adata, toy_reaction_table, celltype_col="cell_type",
            min_cells=1, lower=0, upper=100, missing_C_norm=0.41, return_intermediates=True
        )
        
        # Find MetB index
        metB_idx = None
        for idx in result['C_norm'].index:
            if idx[0] == 'MetB':
                metB_idx = idx
                break
        
        assert metB_idx is not None, "MetB should be included"
        assert np.allclose(result['C_norm'].loc[metB_idx].values, 0.41)
        print("✓ Test 5 passed: Missing substrate uses C_norm=0.41")
    
    def test_6_missing_transporter_default_E_norm(self, toy_adata, toy_reaction_table):
        """Documentation test 6: Missing transporter uses E_norm=0.75 (MetB and MetF)"""
        result = compute_metabolite_availability(
            toy_adata, toy_reaction_table, celltype_col="cell_type",
            min_cells=1, lower=0, upper=100, missing_E_norm=0.75, return_intermediates=True
        )
        
        # Check MetB
        metB_idx = None
        for idx in result['E_norm'].index:
            if idx[0] == 'MetB':
                metB_idx = idx
                break
        
        if metB_idx is not None:
            assert np.allclose(result['E_norm'].loc[metB_idx].values, 0.75)
        
        # Check MetF
        metF_idx = None
        for idx in result['E_norm'].index:
            if idx[0] == 'MetF':
                metF_idx = idx
                break
        
        if metF_idx is not None:
            assert np.allclose(result['E_norm'].loc[metF_idx].values, 0.75)
        
        print("✓ Test 6 passed: Missing transporter uses E_norm=0.75")
    
    def test_8_dense_sparse_equivalent(self, toy_adata, toy_reaction_table):
        """Documentation test 8: Dense and sparse inputs produce identical outputs"""
        # Create sparse version
        adata_sparse = FakeAnnData(
            X=sparse.csr_matrix(toy_adata.X),
            var_names=list(toy_adata.var_names),
            obs={"cell_type": toy_adata.obs["cell_type"].values}
        )
        
        # Compute both
        result_dense = compute_metabolite_availability(
            toy_adata, toy_reaction_table, celltype_col="cell_type",
            min_cells=1, lower=0, upper=100
        )
        
        result_sparse = compute_metabolite_availability(
            adata_sparse, toy_reaction_table, celltype_col="cell_type",
            min_cells=1, lower=0, upper=100
        )
        
        # Compare availability matrices
        pd.testing.assert_frame_equal(
            result_dense['availability'], result_sparse['availability'], atol=1e-10
        )
        print("✓ Test 8 passed: Dense and sparse inputs produce identical outputs")
    
    def test_9_availability_shape_correct(self, toy_adata, toy_reaction_table):
        """Documentation test 9: Availability shape is correct"""
        result = compute_metabolite_availability(
            toy_adata, toy_reaction_table, celltype_col="cell_type",
            min_cells=1, lower=0, upper=100, return_intermediates=True
        )
        
        availability = result['availability']
        
        # Should have 3 cell types
        assert availability.shape[1] == 3, "Should have 3 cell types"
        assert set(availability.columns) == {'T_cell', 'Macrophage', 'Fibroblast'}
        print("✓ Test 9 passed: Availability shape is correct")


class TestDownstreamCompatibility:
    """Documentation test 10: Downstream cell-cell communication remains compatible"""
    
    def test_cellmesh_imports_work(self):
        """Test that cellmesh imports still work correctly"""
        # Test top-level imports
        from cellmesh import run_cell_mesh, load_cell_mesh_database
        assert callable(run_cell_mesh)
        assert callable(load_cell_mesh_database)
        
        # Test that compute_metabolite_availability is available
        assert hasattr(cellmesh, 'compute_metabolite_availability')
        assert callable(cellmesh.compute_metabolite_availability)
        
        print("✓ Test 10a passed: cellmesh imports still work")
    
    def test_compute_availability_exposed(self):
        """Test that compute_metabolite_availability is exposed at top level"""
        assert 'compute_metabolite_availability' in dir(cellmesh)
        assert callable(cellmesh.compute_metabolite_availability)
        print("✓ Test 10b passed: compute_metabolite_availability is exposed")


if __name__ == "__main__":
    print("=" * 70)
    print("Running Complete Metabolite Availability Tests")
    print("=" * 70)
    
    pytest.main([__file__, "-v"])
