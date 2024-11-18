import anndata
import pytest

from benccchmarker.core.data_generator import DataGenerator

@pytest.mark.parametrize(
        "num_cells, num_genes, num_lr_pairs, num_cell_types",
        [
            (500, 100, 10, 5),
            (1000, 100, 10, 4),
            (5000, 1000, 10, 3),
        ]
)
def test_denovo_simulation_output_shape(num_cells, num_genes, num_lr_pairs, num_cell_types):
    output_path = "benccchmarker/tmp/test_output"
    data_generator = DataGenerator()

    data_generator.simulate_denovo(
        num_cells=num_cells,
        num_genes=num_genes,
        num_lr_pairs=num_lr_pairs,
        num_cell_types=num_cell_types,
        output_path=output_path
    )

    simulated_data_path = f"{output_path}/simulated_data.h5ad"

    adata = anndata.read_h5ad(simulated_data_path)

    assert adata.shape == (num_cells, num_genes)

    # Check the number of cell types in adata
    assert adata.obs["cell_type"].nunique() == num_cell_types

@pytest.mark.parametrize(
        "num_genes, num_lr_pairs",
        [
            (100, 200),
            (100, 300),
        ]
)
def test_denovo_simulation_input_size(num_genes, num_lr_pairs):
    output_path = "benccchmarker/tmp/test_output"
    data_generator = DataGenerator()

    with pytest.raises(ValueError):
        data_generator.simulate_denovo(
            num_cells=500,
            num_genes=num_genes,
            num_lr_pairs=num_lr_pairs,
            output_path=output_path
        )

@pytest.mark.parametrize(
        "mean_expression, dispersion",
        [
            (100, 5),
            (200, 10), #
        ]
)
def test_invalid_distribution_parameters(mean_expression, dispersion):
    data_generator = DataGenerator()
    
    # Test invalid parameters
    with pytest.raises(ValueError, match="Invalid mean expression and dispersion parameters"):
        data_generator._check_distribution_parameters(mean_expression, dispersion)

@pytest.mark.parametrize(
        "num_cells, num_genes",
        [
            (1000, 100),
            (500, 200),
        ]
)
def test_simulation_output_consistency(num_cells, num_genes):
    output_path = "benccchmarker/tmp/test_output"
    data_generator = DataGenerator()

    data_generator.simulate_denovo(
        num_cells=num_cells,
        num_genes=num_genes,
        num_lr_pairs=10,
        output_path=output_path
    )

    simulated_data_path = f"{output_path}/simulated_data.h5ad"

    adata = anndata.read_h5ad(simulated_data_path)

    test_adata = anndata.read_h5ad(f"benccchmarker/test/data/simulated_data_{num_cells}c_{num_genes}g.h5ad")

    # Check if the generated data is consistent with the test data
    assert adata.X.shape == test_adata.X.shape
    assert adata.var_names.tolist() == test_adata.var_names.tolist()
    assert adata.obs_names.tolist() == test_adata.obs_names.tolist()

    # Check if some random numbers in the generated data is consistent with the test data
    random_indices = [[num_cells-num_cells, num_genes-num_genes], [num_cells - int(num_cells / 2), num_genes - int(num_genes/2)], [num_cells - int(num_cells / 3), num_genes - int(num_genes/3)], [num_cells-1, num_genes-1]]
    for random_index in random_indices:
        assert adata.X[random_index[0], random_index[1]] == test_adata.X[random_index[0], random_index[1]]

    # Load and compare if the ground_truth.csv and test_ground_truth.csv are the same

    with open(f"{output_path}/simulated_data_ground_truth.csv", "r") as f:
        ground_truth = f.read()
    
    test_ground_truth_path = f"simulated_data_ground_truth_{num_cells}c_{num_genes}g.csv"

    with open(f"benccchmarker/test/data/{test_ground_truth_path}", "r") as f:
        test_ground_truth = f.read()
    
    assert ground_truth == test_ground_truth

@pytest.mark.parametrize(
        "num_cells, num_genes, num_lr_pairs, overexpression_scale, target_gene, target_cell, gene_to_compare, cell_to_compare, control_gene, expected_output",
        [
            # The number is obtained from me checking the data prior to the test
            # Here I am expecting that ACKR4 in CellType2 to be 2 times higher than 
            # the mean expression value of CX3CR1 in the same type of cell because 
            # of the injected 'fake' communication, same goes for the CellType3
            # as CX3CR1 is overexpressed by 2 times in CellType3
            (1000, 100, 10, 2, "CX3CR1", "CellType2", "ACKR4", "CellType3", "CCDC167", {"mean": 99.21993127147766, "max": 198, "min": 37}),
        ]
)
def test_expression_values(num_cells, num_genes, num_lr_pairs, overexpression_scale, target_gene, target_cell, gene_to_compare, cell_to_compare, control_gene, expected_output):
    output_path = "benccchmarker/tmp/test_output"
    data_generator = DataGenerator()

    data_generator.simulate_denovo(
        num_cells=num_cells,
        num_genes=num_genes,
        num_lr_pairs=num_lr_pairs,
        overexpression_scale=overexpression_scale,
        output_path=output_path
    )

    simulated_data_path = f"{output_path}/simulated_data.h5ad"

    adata = anndata.read_h5ad(simulated_data_path)

    # Check if the mean, min and max expression value for target_gene in target_cell 
    # is approximately equal to the expected output
    target_gene_expression_values = adata.X[adata.obs["cell_type"] == target_cell, adata.var_names == target_gene].flatten()

    assert target_gene_expression_values.mean() == pytest.approx(expected_output["mean"], rel=1e-2)
    assert target_gene_expression_values.max() == expected_output["max"]
    assert target_gene_expression_values.min() == expected_output["min"]

    # Check if the mean expression value between target_gene and gene_to_compare 
    # is approximately equal to overexpression_scale
    gene_to_compare_expression_values = adata.X[adata.obs["cell_type"] == target_cell, adata.var_names == gene_to_compare].flatten()

    assert  gene_to_compare_expression_values.mean() / target_gene_expression_values.mean() == pytest.approx(overexpression_scale, rel=1e-1)

    # Check if the mean expression value between target_gene and gene_to_compare
    # in cell_to_compare is approximately equal to overexpression_scale
    cell_to_compare_expression_values = adata.X[adata.obs["cell_type"] == cell_to_compare, adata.var_names == target_gene].flatten()

    assert  cell_to_compare_expression_values.mean() / target_gene_expression_values.mean() == pytest.approx(overexpression_scale, rel=1e-1)\
    
    # Check if the mean expression value between target_gene and control_gene
    # in target_cell is approximately equal to 1
    control_gene_expression_values = adata.X[adata.obs["cell_type"] == target_cell, adata.var_names == control_gene].flatten()

    assert  control_gene_expression_values.mean() / target_gene_expression_values.mean() == pytest.approx(1, rel=1e-1)

@pytest.mark.parametrize(
        "num_batch",
        [   
            1,
            2,
            3,
        ]
)
def test_batch_effect(num_batch):
    output_path = "benccchmarker/tmp/test_output"
    data_generator = DataGenerator()

    data_generator.simulate_denovo(
        num_cells=1000,
        num_genes=100,
        num_lr_pairs=10,
        num_batch=num_batch,
        output_path=output_path
    )

    simulated_data_path = f"{output_path}/simulated_data.h5ad"

    adata = anndata.read_h5ad(simulated_data_path)

    # Check if the batch number is correctly generated
    assert num_batch == adata.obs["batch"].nunique()

@pytest.mark.parametrize(
    "num_cells, num_genes, test_data_path",
    [
        (5000, 5000, "benccchmarker/test/data/test_data_simulation_from_reference.h5ad"),
    ]       
)
def test_simulation_from_reference_output(num_cells, num_genes, test_data_path):
    output_path = "benccchmarker/tmp/test_output"
    adata = anndata.read_h5ad(test_data_path)

    data_generator = DataGenerator()

    data_generator.simulate_from_reference(
        adata,
        num_lr_pairs=10,
        overexpression_scale=1.25,
        output_path=output_path,
        output_filename="simulated_data",
    )

    adata = anndata.read_h5ad(f"{output_path}/simulated_data.h5ad")

    assert adata.shape == (num_cells, num_genes)

# @pytest.mark.parametrize(
#         "differential_interaction_map, expected_output",
#         [   
#             (False, 1),
#             (True, 2),
#         ]
# )
# def test_differential_interaction_map(differential_interaction_map, expected_output):
#     output_path = "benccchmarker/tmp/test_output"
#     data_generator = DataGenerator()

#     data_generator.simulate_denovo(
#         num_cells=1000,
#         num_genes=100,
#         num_lr_pairs=10,
#         differential_interaction_map=differential_interaction_map,
#         output_path=output_path
#     )

#     simulated_data_path = f"{output_path}/simulated_data.h5ad"

#     adata = anndata.read_h5ad(simulated_data_path)

#     # Check if the batch number is correctly generated
#     assert expected_output == adata.obs["condition"].nunique()


# IGNORE THIS TEST FOR NOW
# @pytest.mark.parametrize(
#         # The logic of this test is the same with the `test_expression_values`
#         "num_cells, num_genes, num_lr_pairs, overexpression_scale, differential_interaction, differential_interaction_map, non_differential_interaction_map, expected_condition_size",
#         [
#             (1000, 100, 10, 1.5, True, {"source": "CellType0", "target": "CellType0", "ligand": "CCL21", "receptor": "ACKR4"}, {"source": "CellType1", "target": "CellType1", "ligand": "TNFSF11", "receptor": "TNFRSF11A"}, 2)
#         ]
# )
# def test_differential_interaction_map_output_consistency(num_cells, num_genes, num_lr_pairs, overexpression_scale, differential_interaction, differential_interaction_map, non_differential_interaction_map, expected_condition_size):
#     output_path = "benccchmarker/tmp/test_output"
#     data_generator = DataGenerator()

#     data_generator.simulate_denovo(
#         num_cells=num_cells,
#         num_genes=num_genes,
#         num_lr_pairs=num_lr_pairs,
#         overexpression_scale=overexpression_scale,
#         differential_interaction=differential_interaction,
#         output_path=output_path
#     )

#     simulated_data_path = f"{output_path}/simulated_data.h5ad"

#     adata = anndata.read_h5ad(simulated_data_path)

#     # Check if the number of conditions is correct
#     assert adata.obs["condition"].nunique() == expected_condition_size

#     # Check if the ratio of the mean expression values between the different conditions is approximately equal to overexpression_scale for the differential interactions

#     condition_zero_adata = adata[adata.obs["condition"] == 0]
#     target_gene_condition_zero_adata = condition_zero_adata.X[condition_zero_adata.obs["cell_type"] == differential_interaction_map["source"], condition_zero_adata.var_names == differential_interaction_map["ligand"]].flatten()
    
#     condition_one_adata = adata[adata.obs["condition"] == 1]
#     target_gene_condition_one_adata = condition_one_adata.X[condition_one_adata.obs["cell_type"] == differential_interaction_map["source"], condition_one_adata.var_names == differential_interaction_map["ligand"]].flatten()

#     assert target_gene_condition_one_adata.mean() / target_gene_condition_zero_adata.mean() == pytest.approx(overexpression_scale, rel=1e-1)

#     # Check if the ratio of the mean expression values between the different conditions is approximately equal to 1 for the non-differential interactions

#     target_gene_condition_zero_adata = condition_zero_adata.X[condition_zero_adata.obs["cell_type"] == non_differential_interaction_map["source"], condition_zero_adata.var_names == non_differential_interaction_map["ligand"]].flatten()

#     target_gene_condition_one_adata = condition_one_adata.X[condition_one_adata.obs["cell_type"] == non_differential_interaction_map["source"], condition_one_adata.var_names == non_differential_interaction_map["ligand"]].flatten()

#     assert target_gene_condition_one_adata.mean() / target_gene_condition_zero_adata.mean() == pytest.approx(1, rel=1e-1)




