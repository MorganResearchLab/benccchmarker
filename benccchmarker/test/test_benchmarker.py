import anndata
import pytest

from benccchmarker.core.benchmarker import Benchmarker

@pytest.mark.parametrize(
        "ground_truth_file, prediction_file, ground_truth_data_params_file, tp, auc, tpr",
        [
            ("benchmarker_simulated_data_ground_truth.csv", "benchmarker_simulated_data_prediction.csv", "benchmarker_simulated_data_ground_truth_data_params.json", 129, 1, 1),
            ("benchmarker_simulated_data_ground_truth.csv", "benchmarker_simulated_data_prediction_half.csv", "benchmarker_simulated_data_ground_truth_data_params.json", 64, 0.75, 0.5),
        ]
)
def test_benchmarker_results(ground_truth_file, prediction_file, ground_truth_data_params_file, tp, auc, tpr):
    output_path = "benccchmarker/test/data"
    
    benchmarker = Benchmarker(
        ground_truth_file_path = f"{output_path}/{ground_truth_file}",
        prediction_file_path = f"{output_path}/{prediction_file}",
        ground_truth_data_params_path = f"{output_path}/{ground_truth_data_params_file}",
        output_dir = output_path
    )

    benchmarker.run()

    benchmarker_result = benchmarker.get_results()
    
    assert benchmarker_result["tp"] == pytest.approx(tp, 0.1)
    assert benchmarker_result["auc"] == pytest.approx(auc, 0.1)
    assert benchmarker_result["tpr"] == pytest.approx(tpr, 0.1)

@pytest.mark.parametrize(
        "ground_truth_file, prediction_file, ground_truth_data_params_file",
        [
            ("benchmarker_simulated_data_ground_truth.csv", "benchmarker_simulated_data_prediction_missing_columns.csv", "benchmarker_simulated_data_ground_truth_data_params.json"),
        ]
)
def test_benchmarker_missing_column_error(ground_truth_file, prediction_file, ground_truth_data_params_file):
    output_path = "benccchmarker/test/data"
    
    benchmarker = Benchmarker(
        ground_truth_file_path = f"{output_path}/{ground_truth_file}",
        prediction_file_path = f"{output_path}/{prediction_file}",
        ground_truth_data_params_path = f"{output_path}/{ground_truth_data_params_file}",
        output_dir = output_path
    )

    with pytest.raises(ValueError):
        benchmarker.run()