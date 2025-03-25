import os
import sys
import subprocess
import time
import tracemalloc
import yaml
import psutil

import anndata
import pandas as pd

tracemalloc.start()

def main():

    # Get arguments from the command line
    start_time = time.time()

    args = sys.argv

    environment = args[5] if len(args) > 5 else None

    adata = anndata.read_h5ad(args[1])
    config = yaml.safe_load(open(args[3]))

    main_output_path = "/".join(args[4].split("/")[:-1])
    output_path = f"{main_output_path}/tmp/disir"

    os.makedirs(output_path, exist_ok=True)

    adata_df = adata.to_df().T

    expression_values_path = f"{output_path}/expression_values.csv"
    gene_names_path = f"{output_path}/gene_names.csv"
    cell_types_path = f"{output_path}/cell_types.csv"

    # Save expression values
    adata_df.to_csv(expression_values_path, index=False, header=False)

    # Save gene names
    adata_df.index.to_series().to_csv(gene_names_path, index=False, header=False)

    # Save cell types
    adata.obs["cell_type"].to_csv(f"{output_path}/cell_types.csv", index=False, header=False)

    end_preprocessing_time = time.time()

    if config["params"]["disir_package"]["subunit_interactions_path"] == None:
        raise ValueError("You're trying to run DISIR, however no path to the subunit interactions file are provided, please provide a path to the subunit interactions file.")
    
    # Start monitoring memory usage
    process = psutil.Process(os.getpid())
    memory_info_before = process.memory_info().rss


    subprocess.run(
        [
            f"{environment}/bin/disir_package" if environment else "disir_package",
            "-tf",
            str(config["params"]["disir_package"]["tf"]),
            "-te",
            str(config["params"]["disir_package"]["te"]),
            "-tp",
            str(config["params"]["disir_package"]["tp"]),
            "-iv",
            str(config["params"]["disir_package"]["iv"]),
            "--scRNA-path",
            expression_values_path,
            "--cell-type-path",
            cell_types_path,
            "--gene-path",
            gene_names_path,
            "--subunit-interactions-path",
            config["params"]["disir_package"]["subunit_interactions_path"],
            "--output-directory-path",
            output_path,
        ], check=True, capture_output=True
    )

    memory_info_after = process.memory_info().rss
    memory_used = (memory_info_after - memory_info_before) / (1024 * 1024)  # Convert to MiB

    tracemalloc.stop()

    end_inference_time = time.time()

    interaction_output_dirs = [i for i in os.listdir(output_path) if os.path.isdir(f"{output_path}/{i}")]

    result_df = pd.DataFrame()

    for interaction_output_dir in interaction_output_dirs:
        ligand, receptor = interaction_output_dir.split("_")
        
        temp_df = pd.read_csv(f"{output_path}/{interaction_output_dir}/heatmaps/Heatmap.csv")

        temp_df = temp_df.melt(id_vars=["Unnamed: 0"], var_name="cell_type", value_name="interaction_strength")
        temp_df.rename(columns={"Unnamed: 0": "source", "cell_type": "target"}, inplace=True)

        # Filter out rows with 0 interaction_strength
        temp_df = temp_df[temp_df["interaction_strength"] != 0]

        temp_df["ligand"] = ligand
        temp_df["receptor"] = receptor

        result_df = pd.concat([result_df, temp_df])

    result_df["label"] = result_df["source"] + "---" + result_df["ligand"] + "---" + result_df["target"] + "---" + result_df["receptor"]

    result_df.to_csv(args[4], index=False)

    _, peak = tracemalloc.get_traced_memory()

    met_path = args[4].replace(".csv", "_met.csv")

    # Create a dataframe to store the memory usage in a column called peak_ram_used_mib
    preprocessing_time = end_preprocessing_time - start_time
    inference_time = end_inference_time - end_preprocessing_time

    met_df = pd.DataFrame(data={"peak_ram_used_mib": [memory_used], "preprocessing_time": [preprocessing_time], "inference_time": [inference_time]})
    met_df.to_csv(met_path, index=False)

if __name__ == "__main__":
    main()