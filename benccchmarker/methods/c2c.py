from pathlib import Path
import sys
import time
import tracemalloc
import yaml

import anndata
import pandas as pd
import scanpy as sc
import cell2cell

def main():

    # Get arguments from the command line
    args = sys.argv

    tracemalloc.start()
    start_time = time.time()

    adata = anndata.read_h5ad(args[1])

    config = yaml.safe_load(open(args[3]))

    lr_pairs = pd.read_csv(f'{config["ligand_receptor_pair_db"]}')
    meta = adata.obs.copy()

    end_preprocessing_time = time.time()

    interactions = cell2cell.analysis.SingleCellInteractions(
        rnaseq_data=adata.to_df().T,
        ppi_data=lr_pairs,
        metadata=meta,
        interaction_columns=('ligand_symbol', 'receptor_symbol'),
        communication_score=config["params"]["SingleCellInteractions"]["communication_score"],
        expression_threshold=config["params"]["SingleCellInteractions"]["expression_threshold"],
        cci_score=config["params"]["SingleCellInteractions"]["cci_score"],
        cci_type=config["params"]["SingleCellInteractions"]["cci_type"],
        aggregation_method=config["params"]["SingleCellInteractions"]["aggregation_method"],
        barcode_col=config["params"]["SingleCellInteractions"]["barcode_col"],
        celltype_col=config["params"]["SingleCellInteractions"]["celltype_col"],
        complex_sep=config["params"]["SingleCellInteractions"]["complex_sep"],
        verbose=config["params"]["SingleCellInteractions"]["verbose"],
    )

    interactions.compute_pairwise_communication_scores()

    result_df = interactions.permute_cell_labels(
        evaluation=config["params"]["permute_cell_labels"]["evaluation"],
        permutations=config["params"]["permute_cell_labels"]["permutations"], 
        fdr_correction=config["params"]["permute_cell_labels"]["fdr_correction"],
        verbose=config["params"]["permute_cell_labels"]["verbose"],
        random_state=config["params"]["permute_cell_labels"]["random_state"],
    )

    end_inference_time = time.time()

    result_df = result_df.reset_index().melt(id_vars='index', var_name='interaction', value_name='p-value')

    # split "index" column (CCL21, ACKR4) into "ligand" and "receptor" columns
    result_df['ligand'] = result_df['index'].str[0]
    result_df['receptor'] = result_df['index'].str[1]
    # Save column name as source and target and add to result_df
    result_df[["source", "target"]] = result_df["interaction"].str.split(";", expand=True).rename(columns={0: "source", 1: "target"})

    result_df["label"] = result_df["source"] + "---" + result_df["ligand"] + "---" + result_df["target"] + "---" + result_df["receptor"]

    result_df.to_csv(args[4], index=False)

    _, peak = tracemalloc.get_traced_memory()

    met_path = args[4].replace(".csv", "_met.csv")

    # Create a dataframe to store the memory usage in a column called peak_ram_used_mib
    preprocessing_time = end_preprocessing_time - start_time
    inference_time = end_inference_time - end_preprocessing_time

    met_df = pd.DataFrame(data={"peak_ram_used_mib": [peak], "preprocessing_time": [preprocessing_time], "inference_time": [inference_time]})
    met_df.to_csv(met_path, index=False)

if __name__ == "__main__":
    main()