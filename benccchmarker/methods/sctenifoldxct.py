import sys
import time
import tracemalloc
import yaml

import anndata
import pandas as pd
import scanpy as sc

import scTenifoldXct as st

def main():

    # Get arguments from the command line
    args = sys.argv

    tracemalloc.start()
    start_time = time.time()

    if args[2] == "hsapiens":
        adata = anndata.read_h5ad(args[1])

        config = yaml.safe_load(open(args[3]))

        # Normalise data
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)

        # Get all unique cell types
        cell_types = sorted(adata.obs["cell_type"].unique())

        end_preprocessing_time = time.time()

        result_df = pd.DataFrame()

        for source_cell_type in cell_types:
            for target_cell_type in cell_types:
                print(source_cell_type, target_cell_type)
                try:
                    xct = st.scTenifoldXct(
                        data = adata, 
                        source_celltype = source_cell_type,
                        target_celltype = target_cell_type,
                        obs_label = 'cell_type',
                        rebuild_GRN = config["params"]["scTenifoldXct"]["rebuild_GRN"],
                        alpha = config["params"]["scTenifoldXct"]["alpha"],
                        mu = config["params"]["scTenifoldXct"]["mu"],
                        scale_w = config["params"]["scTenifoldXct"]["scale_w"],
                        n_dim = config["params"]["scTenifoldXct"]["n_dim"],
                        verbose = config["params"]["scTenifoldXct"]["verbose"]
                    )

                    xct.get_embeds(train = True)

                    xct_pairs = xct.null_test()
                    xct_pairs["source"] = source_cell_type
                    xct_pairs["target"] = target_cell_type

                    result_df = pd.concat([result_df, xct_pairs])
                except:
                    print(f"Error in generating interactions for {source_cell_type} and {target_cell_type}")

        tracemalloc.stop()
        end_inference_time = time.time()
        
        result_df.reset_index(drop=True, inplace=True)
        result_df = result_df[["ligand", "receptor", "source", "target"]]
        result_df["label"] = result_df["source"] + "---" + result_df["ligand"] + "---" + result_df["target"] + "---" + result_df["receptor"]

        result_df.to_csv(args[4], index=False)

        _, peak = tracemalloc.get_traced_memory()

        met_path = args[4].replace(".csv", "_met.csv")

        # Create a dataframe to store the memory usage in a column called peak_ram_used_mib
        preprocessing_time = end_preprocessing_time - start_time
        inference_time = end_inference_time - end_preprocessing_time

        met_df = pd.DataFrame(data={"peak_ram_used_mib": [peak], "preprocessing_time": [preprocessing_time], "inference_time": [inference_time]})
        met_df.to_csv(met_path, index=False)

    else:
        raise ValueError("Species not supported and currently only support 'hsapiens'")

if __name__ == "__main__":
    main()
