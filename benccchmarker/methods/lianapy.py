import os
from pathlib import Path
import sys
import time
import tracemalloc
import yaml

import anndata
import liana as li
import pandas as pd
import scanpy as sc

def main():
    args = sys.argv

    tracemalloc.start()
    start_time = time.time()

    adata = anndata.read_h5ad(args[1])

    config = yaml.safe_load(open(args[3]))

    # Normalise data
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    adata.raw = adata.copy()

    end_preprocessing_time = time.time()

    li.mt.rank_aggregate(
        adata,
        groupby='cell_type',
        resource_name=config["params"]["rank_aggregate"]["resource_name"],
        expr_prop=config["params"]["rank_aggregate"]["expr_prop"],
        min_cells=config["params"]["rank_aggregate"]["min_cells"],
        aggregate_method=config["params"]["rank_aggregate"]["aggregate_method"],
        n_perm=config["params"]["rank_aggregate"]["n_perm"],
        n_jobs=config["params"]["rank_aggregate"]["n_jobs"],
        seed=config["params"]["rank_aggregate"]["seed"],
        verbose=True
    )
    
    end_inference_time = time.time()
    
    result_df = adata.uns['liana_res']
    result_df.rename(columns={"ligand_complex": "ligand", "receptor_complex": "receptor"}, inplace=True)
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