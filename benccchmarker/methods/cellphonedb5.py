import os
from pathlib import Path
import sys
import time
import tracemalloc
import yaml

import anndata
from cellphonedb.src.core.methods import cpdb_analysis_method
import pandas as pd
import scanpy as sc

def main():

    # Get arguments from the command line
    args = sys.argv

    tracemalloc.start()
    start_time = time.time()

    current_file = Path(__file__)

    if args[2] == "hsapiens":
        adata = anndata.read_h5ad(args[1])

        config = yaml.safe_load(open(args[3]))

        # Normalise data
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)

        cpdb_file_path = current_file.parent.parent / "datasets" / "defaults" / "cellphonedb5" / "database" / "cellphonedb.zip"

        main_output_path = "/".join(args[4].split("/")[:-1])
        output_path = f"{main_output_path}/tmp/cellphonedb5"

        os.makedirs(output_path, exist_ok=True)
        
        counts_file_path = f"{output_path}/normalised_adata.h5ad"
        meta_file_path = f"{output_path}/metadata.csv"
        microenvs_file_path = f"{output_path}/microenvironment.csv"

        adata.write_h5ad(counts_file_path)

        # Generate metadata file
        metadata = adata.obs[['cell_type']].reset_index()
        metadata.columns = ['barcode_sample', 'cell_type']

        metadata.to_csv(meta_file_path, index=False)

        # Generate microenvironment file
        # Create a new dataframe from metadata containing only unique cell_type columns
        microenvironment_df = metadata.drop_duplicates(subset=['cell_type'])[["cell_type"]]
        microenvironment_df["microenvironment"] = "Env1"

        microenvironment_df.to_csv(microenvs_file_path, index=False)

        end_preprocessing_time = time.time()

        cpdb_results = cpdb_analysis_method.call(
            cpdb_file_path = cpdb_file_path,           
            meta_file_path = meta_file_path,           
            counts_file_path = counts_file_path,       
            counts_data = 'hgnc_symbol',               
            microenvs_file_path = microenvs_file_path, 
            score_interactions = True,                 
            output_path = output_path,                    
            separator = '|',                           
            threads = config["params"]["cpdb_analysis_method.call"]["threads"],                               
            threshold = config["params"]["cpdb_analysis_method.call"]["threshold"],                           
            result_precision = config["params"]["cpdb_analysis_method.call"]["result_precision"],                      
            debug = False,                             
            output_suffix = None                       
        )

        tracemalloc.stop()

        result_df = cpdb_results["interaction_scores"]

        end_inference_time = time.time()

        result_df = result_df.melt(id_vars=['id_cp_interaction', 'interacting_pair', 'partner_a', 'partner_b', 'gene_a', 'gene_b', 'secreted', 'receptor_a', 'receptor_b', 'annotation_strategy', 'is_integrin', 'directionality', 'classification'], var_name='interaction', value_name='interaction_score')
        result_df[['source', 'target']] = result_df['interaction'].str.split('|', expand=True)
        result_df = result_df[["gene_a", "gene_b", "source", "target", "interaction_score"]].rename(columns={"gene_a": "ligand", "gene_b": "receptor"})
        
        result_df = result_df[result_df["interaction_score"] > 0]
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
        raise ValueError("Currently only hsapiens is supported")
    
if __name__ == "__main__":
    main()