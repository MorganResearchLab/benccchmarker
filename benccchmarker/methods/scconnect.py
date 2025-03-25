import os
from pathlib import Path
import sys
import time
import tracemalloc
import yaml

import anndata
import pandas as pd
import scanpy as sc

import scConnect as cn

def main():
    args = sys.argv

    tracemalloc.start()
    start_time = time.time()

    adata = anndata.read_h5ad(args[1])
    config = yaml.safe_load(open(args[3]))

    species = args[2]

    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)

    cn.database.version = config["params"]["database"]["version"]

    adata = cn.genecall.meanExpression(
        adata, 
        groupby="cell_type", 
        normalization=config["params"]["meanExpression"]["normalization"], 
        use_raw=config["params"]["meanExpression"]["use_raw"], 
        transformation=config["params"]["meanExpression"]["transformation"]
    )
    adata = cn.connect.ligands(adata)
    adata = cn.connect.receptors(adata)

    end_preprocessing_time = time.time()

    adata = cn.connect.specificity(
        adata, 
        n=config["params"]["specificity"]["n"], 
        groupby="cell_type",
        organism=species,
        transformation=config["params"]["specificity"]["transformation"],
        emperical=config["params"]["specificity"]["emperical"],
        merge_dist=config["params"]["specificity"]["merge_dist"]
    )

    edges = cn.connect.interactions(
        emitter=adata, 
        target=adata, 
        self_reference=config["params"]["interactions"]["self_reference"],
        corr_pval=config["params"]["interactions"]["corr_pval"]
    )

    end_inference_time = time.time()

    try:
        result_df = pd.DataFrame(edges)
        result_df.rename(columns={0: "source", 1: "target"}, inplace=True)

        result_df = pd.concat([result_df.drop(columns=2), result_df[2].apply(pd.Series)], axis=1)
        result_df["receptor"] = result_df["receptorgene"].str.extract(r"\[(.*)\]")
        result_df["receptor"] = result_df["receptor"].str.replace("'", "").str.split(", ")

        result_df = result_df.explode("receptorgene")
        result_df["label"] = result_df["source"] + "---" + result_df["ligand"] + "---" + result_df["target"] + "---" + result_df["receptor"]
    except:
        result_df = pd.DataFrame(columns=["source", "target", "ligand", "receptor", "label"])

    result_df.to_csv(args[4], index=False)
        
    _, peak = tracemalloc.get_traced_memory()

    met_path = args[4].replace(".csv", "_met.csv")

    preprocessing_time = end_preprocessing_time - start_time
    inference_time = end_inference_time - end_preprocessing_time

    met_df = pd.DataFrame(data={"peak_ram_used_mib": [peak], "preprocessing_time": [preprocessing_time], "inference_time": [inference_time]})
    met_df.to_csv(met_path, index=False)

if __name__ == "__main__":
    main()