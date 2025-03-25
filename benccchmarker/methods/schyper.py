import pandas as pd
import scanpy as sc
import numpy as np
import torch
import scipy
import matplotlib
import matplotlib.pyplot as plt
import os
from itertools import product
from scHyper import dataprocess as dp

def main():
    args = sys.argv

    tracemalloc.start()
    start_time = time.time()

    adata = anndata.read_h5ad(args[1])

    config = yaml.safe_load(open(args[3]))

    # Normalise data
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    end_preprocessing_time = time.time()

    # adata = dp.meanExpression(adata, type="mean", groupby="cell_type")
    # adata, ligand_receptor_data = dp.process_ligands_receptors(adata, "human", highly_variable=False)
    adata = dp.meanExpression(adata, type="mean", groupby="cell_type")
    adata, ligand_receptor_data = dp.process_ligands_receptors(adata, "human", highly_variable=False)

    interaction_tensor = dp.generate_tensor(adata, ligand_receptor_data)

    triplets, weights, validlrindices = dp.generate_triplets_weights_validlrindices(interaction_tensor)
    validlrs, invalidlrs = dp.generate_validlrs_invalidlrs(validlrindices, ligand_receptor_data)

    validsenderindices, validreceiverindices = dp.generate_validsenderindices_validreceiverindices(interaction_tensor)
    validsenders, invalidsenders, validreceivers, invalidreceivers = dp.generate_validsenders_validreceivers(adata, validsenderindices, validreceiverindices)

    triplets = dp.update_triplets(triplets)
    weights = dp.update_weights(weights)

    train_triplets, test_triplets, train_weights, test_weights = dp.generate_train_test(triplets, weights)
    train_nums_type, test_nums_type = dp.generate_nums_type(train_triplets, test_triplets)

    save_path='out'

    if not os.path.exists(save_path):
        os.makedirs(save_path)

    np.savez(os.path.join(save_path, 'train_data.npz'), nums_type=train_nums_type, train_data=train_triplets, train_weight=train_weights)
    np.savez(os.path.join(save_path, 'test_data.npz'), nums_type=train_nums_type, test_data=test_triplets, test_weight=test_weights)

    use_to_predict = dp.use_to_predict(triplets)
    np.savez(os.path.join(save_path, 'use_to_predict.npz'), use_to_predict=use_to_predict)

    # os.system(f"!cd scHyper/scHyper && python model/main_torch.py --data ../../out")

    df_nn, candidates = dp.genenrate_df_nn_candidates(validlrs, validsenders, validreceivers, triplets, use_to_predict)
    df_enriched, tensor_pval = dp.null_test(df_nn, candidates, pval=0.05, plot=False)

    df_enriched.reset_index(inplace=True)
    df_enriched[['source', 'ligand', 'receptor', 'target']] = df_enriched['index'].str.split('_', expand=True)
    result_df = df_enriched[["source", "ligand", "receptor", "target", "p_val", "prob"]]

    end_inference_time = time.time()
    
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