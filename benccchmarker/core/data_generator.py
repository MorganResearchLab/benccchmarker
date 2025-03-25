import json
import logging
import os
from pathlib import Path
import warnings

# warnings.filterwarnings("ignore")

import anndata
import numpy as np
import pandas as pd
import scanpy as sc
import scipy
from scipy.stats import nbinom, lognorm, bernoulli

from benccchmarker.utils import convert_anndata_to_seurat, get_n_p_from_mean_dispersion

class DataGenerator():
    """
    DataGenerator class for generating both denovo simulated data and simulated data from reference of single-cell RNA-seq data with
    injected cell-to-cell communication.

    Parameters
    ----------
    top_lr_pair_path : str
        Path to the file containing ligand-receptor pairs.
    seed : int  
        Random seed.
    """
    def __init__(
        self,
        top_lr_pair_path="default",
        seed=42,
    ):
        self.current_file = Path(__file__)

        if top_lr_pair_path == "default":
            top_lr_pair_path = self.current_file.parent.parent / "datasets" / "lrdb_top_100_pairs.csv"

        self.top_lr_pair_path = top_lr_pair_path
        self.lr_pair_df = pd.read_csv(top_lr_pair_path)
        self.lr_pair_df[["ligand", "receptor"]] = self.lr_pair_df["key"].str.upper().str.split("---", n=1, expand=True)

        np.random.seed(seed)
        self.seed = seed

        # Initialise empty variables to store data
        self.simulated_data = None
        self.synthetic_data = None

        # Create a logger
        logging.basicConfig(
            level=logging.INFO,
            format="%(asctime)s - %(message)s",
            datefmt="%Y-%m-%d %H:%M:%S",
        )

    def _check_distribution_parameters(self, baseline_mean, dispersion):
        """
        Check if the mean expression and dispersion parameters are valid.
        Raise exception if baseline_mean / (dispersion ** 2) and (baseline_mean ** 2) / ((dispersion ** 2) - baseline_mean) < 0

        Parameters
        ----------
        baseline_mean : float
            Mean expression level.
        dispersion : float
            Dispersion parameter.
        
        Returns
        -------
        bool
            True if the mean expression and dispersion parameters are valid.
        """
        if baseline_mean / (dispersion ** 2) < 0 or (baseline_mean ** 2) / ((dispersion ** 2) - baseline_mean) < 0:
            logging.error("Invalid mean expression and dispersion parameters, baseline_mean / (dispersion ** 2) and (baseline_mean ** 2) / ((dispersion ** 2) - baseline_mean) has to be greater than 0.")
            raise ValueError("Invalid mean expression and dispersion parameters, baseline_mean / (dispersion ** 2) and (baseline_mean ** 2) / ((dispersion ** 2) - baseline_mean) has to be greater than 0.")
        return True

    def _generate_lr_pairs(self, lr_pair_df, num_lr_pairs, output_path, gene_names=None):
        """
        Generate ligand-receptor pairs.

        Parameters
        ----------
        lr_pair_df : pd.DataFrame
            Dataframe containing ligand-receptor pairs.
        num_lr_pairs : int
            Number of ligand-receptor pairs to generate.
        seed : int
            Random seed.

        Returns
        -------
        lr_pairs: dict
            Dictionary containing ligand-receptor pairs.
        """

        # Sample ligand-receptor pairs
        lr_pair_sample_df = lr_pair_df.sample(
            num_lr_pairs, 
            random_state=self.seed
        ).copy()

        lr_pair_sample_df.to_csv(f"{output_path}/ground_truth_lr_pairs.csv", index=False)

        # Generate ligand-receptor pair dictionary
        lr_pairs = {}

        for ligand, receptor in [tuple(x) for x in lr_pair_sample_df[["ligand", "receptor"]].to_numpy()]:
            if gene_names != None:
                if ligand in gene_names and receptor in gene_names:
                    lr_pairs[ligand] = receptor
            else:
                lr_pairs[ligand] = receptor
        return lr_pairs

    def _generate_cell_type_labels(
        self,
        num_cells,
        num_cell_types,
        cell_type_proportions=[]
    ):
        """
        Generate cell type labels.

        Parameters
        ----------
        num_cells : int
            Number of cells.
        num_cell_types : int
            Number of cell types.

        Returns
        -------
        list
            List of cell type labels.
        """

        if cell_type_proportions != []:
            cell_type_labels = np.random.choice(range(num_cell_types), num_cells, p=cell_type_proportions)
            cell_type_labels = [f"CellType{i}" for i in cell_type_labels]

            return cell_type_labels
        else:
            cell_type_labels = np.random.choice(range(num_cell_types), num_cells)
            cell_type_labels = [f"CellType{i}" for i in cell_type_labels]

            return np.array(cell_type_labels)

    def _generate_source_target_pairs(
        self,
        num_cell_types
    ):
        """
        Generate source target pairs.

        Parameters
        ----------
        num_cell_types : int
            Number of cell types.

        Returns
        -------
        np.array
            Array of source target pairs.
        """
        source_target_pairs = np.stack(np.meshgrid(range(num_cell_types), range(num_cell_types)), -1).reshape(-1, 2)
        source_target_pairs = np.array([f"cell_{source_target_pair[0]}-cell_{source_target_pair[1]}" for source_target_pair in source_target_pairs])

        return source_target_pairs

    def _generate_cell_type_receptor_ligand_dict(
        self,
        cell_types,
        lr_pairs
    ):
        """
        Generate cell type receptor ligand dictionary.

        Parameters
        ----------
        cell_types : list
            List of cell types.
        lr_pair_sample_df : pd.DataFrame
            Dataframe containing ligand-receptor pairs.

        Returns
        -------
        dict
            Dictionary containing cell type receptor ligand pairs.
        """
        cell_type_receptor_ligand_dict = {}
        for i, cell_type in enumerate(cell_types):
            ligand_genes = np.array(sorted(list(set(lr_pairs.keys()))))

            cell_type_receptor_ligand_dict[cell_type] = {}
            cell_type_receptor_ligand_dict[cell_type]["ligands"] = ligand_genes

            # receptor_genes = np.array(lr_pair_sample_df["receptor"].to_list())
            # cell_type_receptor_ligand_dict[cell_type]["receptors"] = receptor_genes

        return cell_type_receptor_ligand_dict

    def _generate_interaction_maps(
        self,
        cell_type_receptor_ligand_dict,
        lr_pairs
    ):
        """
        Generate interaction maps.

        Parameters
        ----------
        cell_type_receptor_ligand_dict : dict
            Dictionary containing cell type receptor ligand pairs.
        lr_pairs : dict
            Dictionary containing ligand-receptor pairs.

        Returns
        -------
        interaction_maps : np.array
            Array of interaction maps it looks like this [{"source": "CellType0", "target": "CellType0", "ligand": "CCL21", "receptor": "ACKR4"}].
        """
        interaction_maps = []
        for source_cell_type in cell_type_receptor_ligand_dict.keys():
            for target_cell_type in cell_type_receptor_ligand_dict.keys():
                for ligand in cell_type_receptor_ligand_dict[source_cell_type]["ligands"]:
                    interaction_map = {}

                    interaction_map["source"] = source_cell_type
                    interaction_map["target"] = target_cell_type
                    interaction_map["ligand"] = ligand
                    interaction_map["receptor"] = lr_pairs[ligand]

                    interaction_maps.append(interaction_map)

        interaction_maps = np.array(interaction_maps)

        interaction_map_probas = np.random.uniform(0, 1, len(interaction_maps))
        interaction_maps = interaction_maps[np.argwhere(interaction_map_probas > 0.9).T[0]]

        return interaction_maps

    def _find_lr_genes_intersection(
        self,
        adata,
        lrdb_path="default"
    ):
        """
        Find the intersection of ligand-receptor genes in the anndata object and the ligand-receptor database.

        Parameters
        ----------
        adata : anndata.AnnData
            Anndata object.
        lrdb_path : str
            Path to the ligand-receptor database.
        
        Returns
        -------
        lr_genes_in_matrix : list
            List of ligand-receptor genes in the anndata object.
        """
        if lrdb_path == "default":
            lrdb_path = self.current_file.parent.parent / "datasets" / "lrdb_key.parquet.gzip"

        lrdb = pd.read_parquet(lrdb_path)
        lr_genes = sorted(list(set(lrdb["ligand"]).union(set(lrdb["receptor"]))))
        lr_genes = [lr_gene.upper() for lr_gene in lr_genes]
        lr_genes_in_matrix = sorted(list(set(lr_genes).intersection(set(adata.var_names.to_list()))))
        
        return lr_genes_in_matrix

    def _generate_ground_truth(self, interaction_maps):
        """
        Generate ground truth dataframe from interaction maps.

        Parameters
        ----------
        interaction_maps : np.array
            Array of interaction maps.

        Returns
        -------
        ground_truth_df : pd.DataFrame
            Dataframe containing ground truth of cell-to-cell communication.
        """
        ground_truth_df = pd.DataFrame(interaction_maps.tolist())
        ground_truth_df["label"] = ground_truth_df.apply(lambda row: '---'.join(row), axis=1)

        return ground_truth_df

    def _validate_inputs(self, num_lr_pairs, num_genes, num_cell_types, marker_genes_per_cell_type, output_path):
        """
        Validate the inputs.

        Parameters
        ----------
        num_lr_pairs : int
            Number of ligand-receptor pairs.
        num_genes : int
            Number of genes.
        output_path : str
            Path to save the output file.
        """
        if num_lr_pairs > num_genes:
            logging.error("num_lr_pairs has to be lower than num_genes.")
            raise ValueError("num_lr_pairs has to be lower than num_genes.")

        if (num_genes - num_lr_pairs) < marker_genes_per_cell_type:
            logging.error("You don't have enough genes to assign marker genes for each cell type.")
            raise ValueError("You don't have enough genes to assign marker genes for each cell type.")
        
        if (num_genes - num_lr_pairs) < num_cell_types * marker_genes_per_cell_type:
            logging.warning("There might be multiple overlap between marker genes for each cell type.")
            warnings.warn("There might be multiple overlap between marker genes for each cell type. It is recommended to increase the number of genes.")


        if not os.path.exists(output_path):
            os.makedirs(output_path, exist_ok=True)

    def _update_interactions(self, simulated_df, interaction_maps, scenario):
        if scenario == "lherhe":
            for interaction in interaction_maps:
                simulated_df.loc[simulated_df.index == interaction["source"], interaction["ligand"]] = -999
                simulated_df.loc[simulated_df.index == interaction["target"], interaction["receptor"]] = -999

        elif scenario == "lherme":
            for interaction in interaction_maps:
                simulated_df.loc[simulated_df.index == interaction["source"], interaction["ligand"]] = -999
        else:
            for interaction in interaction_maps:
                simulated_df.loc[simulated_df.index == interaction["target"], interaction["receptor"]] = -999
        return simulated_df

     # Simulate dropout on simulated_df
    def _simulate_dropout(self, df, dropout_rate):
        mask = np.random.rand(*df.shape) < dropout_rate
        df_dropout = df.mask(mask, 0)
        return df_dropout

    def _validate_cell_type_input(self, num_cell_types, cell_type_proportions):
        if cell_type_proportions == []:
            return True

        if num_cell_types != len(cell_type_proportions):
            logging.error("Number of cell types and cell type proportions do not match.")
            raise ValueError("Number of cell types and cell type proportions do not match.")

        if sum(cell_type_proportions) != 1:
            logging.error("Cell type proportions do not sum to 1.")
            raise ValueError("Cell type proportions do not sum to 1.")

        return True

    def simulate_denovo(
        self,
        num_cells=5000,
        num_genes=2500,
        num_lr_pairs=10,
        num_cell_types=5,
        cell_type_proportions=[],
        marker_genes_per_cell_type=10,
        method="overexpression",
        communication_strength=2,
        baseline_mean=10,
        marker_mean=25,
        dispersion=5,
        dropout_a=1.5,
        dropout_b=2,
        lib_size_mean=1,
        lib_size_sigma=0.5,
        num_batch=1,
        batch_factor=0.5,
        batch_genes_per_batch = 30,
        batch_mean=0,
        batch_sd=0.1,
        be_method="shift_by_mean",
        differential_interaction=False,
        de_genes_per_condition = 50,
        condition_fold_change = 2,
        output_path="out",
        output_filename="simulated_data",
        scenario="lherhe",
        add_control=True
    ):

        """
        Simulate single-cell RNA-seq data.

        Parameters
        ----------
        num_cells : int
            Number of cells to simulate.
        num_genes : int
            Number of genes to simulate.
        num_lr_pairs : int
            Number of ligand-receptor pairs to simulate.
        num_cell_types : int
            Number of cell types to simulate.
        method : str
            Method to simulate data. "overexpression" will simulate the overexpression of
            ligand-receptor genes. "zero" will zero out the expression values
            of the non ligand-receptor genes.
        communication_strength : float
            Scale factor for overexpressing ligand-receptor pairs. Required if method is "overexpression".
        baseline_mean : float
            Mean expression level.
        dispersion : float
            Dispersion parameter.
        num_batch : int
            Number of batches to simulate.
        batch_factor : float
            Batch factor to scale the expression values.
        batch_mean : float
            Mean value to shift the expression values.
        batch_sd : float
            Standard deviation value to shift the expression values.
        be_method : str
            Batch effect method. "scale" will scale the expression values. "shift_by_sd" will shift the expression values by a random number sampled from normal distribution. "shift_by_mean" will shift the expression values by a random number sampled from normal distribution.
        output_path : str
            Path to save the output file.
        output_filename : str
            Name of the output file.
        scenario : str
            Name of the scenario to simulate the gene expression. 
            "lherhe" means "Ligand is Highly Expressed, Receptor is Highly Expressed".
            "lmerhe" means "Ligand is Moderately Expressed, Receptor is Highly Expressed". 
            "lherme" means "Ligand is Highly Expressed, Receptor is Moderately Expressed".
        add_control : bool
            Add control genes to the simulated data.
        """

        # Make sure num_lr_pairs is always greater than num_genes, otherwise raise an exception
        self._validate_inputs(num_lr_pairs, num_genes, num_cell_types, marker_genes_per_cell_type, output_path)
        self._validate_cell_type_input(num_cell_types, cell_type_proportions)

        # Check if the mean expression and dispersion parameters are valid
        # self._check_distribution_parameters(baseline_mean, dispersion)

        # Sample ligand-receptor pairs based on the number of ligand-receptor pairs input
        lr_pairs = self._generate_lr_pairs(self.lr_pair_df, num_lr_pairs, output_path)
        lr_genes = np.concatenate([list(lr_pairs.keys()), list(lr_pairs.values())])
        num_lr_genes = len(np.unique(lr_genes))

        cell_types = [f"CellType{i}" for i in range(num_cell_types)]
        cell_barcodes = [f"cell_{i}" for i in range(num_cells)]
        cell_type_labels = self._generate_cell_type_labels(num_cells, num_cell_types)

        # Generate source target pairs
        source_target_pairs = self._generate_source_target_pairs(num_cell_types)

        # Generte combination of cell types and ligand-receptor pairs
        cell_type_receptor_ligand_dict = self._generate_cell_type_receptor_ligand_dict(cell_types, lr_pairs)

        # Generate interaction maps
        interaction_maps = self._generate_interaction_maps(cell_type_receptor_ligand_dict, lr_pairs)
        print(f"Generated {len(interaction_maps)} interactions")
        
        # Generate anndata
        baseline_expressions = np.full((num_cells, num_genes), baseline_mean)

        non_lr_gene_path = self.current_file.parent.parent / "datasets" / "non_lr_genes_sample.csv"
        non_lr_gene_names_df = pd.read_csv(non_lr_gene_path)

        non_lr_genes = np.array(non_lr_gene_names_df.sample(num_genes - num_lr_genes, random_state=self.seed)["gene"].tolist())

        all_genes = np.unique(np.concatenate([lr_genes, non_lr_genes]))

        gene_is_overexpressed = np.zeros(num_genes, dtype=int) - 1

        overexpressed_mean = baseline_mean * communication_strength

        if differential_interaction:
            condition_labels = np.zeros(num_cells, dtype=int)
            for cell_type in np.unique(cell_type_labels):
                cell_type_idx = np.where(cell_type_labels == cell_type)[0]
                np.random.shuffle(cell_type_idx)
                split = len(cell_type_idx) // 2
                condition_labels[cell_type_idx[:split]] = 0
                condition_labels[cell_type_idx[split:]] = 1

            de_genes_control = np.random.choice(all_genes, de_genes_per_condition, replace=False)
            de_genes_treated = np.random.choice(all_genes, de_genes_per_condition, replace=False)

            de_genes_control_indices = [np.where(all_genes == gene)[0][0] for gene in de_genes_control]
            de_genes_treated_indices = [np.where(all_genes == gene)[0][0] for gene in de_genes_treated]

            treated_cells = np.where(condition_labels == 1)[0]
            baseline_expressions = baseline_expressions.astype(float)
            baseline_expressions[treated_cells[:, np.newaxis], de_genes_treated_indices] *= condition_fold_change

            control_cells = np.where(condition_labels == 0)[0]
            baseline_expressions[control_cells[:, np.newaxis], de_genes_control_indices] /= condition_fold_change


        marker_genes = {}
        marker_indices = np.zeros(num_genes)
        for cell_type in cell_types:
            markers = np.random.choice(non_lr_genes, marker_genes_per_cell_type, replace=False)
            for marker in markers:
                marker_index = np.where(all_genes == marker)[0][0]
                cell_type_index = np.where(cell_type_labels == cell_type)[0]
                baseline_expressions[cell_type_index, marker_index] = marker_mean * np.random.uniform(0.9, 1.5)
                gene_is_overexpressed[marker_index] = 1
            marker_genes[cell_type] = markers

        if scenario == "lherhe":
            for interaction_map in interaction_maps:
                ligand_index = np.where(all_genes == interaction_map["ligand"])[0][0]
                receptor_index = np.where(all_genes == interaction_map["receptor"])[0][0]

                source_index = np.where(cell_type_labels == interaction_map["source"])[0]
                target_index = np.where(cell_type_labels == interaction_map["target"])[0]

                gene_is_overexpressed[ligand_index] = 1
                gene_is_overexpressed[receptor_index] = 1

                baseline_expressions[source_index, ligand_index] = overexpressed_mean
                baseline_expressions[target_index, receptor_index] = overexpressed_mean
        elif scenario == "lmerhe":
            for interaction_map in interaction_maps:
                receptor_index = np.where(all_genes == interaction_map["receptor"])[0][0]

                target_index = np.where(cell_type_labels == interaction_map["target"])[0]

                gene_is_overexpressed[ligand_index] = 1

                baseline_expressions[target_index, receptor_index] = overexpressed_mean
        elif scenario == "lherme":
            for interaction_map in interaction_maps:
                ligand_index = np.where(all_genes == interaction_map["ligand"])[0][0]

                source_index = np.where(cell_type_labels == interaction_map["source"])[0]

                gene_is_overexpressed[ligand_index] = 1

                baseline_expressions[source_index, ligand_index] = overexpressed_mean
        else:
            raise ValueError("Invalid scenario. Please choose from 'lherhe', 'lmerhe', 'lherme'.")

        log_lib_size = np.random.normal(loc=0, scale=lib_size_sigma, size=num_cells)
        lib_size_factors = np.exp(log_lib_size)
        lib_size_factors /= np.mean(lib_size_factors)

        cell_means = baseline_expressions * lib_size_factors[:, np.newaxis]

        if num_batch > 1:
            batch_effect_strength = overexpressed_mean * 2
            batch_labels = np.random.choice(num_batch, size=num_cells)

            batch_effect_matrix = np.zeros((num_batch, num_genes))
            for batch in range(num_batch):
                perturbed_genes = np.random.choice(
                    num_genes, 
                    size=batch_genes_per_batch, 
                    replace=False
                )
                batch_effect_matrix[batch, perturbed_genes] = batch_effect_strength

            for i in range(num_cells):
                batch = batch_labels[i]
                cell_means[i, :] += batch_effect_matrix[batch, :]

            cell_means = np.clip(cell_means, a_min=0, a_max=None)

        baseline_mean_genes = np.where(gene_is_overexpressed > 0, overexpressed_mean, baseline_mean)
        dispersion_genes = 1 / (baseline_mean_genes + 1) + 0.1
        shape_genes = 1 / dispersion_genes

        scale_cell_gene = cell_means * dispersion_genes
        lambda_cell_gene = np.random.gamma(shape=shape_genes, scale=scale_cell_gene)
        counts = np.random.poisson(lambda_cell_gene)

        log_mu = np.log(cell_means + 1e-6)
        dropout_probs = 1 / (1 + np.exp(dropout_a * (log_mu - dropout_b)))
        dropout_mask = bernoulli.rvs(dropout_probs).astype(bool)
        counts[dropout_mask] = 0

        var = pd.DataFrame(index=all_genes)
        obs = pd.DataFrame(index=cell_barcodes)

        obs["cell_type"] = cell_type_labels

        adata = anndata.AnnData(
            X=counts.astype(np.float32),  # Use float32 to save memory
            obs=obs,
            var=var
        )
        
        if num_batch > 1:
            adata.obs["batch"] = batch_labels
        else:
            adata.obs["batch"] = 0

        if differential_interaction:
            adata.obs["condition"] = condition_labels
        else:
            adata.obs["condition"] = 0
            
        adata.layers['counts'] = adata.X.copy()

        return adata

        adata.write_h5ad(f'{output_path}/{output_filename}.h5ad')

        print(f"Simulated data saved to {output_path}/{output_filename}.h5ad")

        ground_truth_df = self._generate_ground_truth(interaction_maps)
        ground_truth_df.to_csv(f"{output_path}/{output_filename}_ground_truth.csv", index=False)

        print(f"Ground truth saved to {output_path}/simulation/{output_filename}_ground_truth.csv")

        self.simulated_data = adata

        print(f"Converting AnnData object to Seurat object")

        convert_anndata_to_seurat(
            f"{output_path}/{output_filename}.h5ad",
            f"{output_path}/{output_filename}_seurat.rds"
        )
        
        if os.path.exists(f'{output_path}/{output_filename}_seurat.rds'):
            print(f"Conversion complete. The Seurat object is saved as '{output_path}/{output_filename}_seurat.rds'.")
        else:
            raise Exception("Conversion failed. Please check the output path.")

        # Save the parameters to a new file called data_params.json
        data_params = {
            "num_cells": num_cells,
            "num_genes": num_genes,
            "num_lr_pairs": num_lr_pairs,
            "num_cell_types": num_cell_types,
            "cell_type_proportions": cell_type_proportions,
            "marker_genes_per_cell_type": marker_genes_per_cell_type,
            "method": method,
            "communication_strength": communication_strength,
            "baseline_mean": baseline_mean,
            "marker_mean": marker_mean,
            "dispersion": dispersion,
            "dropout_a": dropout_a,
            "dropout_b": dropout_b,
            "lib_size_mean": lib_size_mean,
            "lib_size_sigma": lib_size_sigma,
            "num_batch": num_batch,
            "batch_factor": batch_factor,
            "batch_genes_per_batch": batch_genes_per_batch,
            "batch_mean": batch_mean,
            "batch_sd": batch_sd,
            "be_method": be_method,
            "differential_interaction": differential_interaction,
            "output_path": output_path,
            "output_filename": output_filename,
            "scenario": scenario,
            "add_control": add_control
        }

        with open(f'{output_path}/data_params.json', 'w') as fp:
            json.dump(data_params, fp)
