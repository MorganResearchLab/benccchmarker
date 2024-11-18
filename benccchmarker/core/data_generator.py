import json
import logging
import os
from pathlib import Path
import warnings

warnings.filterwarnings("ignore")

import anndata
import numpy as np
import pandas as pd
import scanpy as sc
import scipy
from scipy.stats import nbinom

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

    def _check_distribution_parameters(self, mean_expression, dispersion):
        """
        Check if the mean expression and dispersion parameters are valid.
        Raise exception if mean_expression / (dispersion ** 2) and (mean_expression ** 2) / ((dispersion ** 2) - mean_expression) < 0

        Parameters
        ----------
        mean_expression : float
            Mean expression level.
        dispersion : float
            Dispersion parameter.
        
        Returns
        -------
        bool
            True if the mean expression and dispersion parameters are valid.
        """
        if mean_expression / (dispersion ** 2) < 0 or (mean_expression ** 2) / ((dispersion ** 2) - mean_expression) < 0:
            logging.error("Invalid mean expression and dispersion parameters, mean_expression / (dispersion ** 2) and (mean_expression ** 2) / ((dispersion ** 2) - mean_expression) has to be greater than 0.")
            raise ValueError("Invalid mean expression and dispersion parameters, mean_expression / (dispersion ** 2) and (mean_expression ** 2) / ((dispersion ** 2) - mean_expression) has to be greater than 0.")
        return True

    def _generate_lr_pairs(self, lr_pair_df, num_lr_pairs, gene_names=None):
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
        
        cell_type_weights = np.random.randint(1, 10, num_cell_types)
        cell_type_probas = [cell_type_weight/sum(cell_type_weights) for cell_type_weight in cell_type_weights]
        cell_type_labels = np.random.choice(range(num_cell_types), size=num_cells, p=cell_type_probas)
        cell_type_labels = [f"CellType{i}" for i in cell_type_labels]

        return cell_type_labels

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
        interaction_maps = interaction_maps[np.argwhere(interaction_map_probas > 0.5).T[0]]

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

    def simulate_denovo(
        self,
        num_cells=5000,
        num_genes=500,
        num_lr_pairs=10,
        num_cell_types=5,
        method="overexpression",
        overexpression_scale=2,
        mean_expression=100,
        dispersion=30,
        num_batch=1,
        batch_factor=0.5,
        batch_mean=0,
        batch_sd=0.1,
        be_method="shift_by_mean",
        differential_interaction=False,
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
        overexpression_scale : float
            Scale factor for overexpressing ligand-receptor pairs. Required if method is "overexpression".
        mean_expression : float
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
        if num_lr_pairs > num_genes:
            logging.error("num_lr_pairs has to be lower than num_genes.")
            raise ValueError("num_lr_pairs has to be lower than num_genes.")

        if not os.path.exists(output_path):
            os.makedirs(output_path)

        # Check if the mean expression and dispersion parameters are valid
        self._check_distribution_parameters(mean_expression, dispersion)

        # Sample ligand-receptor pairs based on the number of ligand-receptor pairs input
        lr_pairs = self._generate_lr_pairs(self.lr_pair_df, num_lr_pairs)

        lr_genes = np.concatenate([list(lr_pairs.keys()), list(lr_pairs.values())])
        num_lr_genes = len(np.unique(lr_genes))

        cell_types = [f"CellType{i}" for i in range(num_cell_types)]
        cell_labels = [f"cell_{i}" for i in range(num_cells)]
        
        cell_type_labels = self._generate_cell_type_labels(num_cells, num_cell_types)

        # Generate source target pairs
        source_target_pairs = self._generate_source_target_pairs(num_cell_types)

        # Generte combination of cell types and ligand-receptor pairs
        cell_type_receptor_ligand_dict = self._generate_cell_type_receptor_ligand_dict(cell_types, lr_pairs)

        # Generate interaction maps
        interaction_maps = self._generate_interaction_maps(cell_type_receptor_ligand_dict, lr_pairs)
        print(f"Generated {len(interaction_maps)} interactions")
        
        # Generate anndata
        expression_levels = np.zeros((num_cells, num_lr_genes))

        gene_labels = np.unique(np.concatenate([list(lr_pairs.keys()), list(lr_pairs.values())]))
        
        simulated_df = pd.DataFrame(expression_levels, columns=gene_labels, index=cell_type_labels)

        if scenario == "lherhe":
            for interaction_map in interaction_maps:
                simulated_df.loc[interaction_map["source"], interaction_map["ligand"]] = -999
                simulated_df.loc[interaction_map["target"], interaction_map["receptor"]] = -999
        elif scenario == "lhemhe":
            for interaction_map in interaction_maps:
                simulated_df.loc[interaction_map["source"], interaction_map["ligand"]] = -999
        else:
            for interaction_map in interaction_maps:
                simulated_df.loc[interaction_map["target"], interaction_map["receptor"]] = -999

        # IGNORE THIS PART OF THE CODE FOR NOW WILL HAVE IT IN THE NEARBY FUTURE RELEASES
        modified_simulated_df = simulated_df.copy()
        condition_labels = np.random.choice(range(1), len(cell_labels))

        # if differential_interaction:

        #     condition_labels = np.random.choice(range(2), len(cell_labels))
        #     simulated_df["condition"] = condition_labels

        #     non_differential_interaction_map_probas = np.random.uniform(0, 1, len(interaction_maps))
        #     non_differential_interaction_maps = interaction_maps[np.argwhere(non_differential_interaction_map_probas > 0.5).T[0]]

        #     differential_interaction_maps = interaction_maps[np.argwhere(non_differential_interaction_map_probas <= 0.5).T[0]]

        #     # Remove element from non_differential_interaction_maps if both the source and ligand are in the differential_interaction_maps or if both the target and receptor are in the differential_interaction_maps

        #     print(f"{len(non_differential_interaction_maps)} which are differential interactions")

        #     print(non_differential_interaction_maps)

        #     non_differential_interaction_maps = [interaction_map for interaction_map in non_differential_interaction_maps if (interaction_map["source"] not in [interaction["source"] for interaction in differential_interaction_maps] and interaction_map["ligand"] not in [interaction["ligand"] for interaction in differential_interaction_maps]) or (interaction_map["target"] not in [interaction["target"] for interaction in differential_interaction_maps] and interaction_map["receptor"] not in [interaction["receptor"] for interaction in differential_interaction_maps])]

        #     print(f"{len(non_differential_interaction_maps)} which are non differential interactions")

        #     # Save differential_interaction_maps into a .json file
        #     with open(f"{output_path}/{output_filename}_differential_interaction_maps.json", "w") as f:
        #         json.dump(differential_interaction_maps.tolist(), f)

        #     with open(f"{output_path}/{output_filename}_non_differential_interaction_maps.json", "w") as f:
        #         json.dump(non_differential_interaction_maps.tolist(), f)

        #     if scenario == "lherhe":
        #         for interaction_map in differential_interaction_maps:
        #             modified_simulated_df.loc[(simulated_df.index == interaction_map["source"]) & ( simulated_df["condition"] == 0), interaction_map["ligand"]] = 0
        #             modified_simulated_df.loc[(simulated_df.index == interaction_map["target"]) & ( simulated_df["condition"] == 1), interaction_map["receptor"]] = 0
        #     elif scenario == "lhemhe":
        #         for interaction_map in differential_interaction_maps:
        #             modified_simulated_df.loc[(simulated_df.index == interaction_map["source"]) & ( simulated_df["condition"] == 0), interaction_map["ligand"]] = 0
        #     else:
        #         for interaction_map in differential_interaction_maps:
        #             modified_simulated_df.loc[(simulated_df.index == interaction_map["target"]) & ( simulated_df["condition"] == 1), interaction_map["receptor"]] = 0
        
        n, p = get_n_p_from_mean_dispersion(mean_expression, dispersion)
        modified_simulated_df = modified_simulated_df.map(lambda x: nbinom.rvs(n, p) if x == 0 else x)

        mean_altered_expression = mean_expression * overexpression_scale
        n, p = get_n_p_from_mean_dispersion(mean_altered_expression, dispersion)

        modified_simulated_df = modified_simulated_df.map(lambda x: nbinom.rvs(n, p) if x == -999 else x)

        modified_simulated_df["condition"] = condition_labels
        simulated_df = modified_simulated_df.copy()        
        simulated_df.drop(columns=["condition"], inplace=True)

        if add_control:
            non_lr_gene_path = self.current_file.parent.parent / "datasets" / "non_lr_genes_sample.csv"
            non_lr_gene_names_df = pd.read_csv(non_lr_gene_path)

            non_lr_genes = non_lr_gene_names_df.sample(num_genes - num_lr_genes, random_state=self.seed)["gene"].tolist()

            markers = non_lr_genes[:num_cell_types]
            cell_type_markers = {}

            for marker in markers:
                simulated_df[marker] = np.nan

            for i, cell_type in enumerate(cell_types):
                cell_type_markers[cell_type] = markers[i]
            
            # Save cell_type_markers into a .json file
            with open(f"{output_path}/{output_filename}_cell_type_markers.json", "w") as f:
                json.dump(cell_type_markers, f)

            for cell_type in cell_types:
                n, p = get_n_p_from_mean_dispersion(mean_altered_expression, dispersion)
                simulated_df.loc[cell_type, cell_type_markers[cell_type]] = nbinom.rvs(n, p, size=len(simulated_df.loc[cell_type, cell_type_markers[cell_type]]))

            n,p = get_n_p_from_mean_dispersion(mean_expression, dispersion)
            simulated_df = simulated_df.map(lambda x: nbinom.rvs(n, p) if pd.isna(x) else x)

            # create new columns from non_lr_genes[num_cell_types:] and append it to the simulated_df with value of nbinom.rvs(n, p)
            for gene in non_lr_genes[num_cell_types:]:
                simulated_df[gene] = nbinom.rvs(n, p, size=len(simulated_df))            

        if num_batch > 1:
            cell_types = simulated_df.index

            # Generate batch based on num_batch
            batch = np.random.choice(range(num_batch), len(cell_types))

            modified_simulated_df = simulated_df.copy()

            simulated_df["batch"] = batch

            if be_method == "shift_by_sd":
                random_factors = np.random.uniform(0, 1, num_batch) * dispersion
            
            if be_method == "shift_by_mean":
                random_factors = np.random.uniform(0, 1, num_batch) * mean_expression

            random_factors = [int(factor) for factor in random_factors]

            # Save random_factors to a file
            with open(f"{output_path}/{output_filename}_random_factors_for_be.json", "w") as f:
                json.dump(random_factors, f)

            modified_simulated_dfs = []

            for i in range(num_batch):

                # Multiplies the expression values of each batch by a factor (scaling)
                if be_method == "scale":
                    modified_simulated_df.loc[simulated_df["batch"] == i] = modified_simulated_df.loc[simulated_df["batch"] == i].astype(float) * (i + 1) * batch_factor

                # Add by a random number sampled from normal distribution (shifting)
                # Changing sd
                if be_method == "shift_by_sd":

                    random_vector = np.random.normal(batch_mean, random_factors[i], size=(modified_simulated_df.loc[simulated_df["batch"] == i].shape[1]))

                    # Save random_vector to a file
                    with open(f"{output_path}/{output_filename}_random_vector_for_be.json", "w") as f:
                        json.dump(random_vector.tolist(), f)
                    
                    # Convert random_vector to int
                    random_vector = [int(factor) for factor in random_vector]

                    random_matrix = np.tile(random_vector, (modified_simulated_df.loc[simulated_df["batch"] == i].shape[0], 1))

                    modified_simulated_df.loc[simulated_df["batch"] == i] += random_matrix
                    
                # Changing mean
                if be_method == "shift_by_mean":
                    
                    random_vector = np.random.normal(random_factors[i] , batch_sd, size=(modified_simulated_df.loc[simulated_df["batch"] == i].shape[1]))

                    # Save random_vector to a file
                    with open(f"{output_path}/{output_filename}_random_vector_for_be.json", "w") as f:
                        json.dump(random_vector.tolist(), f)
                    
                    # Convert random_vector to int
                    random_vector = [int(factor) for factor in random_vector]

                    random_matrix = np.tile(random_vector, (modified_simulated_df.loc[simulated_df["batch"] == i].shape[0], 1))

                    modified_simulated_df.loc[simulated_df["batch"] == i] += random_matrix

                # Add a constant value (shifting non random)
                if be_method == "shift_by_constant":
                    modified_simulated_df.loc[simulated_df["batch"] == i] = simulated_df.loc[simulated_df["batch"] == i].astype(float) + (i + 1) * batch_factor
            
            simulated_df.to_csv(f"{output_path}/{output_filename}_batches.csv")

            batch_labels = simulated_df["batch"].to_list()
            batch_label_df = pd.DataFrame(simulated_df["batch"])
            batch_label_df.index = cell_labels

            batch_label_df.to_csv(f"{output_path}/{output_filename}_batch_labels.csv")
            simulated_df = modified_simulated_df.copy()

            # Convert any negative values to 0
            simulated_df[simulated_df < 0] = 0

        # Drop batch column
        simulated_df.index = cell_labels

        adata = anndata.AnnData(simulated_df)
        adata.obs["cell_type"] = cell_type_labels

        if num_batch > 1:
            adata.obs["batch"] = batch_labels
        else:
            adata.obs["batch"] = 0

        if differential_interaction:
            adata.obs["condition"] = condition_labels
        else:
            adata.obs["condition"] = 0
            
        adata.layers['counts'] = adata.X.copy()

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
            "method": method,
            "overexpression_scale": overexpression_scale,
            "mean_expression": mean_expression,
            "dispersion": dispersion,
            "num_batch": num_batch,
            "batch_factor": batch_factor,
            "batch_mean": batch_mean,
            "batch_sd": batch_sd,
            "be_method": be_method,
            "differential_interaction": differential_interaction,
            "output_path": output_path,
            "output_filename": output_filename,
            "num_batch": num_batch,
            "scenario": scenario,
            "add_control": add_control
        }

        with open(f'{output_path}/data_params.json', 'w') as fp:
            json.dump(data_params, fp)

    def simulate_from_reference(
        self,
        adata,
        cell_type_labels="cell_type",
        num_cells=None,
        num_lr_pairs=10,
        overexpression_scale=1.25,
        output_path="out",
        output_filename="synthetic_data",
        num_batch=1,
        batch_mean=0,
        batch_sd=0.1,
        batch_factor=0.5,
        be_method="shift_by_mean",
        scenario="default"
    ):
        """
        Synthesise single-cell RNA-seq data.
        
        Parameters
        ----------
        adata : anndata.AnnData
            Anndata object.
        num_cells : int
            Number of cells to synthesize.
        num_lr_pairs : int
            Number of ligand-receptor pairs to synthesize.
        overexpression_scale : float
            Scale factor for overexpressing ligand-receptor pairs.
        output_path : str
            Path to save the output file.
        output_filename : str
            Name of the output file.
        scenario : str
            Name of the scenario to synthesize the gene expression. 
            "default" will synthesize the gene expression based on the ligand-receptor pairs in the anndata object.
        """


        # if scenario == "remove_hvg":
        #     sc.pp.highly_variable_genes(adata)
        #     adata = adata[:, ~adata.var.highly_variable]
            
        # Find ligand and receptor genes in the anndata object
        lr_genes_in_matrix = self._find_lr_genes_intersection(adata)
        print(f"Found {len(lr_genes_in_matrix)} ligand-receptor pairs in anndata")

        gene_names = adata.var_names.to_list()
        
        # Identify ligand and receptor gene indices
        lr_flag_matrix = np.array([gene in lr_genes_in_matrix for gene in gene_names])
        non_lr_gene_indices = np.where(~lr_flag_matrix)[0]
        lr_gene_indices = np.where(lr_flag_matrix)[0]

        np.random.shuffle(non_lr_gene_indices)
        new_lr_gene_indices = non_lr_gene_indices[:len(lr_gene_indices)]
        index_mapping = dict(zip(lr_gene_indices, new_lr_gene_indices))

        for lr_idx, new_idx in index_mapping.items():
            gene_names[lr_idx], gene_names[new_idx] = gene_names[new_idx], gene_names[lr_idx]
            
        adata.var_names = gene_names

        lr_pairs = self._generate_lr_pairs(self.lr_pair_df, num_lr_pairs, gene_names)

        print(f"Found {len(lr_pairs)} eligible ligand-receptor out of {num_lr_pairs} pairs")

        num_cell_types = adata.obs["cell_type"].nunique()

        cell_types = adata.obs["cell_type"].unique()
        cell_type_labels = adata.obs.cell_type.to_list()

        # Generate source target pairs
        source_target_pairs = self._generate_source_target_pairs(num_cell_types)

         # Generte combination of cell types and ligand-receptor pairs
        cell_type_receptor_ligand_dict = self._generate_cell_type_receptor_ligand_dict(cell_types, lr_pairs)

        # Generate interaction maps
        interaction_maps = self._generate_interaction_maps(cell_type_receptor_ligand_dict, lr_pairs)
        print(f"Generated {len(interaction_maps)} interactions")

        synth_adata = adata.copy()
        unique_row_col_index_pairs = set()

        synthetic_gene_indices = set([])
        ineligible_ground_truth_genes = set([])

        for interaction_map in interaction_maps:
            
            source_cells = synth_adata.obs.index[synth_adata.obs["cell_type"] == interaction_map["source"]]
            target_cells = synth_adata.obs.index[synth_adata.obs["cell_type"] == interaction_map["target"]]
            ligand_genes = synth_adata.var.index[synth_adata.var.index == interaction_map["ligand"]]
            receptor_genes = synth_adata.var.index[synth_adata.var.index == interaction_map["receptor"]]

            source_indices = synth_adata.obs.index.get_indexer(source_cells)
            target_indices = synth_adata.obs.index.get_indexer(target_cells)
            ligand_indices = synth_adata.var.index.get_indexer(ligand_genes)
            receptor_indices = synth_adata.var.index.get_indexer(receptor_genes)

            synthetic_gene_indices.add((ligand_genes[0], ligand_indices[0]))
            synthetic_gene_indices.add((receptor_genes[0], receptor_indices[0]))

            for row in source_indices:
                for col in ligand_indices:
                    unique_row_col_index_pairs.add((row, col))

            for row in target_indices:
                for col in receptor_indices:
                    unique_row_col_index_pairs.add((row, col))

        synth_adata_lil = synth_adata.X.tolil()

        for row, col in unique_row_col_index_pairs:
            # There are some values that are having 0 expression so overexpressing them will result in 0 expression
            synth_adata_lil[row, col] *= overexpression_scale

        # Find all the columns where column_name is in lr_genes_in_matrix
        lr_genes_in_matrix = sorted(list(set(lr_genes_in_matrix)))
        lr_gene_indices = [synth_adata.var.index.get_loc(gene) for gene in lr_genes_in_matrix]
        
        # Zero out the expression values for the lr_gene_indices if the row and col index pair is not in unique_row_col_index_pairs
        not_in_col = 0

        for col in lr_gene_indices:
            if col not in [col for _, col in unique_row_col_index_pairs]:
                synth_adata_lil[:, col] = 0
                not_in_col += 1
                
        synth_adata.X = synth_adata_lil.tocsr()

        synthetic_gene_indices = sorted(list(synthetic_gene_indices))
        for gene, idx in synthetic_gene_indices:
            if synth_adata_lil[:, idx].mean() == 0:
                ineligible_ground_truth_genes.add(gene)


        # Get all of the positions where the value is 0 in synth_adata.X
        zero_indices = np.argwhere(synth_adata.X.toarray() == 0)

        print(f"There are {len(zero_indices)} expression values with 0 expression out of {synth_adata.X.shape[0] * synth_adata.X.shape[1]}")

        # Add new columns to synth_adata called batch where it splits the data into num_batch batches based on the cell_type proportionally
        if num_batch > 1:
            cell_type_counts = synth_adata.obs['cell_type'].value_counts()
            batch_labels = np.zeros(synth_adata.shape[0], dtype=int)

            for cell_type, count in cell_type_counts.items():
                indices = synth_adata.obs.index[synth_adata.obs['cell_type'] == cell_type].tolist()
                int_indices = [synth_adata.obs.index.get_loc(idx) for idx in indices]
                batch_labels_for_cell_type = np.random.choice(range(num_batch), size=count, replace=True)
                batch_labels[int_indices] = batch_labels_for_cell_type

            synth_adata.obs['batch'] = batch_labels

            # random_factors = np.random.uniform(0, 1, num_batch)
            random_factors = np.random.uniform(0, 9, num_batch)

            for i in range(num_batch):
                # # Scaling
                if be_method == "scale":
                    synth_adata.X[synth_adata.obs["batch"] == i] *= (i + 1) * batch_factor

                # # Shifting
                # Changing sd
                if be_method == "shift_by_sd":
                    random_vector = np.random.normal(batch_mean , random_factors[i], size=(synth_adata.X[synth_adata.obs["batch"] == i].shape[1]))

                    random_matrix = np.tile(random_vector, (synth_adata.X[synth_adata.obs["batch"] == i].shape[0], 1))

                    synth_adata.X[synth_adata.obs["batch"] == i] += random_matrix


                # # Add a constant value (shifting non random)
                if be_method == "shift_by_constant":
                    if scipy.sparse.issparse(synth_adata.X):
                        synth_adata.X = synth_adata.X.toarray() 
                        synth_adata.X[synth_adata.obs["batch"] == i] += (i + 1) * batch_factor
                        synth_adata.X = scipy.sparse.csr_matrix(synth_adata.X)
                # synth_adata.X[synth_adata.obs["batch"] == i] += (i + 1) * batch_factor
        
        # Update synth_adata.X to 0 using the zero_indices
        synth_adata.X[zero_indices[:, 0], zero_indices[:, 1]] = 0

        # New ground truth
        interaction_maps = np.array([interaction_map for interaction_map in interaction_maps if interaction_map["ligand"] not in ineligible_ground_truth_genes and interaction_map["receptor"] not in ineligible_ground_truth_genes])

        if not os.path.exists(output_path):
            os.makedirs(output_path)
    
        with open(f'{output_path}/synthetic_ground_truth.json', 'w') as fp:
            json.dump(interaction_maps.tolist(), fp)

        synth_adata.write_h5ad(f'{output_path}/{output_filename}.h5ad')

        self.synthetic_data = synth_adata

        ground_truth_df = self._generate_ground_truth(interaction_maps)
        ground_truth_df.to_csv(f"{output_path}/{output_filename}_ground_truth.csv", index=False)

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
            "reference_adata": adata.obs["cell_type"].nunique(),
            "num_cells": num_cells,
            "num_lr_pairs": num_lr_pairs,
            "overexpression_scale": overexpression_scale,
            "output_path": output_path,
            "output_filename": output_filename,
            "num_batch": num_batch,
            "batch_mean": batch_mean,
            "batch_sd": batch_sd,
            "batch_factor": batch_factor,
            "be_method": be_method,
            "scenario": scenario
        }

        with open(f'{output_path}/data_params.json', 'w') as fp:
            json.dump(data_params, fp)

        print(f"Conversion complete. The Seurat object is saved as '{output_path}/{output_filename}_seurat.rds'.")
