import os
import subprocess
import yaml
from pathlib import Path

from anndata import read_h5ad
import pandas as pd
import scanpy as sc

def _read_h5ad(file_path: str) -> pd.DataFrame:
    """
    Reads an AnnData file (.h5ad) and converts it to a pandas DataFrame.

    Args:
        file_path (str): Path to the AnnData file.

    Returns:
        pd.DataFrame: DataFrame containing the AnnData.
    """
    adata = read_h5ad(file_path)
    df = adata.to_df().T
    return df

def _read_h5(file_path: str) -> pd.DataFrame:
    """
    Reads an HDF5 file (.h5) and returns a pandas DataFrame.

    Args:
        file_path (str): Path to the HDF5 file.

    Returns:
        pd.DataFrame: DataFrame containing the data.
    """
    df = pd.read_hdf(file_path)
    return df

def _read_dir(dir_path: str) -> pd.DataFrame:
    """
    Placeholder function for reading data from a directory.
    Currently not implemented.

    Args:
        dir_path (str): Path to the directory.

    Returns:
        pd.DataFrame: None.
    """
    pass

def read_file(
    path: str,
) -> pd.DataFrame:
    """
    Reads various types of files and returns a DataFrame.

    Args:
        path (str): Path to the file.

    Returns:
        pd.DataFrame: DataFrame containing the data from the file.

    Raises:
        FileNotFoundError: If the file is not found.
    """
    file_extension = os.path.splitext(path)[-1]

    if os.path.isfile(path):
        if file_extension == '.h5ad':
            return _read_h5ad(path)
        if file_extension == '.h5':
            return _read_h5(path)
        if file_extension == '.csv':
            try:
                return sc.read_csv(path)
            except:
                return pd.read_csv(path)
        if file_extension == '.tsv':
            return pd.read_csv(path, sep='\t')
        if file_extension == '.yaml':
            stream = open(path, "r")
            methods_setup = yaml.load(stream, Loader=yaml.Loader)

            return methods_setup
    else:
        raise FileNotFoundError

def read_file_and_append(
    path:str,
    main_file_list: dict,
    key: str
) -> None:
    """
    Reads a file and appends it to a dictionary under the given key.

    Args:
        path (str): Path to the file.
        main_file_list (dict): Dictionary to which the file will be appended.
        key (str): Key under which the file will be appended.
    """
    main_file_list[key] = read_file(path)

def load_config(config_file):
    """
    Load a configuration file in YAML format.

    Args:
        config_file (str): Path to the configuration file.
    """
    with open(config_file, 'r') as file:
        config = yaml.load(file, Loader=yaml.FullLoader)
    return config

def convert_anndata_to_seurat(adata_path, output_dir):
    """
    Convert an AnnData object to Seurat format.

    Args:
        adata_path (str): Path to the AnnData object.
        output_dir (str): Path to the output directory.
    """

    current_file = Path(__file__)
    converter_path = current_file.parent / "preprocess" / "convert_anndata_to_seurat.r"
    cmd = f'Rscript --vanilla {converter_path} {adata_path} {output_dir}'
    subprocess.run(cmd, shell=True, capture_output=True)

def get_n_p_from_mean_dispersion(mean_expression, dispersion):
    """
    Calculate the n and p parameters for the negative binomial distribution from the mean expression and dispersion.

    Parameters
    ----------
    mean_expression : float
        Mean expression level.
    dispersion : float
        Dispersion parameter.

    Returns
    -------
    tuple
        Tuple containing the n and p parameters.
    """
    p = mean_expression / (dispersion ** 2)
    n = (mean_expression ** 2) / ((dispersion ** 2) - mean_expression)

    return n, p