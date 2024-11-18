import importlib
import subprocess
import os
from pathlib import Path
import yaml

import anndata
import pandas as pd
import rpy2.robjects.packages as rpackages

class Comparator():
    """
    Comparator allows to run and compare the result of different algorithms on the 
    user provided dataset.

    Parameters
    ----------
    output_file_path : str
        Path to the output file.
    species : str
        Species of the dataset.
    file_paths : dict
        A dictionary containing the paths to the files required for the comparison.
    """

    FILE_KEYS = {

        # The following will be used as Comparator input
        'counts': 'counts_file_path',
        'metadata': 'metadata_file_path',
        'lr': 'lr_file_path',
        'microenvironment': 'microenvironment_file_path',
        'degs': 'degs_file_path',
        'protein_complex': 'protein_complex_file_path',
        'transcription_factor': 'transcription_factor_file_path',
        'anndata': 'anndata_file_path',
        'seurat_object': 'seurat_object_file_path',

        # Will be used to specify which methods are tried to be compared
        'methods_to_compare': 'methods_to_compare_path'

    }

    def __init__(
        self,
        output_file_path,
        species="hsapiens",
        **file_paths
    ):
        self.current_file = Path(__file__)

        self.file_paths = file_paths
        self.species = species
        self.output_file_path = output_file_path

        if 'anndata' in file_paths:
            self.file_paths["anndata"] = file_paths["anndata"]

            self.adata = anndata.read_h5ad(file_paths["anndata"])
        elif 'seurat_object' in file_paths:
            self.file_paths["seurat_object"] = file_paths["seurat_object"]
        else:
            raise ValueError("Please provide either anndata or seurat_object file")


        # Open json file containing the methods dict json
        with open(self.current_file.parent.parent / "datasets" / "comparator_methods_dict.json", "r") as file:
            self.methods_dict = json.load(file)

        if 'methods_to_compare' not in self.file_paths:
            print("No methods setup file is specified, will try to run the comparison on all of the methods\n")

            default_methods_to_compare = self._generate_methods_to_compare(self.current_file.parent.parent / "datasets" / "methods.yaml")

            self.methods_to_compare = default_methods_to_compare
        else:
            self.methods_to_compare = self._generate_methods_to_compare(self.file_paths["methods_to_compare"])

    def _check_install(self):
        """
        Checks the availability of required Python and R packages.

        Parameters
        ----------
        methods_to_compare : list
            A list of methods to compare.
        
        Returns
        -------
        list
            A list of installed packages, not installed packages, and not supported packages.
        """

        not_installed_packages = []
        not_supported_packages = []
        for method in self.methods_to_compare:
            if method["name"] in self.methods_dict:
                if self.methods_dict[method["name"]]["language"] == "Python":
                    python_package = self.methods_dict[method["name"]]["library"]

                    if method["environment"] == None:
                        spec = importlib.util.find_spec(python_package)
                    else:
                        spec = subprocess.run(
                                f"{method['environment']}/bin/pip show {method['name']}",
                                shell=True,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE,
                                text=True
                            ).stdout

                    if spec is None or spec == '':
                        not_installed_packages.append(method)
                elif self.methods_dict[method["name"]]["language"] == "R":
                    r_package = self.methods_dict[method["name"]]["library"]

                    if method["environment"] == None:
                        spec = rpackages.isinstalled(r_package)
                    else:
                        spec = subprocess.run(
                                [f"{method['environment']}/bin/Rscript", "-e", f"if (!requireNamespace('{package_name}', quietly = TRUE)) {{ quit(status = 1) }}"],
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE
                            ).returncode == 0

                    if spec == False:
                        not_installed_packages.append(method)
            else:
                # warnings.warn(f"Package {method} isn't supported. Currently only supporting some of the listed Python and R packages, to include more packages, please contact the maintainer.")
                not_supported_packages.append(method)

        installed_packages = list(set((tuple(method.items()) for method in self.methods_to_compare)) - set(tuple(method.items()) for method in not_installed_packages+not_supported_packages))
        installed_packages = [dict(installed_package) for installed_package in installed_packages]

        return installed_packages, not_installed_packages, not_supported_packages
    
    def _generate_methods_to_compare(self, methods_to_compare_path):
        """
        Generates the methods to compare list.

        Parameters
        ----------
        methods_to_compare_path : str
            Path to the methods to compare file.
        
        Returns
        -------
        list
            A list of methods to compare.
        """

        methods_to_compare = yaml.safe_load(open(methods_to_compare_path, "r"))["methods"]

        default_values = {
            "version": None,
            "environment": None,
            "scoring_method": None
        }

        for method in methods_to_compare:
            for key, value in default_values.items():
                method.setdefault(key, value)

        return methods_to_compare

    def _compare(self, installed_packages):
        """
        Performs the comparison among installed packages.

        Parameters
        ----------
        installed_packages : list
            A list of installed
        
        Returns
        -------
        None
        """

        # Current file path
        current_path = os.path.realpath(__file__)

        dir_path = os.path.dirname(os.path.dirname(current_path))
        if not os.path.exists(self.output_file_path):
            os.makedirs(self.output_file_path)

        for i, package in enumerate(installed_packages):
            print(f"Running calculation for {package['name']} ({i+1}/{len(installed_packages)})")

            if package["name"] in self.methods_dict:

                if self.methods_dict[package["name"]]["language"] == "r":
                    method_file_path = f"{dir_path}/methods/{self.methods_dict[package['name']]['file']}.r"
                    method_config_file_path = f"{dir_path}/params/{self.methods_dict[package['name']]['file']}.yaml"
                    method_output_file_path = f"{self.output_file_path}/{self.methods_dict[package['name']]['file']}.csv"

                    if package["environment"] != None:
                        try:
                            subprocess.run([f"{package['environment']}/bin/Rscript", "--vanilla", method_file_path, self.file_paths["seurat_object"], self.species, method_config_file_path, method_output_file_path])
                        except:
                            raise ValueError(f"Error in running the method {package['name']}")

                    else:
                        subprocess.run(["Rscript", "--vanilla", method_file_path, self.file_paths["seurat_object"], self.species, method_config_file_path, method_output_file_path])
                elif self.methods_dict[package["name"]]["language"] == "py":
                    method_file_path = f"{dir_path}/methods/{self.methods_dict[package['name']]['file']}.py"
                    method_config_file_path = f"{dir_path}/params/{self.methods_dict[package['name']]['file']}.yaml"
                    method_output_file_path = f"{self.output_file_path}/{self.methods_dict[package['name']]['file']}.csv"

                    if package["environment"] != None:
                        try:
                            subprocess.run([f"{package['environment']}/bin/python", method_file_path, self.file_paths["anndata"], self.species, method_config_file_path, method_output_file_path])
                        except:
                            raise ValueError(f"Error in running the method {package['name']}")

                    else:
                        subprocess.run(["python", method_file_path, self.file_paths["anndata"], self.species, method_config_file_path, method_output_file_path])

    def run(self, force_yes=False):
        """
        Initiates the comparison process.

        Parameters
        ----------
        force_yes : bool
            A flag to force the comparison process regardless of the package availability.
        
        Returns
        -------
        None
        """
        print("Checking the package availability...")

        installed_packages, not_installed_packages, not_supported_packages = self._check_install()

        print("You're trying to run the comparison for the following methods")
        for i, method in enumerate(sorted(self.methods_to_compare)):
            print(f"{i+1}. {method['name']} ({method['version'] if method['version'] else 'default'})")

        if not_supported_packages:
            print("\nThe following packages are not supported in the current version of benCCChmarker")
            for i, method in enumerate(sorted(not_supported_packages)):
                print(f"{i+1}. {method}")

        if not_installed_packages:
            print("\nThe following packages are not installed")
            for i, method in enumerate(sorted(not_installed_packages)):
                print(f"{i+1}. {method}")

        self._compare(installed_packages)
        
