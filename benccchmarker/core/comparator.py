import importlib
import json
import os
from pathlib import Path
import platform
import subprocess
import yaml

import anndata
import pandas as pd

from benccchmarker.utils import convert_anndata_to_seurat, check_conda_env

class Comparator():

    FILE_KEYS = {

        # The following will be used as Comparator input
        'counts_file_path': 'counts_file_path',
        'metadata_file_path': 'metadata_file_path',
        'lr_file_path': 'lr_file_path',
        'microenvironment_file_path': 'microenvironment_file_path',
        'degs_file_path': 'degs_file_path',
        'protein_complex_file_path': 'protein_complex_file_path',
        'transcription_factor_file_path': 'transcription_factor_file_path',
        'anndata_file_path': 'anndata_file_path',
        'seurat_file_path': 'seurat_file_path',

        # Will be used to specify which methods are tried to be compared
        'comparator_methods_config_path': 'comparator_methods_config_path'

    }

    def __init__(
        self,
        output_file_path,
        species="hsapiens",
        num_nodes=4,
        **file_paths
    ):
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
        self.current_file = Path(__file__)

        self.file_paths = file_paths

        # Check if file paths are in FILE_KEYS raise error saying that the file path is not supported
        for key in self.file_paths:
            if key not in self.FILE_KEYS:
                raise ValueError(f"'{key}' is not supported. Supported file paths are {', '.join(list(self.FILE_KEYS.values()))}")

        self.species = species
        self.output_file_path = os.path.realpath(output_file_path)
        self.num_nodes = num_nodes

        self.comparator_run_info = {
            "methods": [],
            "start_time": None,
            "end_time": None,
            "species": self.species,
            "output_file_path": self.output_file_path,
            "platform": platform.platform(),
            "num_nodes": self.num_nodes
        }

        os.makedirs(self.output_file_path, exist_ok=True)

        if 'anndata_file_path' in file_paths and 'seurat_file_path' in file_paths:
            raise ValueError("Please provide either anndata or seurat_object file")

        if 'anndata_file_path' in file_paths:
            print("Loading anndata file and converting it to Seurat object...")
            self.file_paths['anndata_file_path'] = file_paths['anndata_file_path']

            self.comparator_run_info["input_file_path"] = file_paths['anndata_file_path']

            os.makedirs(f"{self.output_file_path}/tmp", exist_ok=True)
            seurat_file_path = f"{self.output_file_path}/tmp/seurat_object.rds"

            # Check if the anndata file contains "cell_type" obs
            anndata_file = anndata.read_h5ad(self.file_paths['anndata_file_path'])
            if "cell_type" not in anndata_file.obs.columns:
                raise ValueError("Please provide an anndata file with 'cell_type' in obs")

            convert_anndata_to_seurat(
                self.file_paths['anndata_file_path'],
                seurat_file_path
            )

            if os.path.exists(f'{seurat_file_path}'):
                print(f"Conversion complete. The Seurat object is saved as '{seurat_file_path}'")
            else:
                raise Exception("Conversion failed. Please check the output path.")

            self.file_paths["seurat_object"] = seurat_file_path

        elif 'seurat_file_path' in file_paths:
            self.file_paths["seurat_object"] = file_paths["seurat_object"]

            self.comparator_run_info["input_file_path"] = file_paths["seurat_object"]
        else:
            raise ValueError("Please provide either anndata or seurat_object file")

        # Open json file containing the methods dict json
        with open(self.current_file.parent.parent / "datasets" / "comparator_methods_dict.json", "r") as file:
            self.methods_dict = json.load(file)

        if 'comparator_methods_config_path' not in self.file_paths:
            print("No methods setup file is specified, will try to run the comparison on all of the methods\n")

            default_comparator_methods_config_path = self._generate_comparator_methods_config_path(self.current_file.parent.parent / "datasets" / "methods.yaml")

            self.comparator_methods_config_path = default_comparator_methods_config_path
        else:
            self.comparator_methods_config_path = self._generate_comparator_methods_config_path(self.file_paths["comparator_methods_config_path"])

    def _check_install(self):
        """
        Checks the availability of required Python and R packages.

        Parameters
        ----------
        comparator_methods_config_path : list
            A list of methods to compare.
        
        Returns
        -------
        list
            A list of installed packages, not installed packages, and not supported packages.
        """

        not_installed_packages = []
        not_supported_packages = []
        for method in self.comparator_methods_config_path:
            if method["name"] in self.methods_dict:
                if self.methods_dict[method["name"]]["language"] == "python":
                    python_package = self.methods_dict[method["name"]]["library"]
                    if python_package != None:
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

                        # if spec is None or spec == '':
                        #     not_installed_packages.append(method)
                    else:
                        print(f"{method['name']} can not be checked, because it's not installed as a package or a package name is not provided.")
                elif self.methods_dict[method["name"]]["language"] == "R":
                    r_package = self.methods_dict[method["name"]]["library"]

                    if method["environment"] == None:
                        spec = subprocess.run(
                                [f"Rscript", "-e", f"if (!requireNamespace('{package_name}', quietly = TRUE)) {{ quit(status = 1) }}"],
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE
                            ).returncode == 0
                    else:
                        spec = subprocess.run(
                                [f"{method['environment']}/bin/Rscript", "-e", f"if (!requireNamespace('{package_name}', quietly = TRUE)) {{ quit(status = 1) }}"],
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE
                            ).returncode == 0

                    # if spec == False:
                    #     not_installed_packages.append(method)
            else:
                # warnings.warn(f"Package {method} isn't supported. Currently only supporting some of the listed Python and R packages, to include more packages, please contact the maintainer.")
                not_supported_packages.append(method)

        installed_packages = list(set((tuple(method.items()) for method in self.comparator_methods_config_path)) - set(tuple(method.items()) for method in not_installed_packages+not_supported_packages))
        installed_packages = [dict(installed_package) for installed_package in installed_packages]

        return installed_packages, not_installed_packages, not_supported_packages
    
    def _generate_comparator_methods_config_path(self, comparator_methods_config_path_path):
        """
        Generates the methods to compare list.

        Parameters
        ----------
        comparator_methods_config_path_path : str
            Path to the methods to compare file.
        
        Returns
        -------
        list
            A list of methods to compare.
        """

        comparator_methods_config_path = yaml.safe_load(open(comparator_methods_config_path_path, "r"))["methods"]

        default_values = {
            "version": None,
            "environment": None,
            "scoring_method": None
        }

        for method in comparator_methods_config_path:
            for key, value in default_values.items():
                method.setdefault(key, value)

        return comparator_methods_config_path

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

                    # try:
                    if package["environment"] != None:
                        # try:
                        subprocess.run([f"{package['environment']}/bin/Rscript", "--vanilla", method_file_path, self.file_paths["seurat_object"], self.species, method_config_file_path, method_output_file_path, str(self.num_nodes)])

                        self.comparator_run_info["methods"].append(
                            {
                                "name": package["name"],
                                "language": "r",
                                "config_file": method_config_file_path,
                                "output_file": method_output_file_path,
                        
                            }
                        )
                        # except:
                        #     raise ValueError(f"Error in running the method {package['name']}")

                    else:
                        
                        # try:
                        subprocess.run(["Rscript", "--vanilla", method_file_path, self.file_paths["seurat_object"], self.species, method_config_file_path, method_output_file_path, str(self.num_nodes)])

                        self.comparator_run_info["methods"].append(
                            {
                                "name": package["name"],
                                "language": "r",
                                "config_file": method_config_file_path,
                                "output_file": method_output_file_path,
                        
                            }
                        )
                        # except:
                        #     raise ValueError(f"Error in running the method {package['name']}")
                    # except:
                    #     print(f"Error in running the method {package['name']}, this can be because of the method is not installed, misconfiguration of the method environment or the method is not compatible with the current dataset.")
                    #     print("Running the next method...")

                elif self.methods_dict[package["name"]]["language"] == "python":
                    method_file_path = f"{dir_path}/methods/{self.methods_dict[package['name']]['file']}.py"
                    method_config_file_path = f"{dir_path}/params/{self.methods_dict[package['name']]['file']}.yaml"
                    method_output_file_path = f"{self.output_file_path}/{self.methods_dict[package['name']]['file']}.csv"

                    # try:
                    if package["environment"] != None:
                        # try:
                        subprocess.run([f"{package['environment']}/bin/python", method_file_path, self.file_paths['anndata_file_path'], self.species, method_config_file_path, method_output_file_path, package['environment'], str(self.num_nodes)])

                        self.comparator_run_info["methods"].append(
                            {
                                "name": package["name"],
                                "language": "python",
                                "config_file": method_config_file_path,
                                "output_file": method_output_file_path,
                        
                            }
                        )
                        # except:
                        #     raise ValueError(f"Error in running the method {package['name']}")
                    else:
                        # try:
                        subprocess.run(["python", method_file_path, self.file_paths['anndata_file_path'], self.species, method_config_file_path, method_output_file_path, str(self.num_nodes)])

                        self.comparator_run_info["methods"].append(
                            {
                                "name": package["name"],
                                "language": "python",
                                "config_file": method_config_file_path,
                                "output_file": method_output_file_path,
                        
                            }
                        )
                        # except:
                        #     raise ValueError(f"Error in running the method {package['name']}")
                    # except:
                    #     print(f"Error in running the method {package['name']}, this can be because of the method is not installed, misconfiguration of the method environment or the method is not compatible with the current dataset.")
                    #     print("Running the next method...")

                    # Check if the output file is generated
                    if not os.path.exists(method_output_file_path):
                        raise ValueError(f"Error in running the method {package['name']}")
                    else:
                        print(f"Output file for {package['name']} is generated at {method_output_file_path}")

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

        self.comparator_run_info["start_time"] = pd.Timestamp.now()

        print("Checking the package availability...")

        installed_packages, not_installed_packages, not_supported_packages = self._check_install()

        print("You're trying to run the comparison for the following methods")
        
        for i, method in enumerate(sorted(self.comparator_methods_config_path, key=lambda x: x['name'])):
            print(f"{i+1}. {method['name']} ({method['version'] if method['version'] else 'default'})")

        if not_supported_packages:
            print("\nThe following packages are not supported in the current version of benCCChmarker")
            for i, method in enumerate(sorted(not_supported_packages)):
                print(f"{i+1}. {method['name']}")

        if not_installed_packages:
            print("\nThe following packages are not installed")
            for i, method in enumerate(sorted(not_installed_packages)):
                print(f"{i+1}. {method}")

        self._compare(installed_packages)

        print("Comparison is done!")

        self.comparator_run_info["end_time"] = pd.Timestamp.now()

        print("Generating the comparison report...")
        self._generate_report()

        print("Comparison report is generated at", f"{self.output_file_path}/comparator_run_summary.html")

    def get_comparator_run_info(self):
        """
        Returns the comparison results.

        Returns
        -------
        dict
            A dictionary containing the comparison results.
        """

        return self.comparator_run_info

    def _generate_report(self):
        """
        Generates the comparison report.

        Returns
        -------
        None
        """

        comparator_run_info = self.get_comparator_run_info()

        with open(self.current_file.parent.parent / "template" / "index.html", "r") as f:
            comparator_html_template = f.read()

        comparator_html_template = comparator_html_template.replace("$START_TIME", comparator_run_info["start_time"].strftime("%Y-%m-%d %H:%M:%S"))
        comparator_html_template = comparator_html_template.replace("$END_TIME", comparator_run_info["end_time"].strftime("%Y-%m-%d %H:%M:%S"))
        comparator_html_template = comparator_html_template.replace("$SPECIES", comparator_run_info["species"])
        comparator_html_template = comparator_html_template.replace("$PLATFORM", comparator_run_info["platform"])
        comparator_html_template = comparator_html_template.replace("$OUTPUT_PATH", comparator_run_info["output_file_path"])

        runtime_df = pd.DataFrame()
        inference_df = pd.DataFrame()

        for method in comparator_run_info["methods"]:
            method_runtime_df = pd.read_csv(method["output_file"].replace(".csv", "_met.csv"))
            method_inference_df = pd.read_csv(method["output_file"])

            method_inference_df = method_inference_df[["source", "ligand", "target", "receptor"]]
            method_inference_df.drop_duplicates(["source", "ligand", "target", "receptor"], inplace=True)
            method_inference_df["method"] = method["name"]

            method_runtime_df["method"] = method["name"]
            method_runtime_df["count_inferred_interactions"] = method_inference_df.shape[0]


            runtime_df = pd.concat([runtime_df, method_runtime_df])
            inference_df = pd.concat([inference_df, method_inference_df])

        method_runtimes = runtime_df.to_dict(orient="records")
        method_runtimes_html = ""

        for method_runtime in method_runtimes:
            method_name = method_runtime["method"]
            method_preprocessing_time = method_runtime["preprocessing_time"]
            method_inference_time = method_runtime["inference_time"]
            method_peak_ram_used = method_runtime["peak_ram_used_mib"]
            method_count_inferred_interactions = method_runtime["count_inferred_interactions"]

            method_runtimes_html += f"<tr><td>{method_name}</td><td>{method_preprocessing_time:.2f}</td><td>{method_inference_time:.2f}</td><td>{method_peak_ram_used:.2f}</td><td>{method_count_inferred_interactions}</td></tr>"

        inference_sltr_level_count_df = inference_df.groupby(["source", "ligand", "target", "receptor"]).count().reset_index().rename(columns={"method": "count"}).sort_values(by="count", ascending=False)
        inference_cell_level_count_df = inference_df.groupby(["source", "target"]).count().reset_index().sort_values(by="method", ascending=False)[["source", "target", "method"]].rename(columns={"method": "count"})
        inference_cell_level_count_df["count"] = inference_cell_level_count_df.groupby("source")["count"].transform(lambda x: x / x.max())

        top_10_inference_results = inference_sltr_level_count_df.head(10).to_dict(orient="records")

        top_10_inference_results_html = ""
        for top_10_inference_result in top_10_inference_results:
            source = top_10_inference_result["source"]
            ligand = top_10_inference_result["ligand"]
            target = top_10_inference_result["target"]
            receptor = top_10_inference_result["receptor"]
            count = top_10_inference_result["count"]

            top_10_inference_results_html += f"<tr><td>{source}</td><td>{ligand}</td><td>{target}</td><td>{receptor}</td><td>{count}</td></tr>"

        comparator_html_template = comparator_html_template.replace("$TOP_10_INFERENCE_RESULTS", top_10_inference_results_html)
        comparator_html_template = comparator_html_template.replace("$METHOD_RUNTIMES", method_runtimes_html)

        with open(f"{self.output_file_path}/comparator_run_summary.html", "w") as f:
            f.write(comparator_html_template)

    def add_method(self, method_name, method_language, method_library, method_file):
        """
        Adds a new method to the comparator.

        Parameters
        ----------
        method_name : str
            Name of the method.
        method_language : str
            Language of the method.
        method_library : str
            Library of the method.
        method_file : str
            File of the method.
        
        Returns
        -------
        None
        """ 
        self.methods_dict[method_name] = {
            "language": method_language,
            "library": method_library,
            "file": method_file
        }

        with open(self.current_file.parent.parent / "datasets" / "comparator_methods_dict.json", "w") as f:
            json.dump(self.methods_dict, f)

    def remove_method(self, method_name):
        """
        Removes a method from the comparator.

        Parameters
        ----------
        method_name : str
            Name of the method.
        
        Returns
        -------
        None
        """

        self.methods_dict.pop(method_name)

        with open(self.current_file.parent.parent / "datasets" / "comparator_methods_dict.json", "w") as f:
            json.dump(self.methods_dict, f)

    def update_method(self, method_name, method_language, method_library, method_file):
        """
        Updates a method in the comparator.

        Parameters
        ----------
        method_name : str
            Name of the method.
        method_language : str
            Language of the method.
        method_library : str
            Library of the method.
        method_file : str
            File of the method.
        
        Returns
        -------
        None
        """

        self.methods_dict[method_name] = {
            "language": method_language,
            "library": method_library,
            "file": method_file
        }

        with open(self.current_file.parent.parent / "datasets" / "comparator_methods_dict.json", "w") as f:
            json.dump(self.methods_dict, f)

    def get_method(self, method_name):
        """
        Returns the details of a method.

        Parameters
        ----------
        method_name : str
            Name of the method.
        
        Returns
        -------
        dict
            A dictionary containing the details of the method.
        """

        return self.methods_dict[method_name]
    
    def list_methods(self):
        print("The following methods are available for comparison:")
        for i, method in enumerate(self.methods_dict):
            print(f"{i+1}. {method}")
            
    def install_method(self, method_name):
        """
        Installs a method by creating a new environment and installing the required packages

        Parameters
        ----------
        method_name : str
            Name of the method.
        
        Returns
        -------
        None
        """

        if not check_conda_env(f"benccchmarker-{method_name}"):
            if method_name not in self.methods_dict:
                raise ValueError(f"{method_name} is not available in the current version of benCCChmarker")

            installation_script_path = self.current_file.parent.parent / "install" / f"benccchmarker-{method_name}.sh"

            try:
                subprocess.run(["bash", installation_script_path])

                if check_conda_env(f"benccchmarker-{method_name}"):
                    print(f"{method_name} is successfully installed in 'benccchmarker-{method_name}' environment")

                    # Update self.comparator_methods_config_path
                    self.comparator_methods_config_path.append(
                        {
                            "name": method_name,
                            "version": None,
                            "environment": f"{os.environ['CONDA_PREFIX']}/envs/benccchmarker-{method_name}",
                            "scoring_method": None
                        }
                    )

                else:
                    raise ValueError(f"Error occurred while installing {method_name}")
            except subprocess.CalledProcessError as e:
                print(f"Error occurred while installing {method_name}: {e}")
        else:
            print(f"{os.environ['CONDA_PREFIX']}/envs/benccchmarker-{method_name}")
            print(f"{method_name} is already installed in 'benccchmarker-{method_name}' environment")
