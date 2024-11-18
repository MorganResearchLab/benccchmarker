import os
import yaml

from benccchmarker._utils import read_file_and_append

class FileLoader:
    """
    A base class for loading and storing data from files.

    Attributes
    ----------
    counts_file_path : str
        Path to the file containing the raw gene count matrix
    metadata_file_path : str
        Path to the file containing the metadata
    microenvironment_file_path : str
        Path to the file containing the microenvironment information this one is
        specific for the successful running of 'CellPhoneDB'

    simulate: bool
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
        # 'methods_to_compare': 'methods_to_compare_path'

    }

    def __init__(self, simulate: bool = False, **file_paths):
        self.input_files = {}
        self.simulate = simulate
        self._load_files(**file_paths)

        if 'seurat_object' not in self.input_files:
            if not set(['counts', 'metadata']).issubset(set(self.input_files.keys())):
                raise TypeError("Required at least the count matrix and metadata")

            print("Generating Seurat object for R packages comparison")

            adata = self.input_files["counts"].T
            meta = self.input_files["metadata"].set_index("Cell")
            adata.obs = meta

            if not os.path.exists("tmp"):
                os.makedirs("tmp")

            adata.write_h5ad(f"tmp/converted_seurat_object.h5ad")
            self.input_files["seurat_object"] = "tmp/converted_seurat_object.h5ad"
            
        else:
            pass

    def _load_files(self, **file_paths):
        ALLOWED_EXTENSIONS = [
            '.txt',
            '.csv',
            '.h5ad',
            '.tsv',
            '.yaml',
            '.yml'
        ]

        for key, path_attr in self.FILE_KEYS.items():
            file_path = file_paths.get(path_attr)
            if file_path:
                try:
                    if os.path.isfile(file_path):
                        _, file_extension = os.path.splitext(file_path)

                        if file_extension.lower() in ALLOWED_EXTENSIONS:
                            read_file_and_append(file_path, main_file_list=self.input_files, key=key)
                        else:
                            raise ValueError(f"File extension '{file_extension}' cannot be used as input")
                    else:
                        raise FileNotFoundError(f"File '{file_path}' not found.")
                except Exception as e:
                    print(f"Error loading file '{file_path}': {e}\n")


    def _test_successfully_loaded_file_message(self):
        print("Loaded the files successfully!\n")
