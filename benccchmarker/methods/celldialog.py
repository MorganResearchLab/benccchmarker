import os
from pathlib import Path
import sys
import time
import tracemalloc
import yaml

import pandas as pd
import numpy as np
import anndata

import scanpy as sc

def main():
    args = sys.argv

    tracemalloc.start()
    start_time = time.time()

    current_file = Path(__file__)

    adata = anndata.read_h5ad(args[1])
    sc.pp.normalize_total(adata)

    celldialog_db_file_path = current_file.parent.parent / "datasets" / "defaults" / "celldialog" / "database" / "LRI.csv"

    main_output_path = "/".join(args[4].split("/")[:-1])
    output_path = f"{main_output_path}/tmp/celldialog"

    os.makedirs(output_path, exist_ok=True)

    cell_types = adata.obs['cell_type'].reset_index(drop=True)
    unique_cell_types = sorted(cell_types.unique())
    cell_type_dict = {}
    for i, cell_type in enumerate(unique_cell_types):
        cell_type_dict[cell_type] = i

    gene_data = adata.X.T
    gene = np.array(adata.var_names).reshape(-1, 1)
    cell_type_indices = {}

    for cell_type, key in cell_type_dict.items():
        cell_type_indices[cell_type] = cell_types[cell_types == cell_type].index

    row_sums = np.sum(gene_data, axis=1)/gene_data.shape[1]
    gene_data_float = gene_data.astype(float)
    std = np.std(gene_data_float, axis=1)
    total = row_sums + std

    cell_type_means = {}
    cell_type_arrs = {}

    for cell_type, key in cell_type_dict.items():
        cell_type_means[cell_type] = np.mean(gene_data[:, cell_type_indices[cell_type]], axis=1)
        cell_type_arrs[cell_type] = [1 if x > y else 0 for x, y in zip(cell_type_means[cell_type], total)]

    value_thre = np.append(gene, np.expand_dims(cell_type_arrs[unique_cell_types[0]], axis=0).T, axis=1)

    for unique_cell_type in unique_cell_types[1:]:
        value_thre = np.append(value_thre, np.expand_dims(cell_type_arrs[unique_cell_type], axis=0).T, axis=1)

    value_thre = pd.DataFrame(value_thre)

    cell_expression_thresholding_data_path = f'{output_path}/The_expression_thresholding_data.csv'
    value_thre.to_csv(cell_expression_thresholding_data_path, index=False, header=False)

    value_pro = np.append(gene, np.expand_dims(cell_type_means[unique_cell_types[0]], axis=0).T, axis=1)

    for unique_cell_type in unique_cell_types[1:]:
        value_pro = np.append(value_pro, np.expand_dims(cell_type_means[unique_cell_type], axis=0).T, axis=1)

    value_pro = pd.DataFrame(value_pro)
    cell_expression_product_data_path = f'{output_path}/The_expression_product_data.csv'
    value_pro.to_csv(cell_expression_product_data_path, index=False, header=False)

    cell_type_than_zero = {}
    cell_type_sum = {}

    for cell_type, key in cell_type_dict.items():
        cell_type_than_zero[cell_type] = (gene_data[:, cell_type_indices[cell_type]] > 0)
        cell_type_sum[cell_type] = np.sum(cell_type_than_zero[cell_type], axis=1)/gene_data[:, cell_type_indices[cell_type]].shape[1]

    value_cell = np.append(gene, np.expand_dims(cell_type_sum[unique_cell_types[0]], axis=0).T, axis=1)

    for unique_cell_type in unique_cell_types[1:]:
        value_cell = np.append(value_cell, np.expand_dims(cell_type_sum[unique_cell_type], axis=0).T, axis=1)

    value_cell = pd.DataFrame(value_cell)
    cell_expression_data_path = f'{output_path}/The_cell_expression_data.csv'
    value_cell.to_csv(cell_expression_data_path, index=False, header=False)

    end_preprocessing_time = time.time()

    folders_to_create = [
        f"{output_path}/inference/expression_thresholding",
        f"{output_path}/inference/expression_product",
        f"{output_path}/inference/cell_expression",
        f"{output_path}/inference/Three",
        f"{output_path}/inference/Three/TOP",
        f"{output_path}/inference"
    ]

    for folder in folders_to_create:
        os.makedirs(folder, exist_ok=True)

    inference = f"{output_path}/inference"
    cell_type = len(unique_cell_types) - 1 
    thr = pd.read_csv(cell_expression_thresholding_data_path, header=None, index_col=None).to_numpy()
    pro = pd.read_csv(cell_expression_product_data_path, header=None, index_col=None).to_numpy()
    cell = pd.read_csv(cell_expression_data_path, header=None, index_col=None).to_numpy()
    LRI = pd.read_csv(celldialog_db_file_path, header=None, index_col=None).to_numpy()
    
    for w in range(0, LRI.shape[0]):
        zhibiao1 = 0
        zhibiao2 = 0
        flag = False
        for x in range(0, thr.shape[0]):
            if LRI[w][0] == thr[x][0]:
                b11 = thr[x]
                zhibiao1 = 1
            if LRI[w][1] == thr[x][0]:
                b12 = thr[x]
                zhibiao2 = 1
            if zhibiao1 == 1 and zhibiao2 == 1:
                flag = True
                break
        if not flag:
            continue
        b11_l = b11[0]
        b11_r = b12[0]
        LRI_com = "{}_{}".format(b11_l, b11_r)
        b11 = np.delete(b11, 0)
        b12 = np.delete(b12, 0)
        for i in range(0, cell_type + 1):
            for j in range(0, cell_type + 1):
                value0 = 0
                value1 = 1
                if b11[i] == b12[j] and b11[i] == 1 and b12[j] == 1:
                    with open(inference + "/expression_thresholding/" + str(i) + str(j) + ".csv", mode="a") as f:
                        f.write("{},{}\n".format(LRI_com, value1))
                    f.close()
                else:
                    with open(inference + "/expression_thresholding/" + str(i) + str(j) + ".csv", mode="a") as f:
                        f.write("{},{}\n".format(LRI_com, value0))
                    f.close()

    for i in range(0, cell_type + 1):
        for j in range(0, cell_type + 1):
            a1 = pd.read_csv(inference + "/expression_thresholding/" + str(i) + str(j) + ".csv", header=None,
                            index_col=None).to_numpy()  # Obtain the result
            SUM = sum(a1[:, 1])
            with open(inference + "/thresholding_result.csv", mode="a") as f:
                f.write(str(i))
                f.write(str(j))
                f.write('____')
                f.write(str(SUM))
                f.write('\n')
            f.close()

    #  The expression product calculation method code

    for w in range(0, LRI.shape[0]):
        zhibiao1 = 0
        zhibiao2 = 0
        flag = False
        for x in range(0, pro.shape[0]):
            if LRI[w][0] == pro[x][0]:
                b11 = pro[x]
                zhibiao1 = 1
            if LRI[w][1] == pro[x][0]:
                b12 = pro[x]
                zhibiao2 = 1
            if zhibiao1 == 1 and zhibiao2 == 1:
                flag = True
                break
        if not flag:
            continue
        b11_l = b11[0]
        b11_r = b12[0]
        LRI_com = "{}_{}".format(b11_l, b11_r)
        b11 = np.delete(b11, 0)
        b12 = np.delete(b12, 0)
        for i in range(0, cell_type + 1):
            for j in range(0, cell_type + 1):
                Cheng = b11[i] * b12[j]
                with open(inference + "/expression_product/" + str(i) + str(j) + ".csv", mode="a") as f:
                    f.write("{},{}\n".format(LRI_com, Cheng))
                f.close()

    for i in range(0, cell_type + 1):
        for j in range(0, cell_type + 1):
            a1 = pd.read_csv(inference + "/expression_product/" + str(i) + str(j) + ".csv", header=None, index_col=None).to_numpy()  # Obtain the result
            SUM = sum(a1[:, 1])
            with open(inference + "/product_result.csv", mode="a") as f:
                f.write(str(i))
                f.write(str(j))
                f.write('____')
                f.write(str(SUM))
                f.write('\n')
            f.close()

    # The cell expression calculation method code

    for w in range(0, LRI.shape[0]):
        zhibiao1 = 0
        zhibiao2 = 0
        flag = False
        for x in range(0, cell.shape[0]):
            if LRI[w][0] == cell[x][0]:
                b11 = cell[x]
                zhibiao1 = 1
            if LRI[w][1] == cell[x][0]:
                b12 = cell[x]
                zhibiao2 = 1
            if zhibiao1 == 1 and zhibiao2 == 1:
                flag = True
                break
        if not flag:
            continue
        b11_l = b11[0]
        b11_r = b12[0]
        LRI_com = "{}_{}".format(b11_l, b11_r)
        b11 = np.delete(b11, 0)
        b12 = np.delete(b12, 0)
        for i in range(0, cell_type + 1):
            for j in range(0, cell_type + 1):
                cell_ = b11[i] * b12[j]
                with open(inference + "/cell_expression/" + str(i) + str(j) + ".csv", mode="a") as f:
                    f.write("{},{}\n".format(LRI_com, cell_))
                f.close()

    for i in range(0, cell_type + 1):
        for j in range(0, cell_type + 1):
            a1 = pd.read_csv(inference + "/cell_expression/" + str(i) + str(j) + ".csv", header=None, index_col=None).to_numpy()  # Obtain the result
            SUM = sum(a1[:, 1])
            with open(inference + "/cell_result.csv", mode="a") as f:
                f.write(str(i))
                f.write(str(j))
                f.write('____')
                f.write(str(SUM))
                f.write('\n')
            f.close()

    # Processing data
    sum_data = 0
    data = pd.read_csv(inference + "/thresholding_result.csv", header=None,index_col=None).to_numpy()
    data2 = []
    for j in range(len(data)):
        data11 = data[j][0]
        data1 = float(data11[6:])
        a1 = data11[:2]
        if a1[0] == a1[1]:
            data1 = 0.0
        data2.append(data1)
    # normalization
    _range = np.max(data2) - np.min(data2)
    aaa = (data2 - np.min(data2)) / _range
    sum_data = pd.DataFrame(aaa)
    totol = sum_data.values
    result = totol.reshape((cell_type + 1, cell_type + 1))
    result1 = pd.DataFrame(result)

    sum_data = 0
    data = pd.read_csv(inference + "/product_result.csv", header=None,index_col=None).to_numpy()
    data2 = []
    for j in range(len(data)):
        data11 = data[j][0]
        data1 = float(data11[6:])
        a1 = data11[:2]
        if a1[0] == a1[1]:
            data1 = 0.0
        data2.append(data1)
    # normalization
    _range = np.max(data2) - np.min(data2)
    aaa = (data2 - np.min(data2)) / _range
    sum_data = pd.DataFrame(aaa)
    totol = sum_data.values
    result = totol.reshape((cell_type + 1, cell_type + 1))
    result2 = pd.DataFrame(result)

    sum_data = 0
    data = pd.read_csv(inference + "/cell_result.csv", header=None,index_col=None).to_numpy()
    data2 = []
    for j in range(len(data)):
        data11 = data[j][0]
        data1 = float(data11[6:])
        a1 = data11[:2]
        if a1[0] == a1[1]:
            data1 = 0.0
        data2.append(data1)
    # normalization
    _range = np.max(data2) - np.min(data2)
    aaa = (data2 - np.min(data2)) / _range
    sum_data = pd.DataFrame(aaa)
    totol = sum_data.values
    result = totol.reshape((cell_type + 1, cell_type + 1))
    result3 = pd.DataFrame(result)

    # The three-point estimation method

    result_max = np.maximum(result1, result2)
    result_med = np.median([result1, result2, result3], axis=0)
    result_min = np.minimum(result1, result2)
    result_matrix = np.maximum(result_max, result3) + result_med * 4 + np.minimum(result_min, result3)
    result_matrix /= 6
    result_matrix = pd.DataFrame(result_matrix)

    sum_data = 0
    for i in range(0, cell_type + 1):
        for j in range(0, cell_type + 1):
            list1 = pd.read_csv(inference + "/expression_product/" + str(i) + str(j) + ".csv", header=None, index_col=None).to_numpy()  # 得分
            _range1 = np.max(list1[:, 1]) - np.min(list1[:, 1])
            aaa1 = (list1[:, 1] - np.min(list1[:, 1])) / _range1
            count = np.sum(aaa1 > 0.05)
            with open(inference + "/pro_LRi_num.csv",mode="a") as f:
                f.write(str(i))
                f.write(str(j))
                f.write('____')
                f.write(str(count))
                f.write('\n')
            f.close()
    for i in range(0, cell_type + 1):
        for j in range(0, cell_type + 1):
            list2 = pd.read_csv(inference + "/cell_expression/" + str(i) + str(j) + ".csv", header=None, index_col=None).to_numpy()  # 得分
            _range2 = np.max(list2[:, 1]) - np.min(list2[:, 1])
            try:
                # Mat have zero division error
                
                aaa2 = (list2[:, 1] - np.min(list2[:, 1])) / _range2
            except:
                aaa2 = 0
            count = np.sum(aaa2 > 0.05)
            with open(inference + "/cell_LRi_num.csv",mode="a") as f:
                f.write(str(i))
                f.write(str(j))
                f.write('____')
                f.write(str(count))
                f.write('\n')
            f.close()
    for i in range(0, cell_type + 1):
        for j in range(0, cell_type + 1):
            list1 = pd.read_csv(inference + "/expression_product/" + str(i) + str(j) + ".csv", header=None,
                                index_col=None).to_numpy()
            _range1 = np.max(list1[:, 1]) - np.min(list1[:, 1])
            aaa1 = (list1[:, 1] - np.min(list1[:, 1])) / _range1

            list2 = pd.read_csv(inference + "/cell_expression/" + str(i) + str(j) + ".csv", header=None,
                                index_col=None).to_numpy()
            _range2 = np.max(list2[:, 1]) - np.min(list2[:, 1])
            try:
                # Mat have zero division error
                aaa2 = (list2[:, 1] - np.min(list2[:, 1])) / _range2
            except:
                aaa2 = 0
            aaa = (aaa1 + aaa2) / 2
            gene = list1[:, 0]
            gene_val = np.column_stack((gene, aaa))
            df = pd.DataFrame(gene_val)
            output_path = inference + "/Three"
            output_name = '{}{}.csv'
            output_file = os.path.join(output_path, output_name.format(i, j))
            df.to_csv(output_file, index=False, header=False)
            count = np.sum(aaa > 0.05)
            with open(inference + "/Three_LRi_num.csv", mode="a") as f:
                f.write(str(i))
                f.write(str(j))
                f.write('____')
                f.write(str(count))
                f.write('\n')
            f.close()
    data = pd.read_csv(inference + "/Three_LRi_num.csv", header=None, index_col=None).to_numpy()
    data2 = []
    for j in range(len(data)):
        data11 = data[j][0]
        data1 = float(data11[6:])
        a1 = data11[:2]
        if a1[0] == a1[1]:
            data1 = 0.0
        data2.append(data1)
    sum_data = pd.DataFrame(data2)

    totol = sum_data.values
    result = totol.reshape((cell_type + 1, cell_type + 1))
    result = pd.DataFrame(result)

    end_inference_time = time.time()

    # Change key to value, value to key of cell_type_dict
    cell_type_dict = {v: k for k, v in cell_type_dict.items()}
    cell_type_dict
    path = inference + "/Three"

    files = os.listdir(path)
    files = [f for f in files if f.endswith('.csv')]

    result_dfs = []
    for file in files:
        result_df = pd.read_csv(os.path.join(path, file), header=None)
        source = file[0]
        target = file[1]
        result_df['source'] = source
        result_df['target'] = target
        result_df['source'] = result_df['source'].apply(lambda x: cell_type_dict[int(x)]) 
        result_df['target'] = result_df['target'].apply(lambda x: cell_type_dict[int(x)])

        result_dfs.append(result_df)

    result_df = pd.concat(result_dfs, axis=0)
    # Split column 0 into two columns ligand and receptor
    lrresult_df = result_df[0].str.split('_', expand=True).rename(columns={0: 'ligand', 1: 'receptor'})

    # Concatenate the new columns to the original dataframe
    result_df = pd.concat([result_df, lrresult_df], axis=1)
    result_df.drop(columns=0, inplace=True)

    result_df.rename(columns={1: 'value'}, inplace=True)
    result_df["label"] = result_df["source"] + "---" + result_df["ligand"] + "---" + result_df["target"] + "---" + result_df["receptor"]
    result_df = result_df.query("value > 0")

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