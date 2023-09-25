import os
import pandas as pd
import numpy as np

pathRead = "../data/"
pathSave = "../../Results/"

dataClass = ["count-tcga-t", "fpkm-tcga-t",
             "fpkm-tcga-t_normlized", "fpkm-tcga-t_unnormlized"]

def fix_colnames(data):
    return [col[:15] for col in data.columns]


def save_file(data, fileName, pathWay):
    data.to_csv(os.path.join(pathWay, fileName), sep="\t", index=False)


def fpkm_to_tpm(fpkm):
    return np.exp(np.log(fpkm) - np.log(np.sum(fpkm)) + np.log(1e6))


for i in range(len(dataClass)):
    tcga_luad = pd.read_csv(os.path.join(
        pathRead, f"luad-rsem-{dataClass[i]}.txt.gz"), sep="\t", header=0, index_col=1)
    tcga_luad_rownames = tcga_luad.iloc[:, 0]
    tcga_luad = tcga_luad.iloc[:, 1:]

    if "fpkm" in dataClass[i]:
        tcga_luad_tpm = fpkm_to_tpm(tcga_luad)
        tcga_luad.insert(loc= 0 , column="Gene", value= tcga_luad_rownames)
        save_file(tcga_luad, f"TCGA-LUAD-{dataClass[i]}.txt", pathSave)
        tcga_luad_tpm.insert(loc= 0 , column="Gene", value= tcga_luad_rownames)
        save_file(
            tcga_luad, f"TCGA-LUAD-{dataClass[i].replace('fpkm', 'tpm')}.txt", pathSave)
    else:
        tcga_luad.insert(loc= 0 , column="Gene", value= tcga_luad_rownames)
        save_file(tcga_luad, "TCGA-LUAD.txt", pathSave)
