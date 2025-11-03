import os
import sys
import pandas as pd
import numpy as np
import urllib.request
import gzip
import shutil


def fix_colnames(data):
    data.columns = [col[:15] for col in data.columns]
    return data


def save_csv(data, fileName, pathWay):
    os.makedirs(pathWay, exist_ok=True)
    filePath = os.path.join(
        pathWay, fileName.replace(".Rdata", ".csv").replace(".rdata", ".csv")
    )
    data.to_csv(filePath, index=True)
    return filePath


def fpkm_to_tpm(fpkm_matrix):
    fpkm = fpkm_matrix.values.astype(float)
    fpkm = np.maximum(fpkm, 1e-10)
    column_sums = np.sum(fpkm, axis=0, keepdims=True)
    tpm = fpkm / column_sums * 1e6
    tpm_df = pd.DataFrame(tpm, index=fpkm_matrix.index, columns=fpkm_matrix.columns)
    return tpm_df


def download_file(url, filepath):
    if os.path.exists(filepath):
        print(f"File already exists: {filepath}")
        return
    print(f"Downloading {url} to {filepath}...")
    os.makedirs(
        os.path.dirname(filepath) if os.path.dirname(filepath) else ".", exist_ok=True
    )
    urllib.request.urlretrieve(url, filepath)
    print(f"Downloaded: {filepath}")


def process_gtf(gtf_path):
    import gffutils

    db_path = gtf_path.replace(".gtf.gz", ".db")

    if not os.path.exists(db_path):
        print(f"Creating GTF database...")
        db = gffutils.create_db(
            gtf_path,
            dbfn=db_path,
            force=True,
            keep_order=True,
            merge_strategy="merge",
            sort_attribute_values=True,
        )
    else:
        db = gffutils.FeatureDB(db_path)

    geneinfo_list = []
    for feature in db.features_of_type("transcript"):
        geneinfo_list.append(
            {
                "seqnames": feature.chrom,
                "start": feature.start,
                "end": feature.stop,
                "width": feature.stop - feature.start + 1,
                "strand": feature.strand,
                "type": feature.featuretype,
                "gene_id": feature.attributes.get("gene_id", [""])[0],
                "transcript_id": feature.id,
                "gene_name": feature.attributes.get(
                    "gene_name", feature.attributes.get("Name", [""])[0]
                ),
            }
        )

    geneinfo_df = pd.DataFrame(geneinfo_list)
    return geneinfo_df


def main():
    pathRead = os.path.abspath(os.path.join(os.path.dirname(__file__), "data/"))
    pathSave = os.path.abspath(os.path.join(os.path.dirname(__file__), "results/"))
    os.makedirs(pathSave, exist_ok=True)

    dataClass = [
        "count-tcga-t",
        "fpkm-tcga-t",
        "fpkm-tcga-t_normlized",
        "fpkm-tcga-t_unnormlized",
    ]

    tcga_luad_final = None

    for i in range(len(dataClass)):
        print(f"Processing TCGA data: {dataClass[i]}...")
        filePath = os.path.join(pathRead, f"luad-rsem-{dataClass[i]}.txt.gz")

        tcga_luad = pd.read_csv(
            filePath, sep="\t", header=0, index_col=0, low_memory=False
        )
        tcga_luad = tcga_luad.iloc[:, 1:]
        tcga_luad = fix_colnames(tcga_luad)

        if "fpkm" in dataClass[i]:
            save_csv(tcga_luad, f"TCGA-LUAD-{dataClass[i]}.Rdata", pathSave)

            tcga_luad_tpm = fpkm_to_tpm(tcga_luad)
            save_csv(
                tcga_luad_tpm,
                f"TCGA-LUAD-{dataClass[i].replace('fpkm', 'tpm')}.Rdata",
                pathSave,
            )
        else:
            save_csv(tcga_luad, "TCGA-LUAD.Rdata", pathSave)
            tcga_luad_final = tcga_luad

    dataClass = ["count-gtex", "fpkm-gtex_normlized", "fpkm-gtex_unnormlized"]

    for i in range(len(dataClass)):
        print(f"Processing GTEx data: {dataClass[i]}...")
        filePath = os.path.join(pathRead, f"lung-rsem-{dataClass[i]}.txt.gz")

        gtex_luad = pd.read_csv(
            filePath, sep="\t", header=0, index_col=0, low_memory=False
        )
        gtex_luad = gtex_luad.iloc[:, 1:]
        gtex_luad = fix_colnames(gtex_luad)

        if "fpkm" in dataClass[i]:
            save_csv(gtex_luad, f"GTEx-LUAD-{dataClass[i]}.Rdata", pathSave)

            gtex_luad_tpm = fpkm_to_tpm(gtex_luad)
            save_csv(
                gtex_luad_tpm,
                f"GTEx-LUAD-{dataClass[i].replace('fpkm', 'tpm')}.Rdata",
                pathSave,
            )
        else:
            save_csv(gtex_luad, "GTEx-LUAD.Rdata", pathSave)

    print("Processing CNV data...")
    cnv_file = os.path.join(pathSave, "all_thresholded.by_genes_whitelisted.tsv")
    if not os.path.exists(cnv_file):
        download_file(
            "https://api.gdc.cancer.gov/data/7d64377f-2cea-4ee3-917f-8fcfbcd999e7",
            cnv_file,
        )

    tcga_cnv = pd.read_csv(cnv_file, sep="\t", header=0, index_col=0, low_memory=False)
    tcga_cnv = fix_colnames(tcga_cnv)

    if tcga_luad_final is not None:
        common_genes = list(set(tcga_luad_final.index) & set(tcga_cnv.index))
        common_samples = list(set(tcga_luad_final.columns) & set(tcga_cnv.columns))
        tcga_cnv = tcga_cnv.loc[common_genes, common_samples]

    save_csv(tcga_cnv, "CNV-LUAD.Rdata", pathSave)

    print("Processing Peak data...")
    peak_file = os.path.join(pathSave, "TCGA_ATAC_peak_Log2Counts_dedup_sample.gz")
    if not os.path.exists(peak_file):
        download_file(
            "https://tcgaatacseq.s3.us-east-1.amazonaws.com/download/TCGA_ATAC_peak_Log2Counts_dedup_sample.gz",
            peak_file,
        )

    peak_all = pd.read_csv(peak_file, sep="\t", header=0, index_col=0, low_memory=False)
    peak_all = fix_colnames(peak_all)

    luad_rows = [idx for idx in peak_all.index if "LUAD" in str(idx)]
    common_samples_peak = (
        list(set(tcga_luad_final.columns) & set(peak_all.columns))
        if tcga_luad_final is not None
        else peak_all.columns
    )
    peak_luad = peak_all.loc[luad_rows, common_samples_peak]

    if tcga_luad_final is not None:
        peak_samples = list(set(tcga_luad_final.columns) & set(peak_luad.columns))
        tcga_luad_peak_samples = tcga_luad_final.loc[:, peak_samples]
        save_csv(tcga_luad_peak_samples, "Peak-LUAD.Rdata", pathSave)

    print("Processing ATAC probe map...")
    probe_file = os.path.join(pathSave, "TCGA_ATAC_peak.all.probeMap")
    if not os.path.exists(probe_file):
        download_file(
            "https://tcgaatacseq.s3.us-east-1.amazonaws.com/download/TCGA_ATAC_peak.all.probeMap",
            probe_file,
        )

    genes_adj_peak = pd.read_csv(
        probe_file, sep="\t", header=0, index_col=0, low_memory=False
    )
    luad_peak_rows = [idx for idx in genes_adj_peak.index if "LUAD" in str(idx)]
    genes_adj_peak = genes_adj_peak.loc[luad_peak_rows, :]
    save_csv(genes_adj_peak, "ATAC-LUAD.Rdata", pathSave)

    print("Processing GTF file...")
    gtf_file = os.path.join(pathSave, "hg38.ncbiRefSeq.gtf.gz")
    if not os.path.exists(gtf_file):
        download_file(
            "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.ncbiRefSeq.gtf.gz",
            gtf_file,
        )

    try:
        geneinfo_df = process_gtf(gtf_file)
    except ImportError:
        print("Warning: gffutils not available, using simple GTF parser...")
        geneinfo_df = pd.read_csv(
            gtf_file,
            sep="\t",
            comment="#",
            header=None,
            low_memory=False,
            names=[
                "seqnames",
                "source",
                "type",
                "start",
                "end",
                "score",
                "strand",
                "phase",
                "attributes",
            ],
        )
        geneinfo_df = geneinfo_df[geneinfo_df["type"] == "transcript"].copy()
        geneinfo_df["width"] = geneinfo_df["end"] - geneinfo_df["start"] + 1

        def parse_attributes(attr_str):
            attrs = {}
            for item in attr_str.split(";"):
                item = item.strip()
                if "=" in item:
                    key, value = item.split("=", 1)
                    attrs[key.strip()] = value.strip().strip('"')
            return attrs

        geneinfo_df["attrs_dict"] = geneinfo_df["attributes"].apply(parse_attributes)
        geneinfo_df["gene_id"] = geneinfo_df["attrs_dict"].apply(
            lambda x: x.get("gene_id", "")
        )
        geneinfo_df["transcript_id"] = geneinfo_df["attrs_dict"].apply(
            lambda x: x.get("transcript_id", "")
        )
        geneinfo_df["gene_name"] = geneinfo_df["attrs_dict"].apply(
            lambda x: x.get("gene_name", x.get("Name", ""))
        )

        geneinfo_df = geneinfo_df[
            [
                "seqnames",
                "start",
                "end",
                "width",
                "strand",
                "type",
                "gene_id",
                "transcript_id",
                "gene_name",
            ]
        ]

    save_csv(geneinfo_df, "Geneinfo_df.Rdata", pathSave)

    print("All data processing completed!")
    return 0


if __name__ == "__main__":
    sys.exit(main())
