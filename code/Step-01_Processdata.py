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


def download_file(url, filepath, timeout=600):
    """Download file with timeout and error handling"""
    if os.path.exists(filepath):
        print(f"File already exists: {filepath}")
        return
    
    print(f"Downloading {url} to {filepath}...")
    print(f"Note: This may take a while for large files (timeout: {timeout}s)")
    
    os.makedirs(
        os.path.dirname(filepath) if os.path.dirname(filepath) else ".", exist_ok=True
    )
    
    try:
        # Use urlretrieve with timeout
        import socket
        socket.setdefaulttimeout(timeout)
        urllib.request.urlretrieve(url, filepath)
        print(f"Downloaded: {filepath}")
    except Exception as e:
        print(f"Error downloading {url}: {e}")
        print(f"You may need to download manually or increase timeout.")
        raise


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
    # Get the parent directory (project root)
    script_dir = os.path.dirname(os.path.abspath(__file__))
    base_dir = os.path.dirname(script_dir)  # Project root directory
    pathRead = os.path.join(base_dir, "data")
    path_save = os.path.join(base_dir, "results")
    os.makedirs(path_save, exist_ok=True)

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
            save_csv(tcga_luad, f"TCGA-LUAD-{dataClass[i]}.Rdata", path_save)

            tcga_luad_tpm = fpkm_to_tpm(tcga_luad)
            save_csv(
                tcga_luad_tpm,
                f"TCGA-LUAD-{dataClass[i].replace('fpkm', 'tpm')}.Rdata",
                path_save,
            )
        else:
            save_csv(tcga_luad, "TCGA-LUAD.Rdata", path_save)
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
            save_csv(gtex_luad, f"GTEx-LUAD-{dataClass[i]}.Rdata", path_save)

            gtex_luad_tpm = fpkm_to_tpm(gtex_luad)
            save_csv(
                gtex_luad_tpm,
                f"GTEx-LUAD-{dataClass[i].replace('fpkm', 'tpm')}.Rdata",
                path_save,
            )
        else:
            save_csv(gtex_luad, "GTEx-LUAD.Rdata", path_save)

    print("Processing CNV data...")
    cnv_file = os.path.join(path_save, "all_thresholded.by_genes_whitelisted.tsv")
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

    save_csv(tcga_cnv, "CNV-LUAD.Rdata", path_save)

    print("Processing Peak data...")
    peak_file = os.path.join(path_save, "TCGA_ATAC_peak_Log2Counts_dedup_sample.gz")
    
    # Check if file exists and is valid (not corrupted)
    import gzip
    file_valid = False
    expected_size_gb = 1.9  # Expected file size in GB
    
    if os.path.exists(peak_file):
        file_size_gb = os.path.getsize(peak_file) / (1024**3)
        print(f"Peak file exists: {peak_file} ({file_size_gb:.2f} GB)")
        
        # Check if file size is reasonable (at least 80% of expected)
        if file_size_gb < expected_size_gb * 0.8:
            print(f"Warning: File size ({file_size_gb:.2f} GB) is much smaller than expected ({expected_size_gb} GB)")
            print(f"File may be incomplete. Removing and re-downloading...")
            os.remove(peak_file)
            file_valid = False
        else:
            try:
                # Try to open and read a small portion to check if file is valid
                with gzip.open(peak_file, 'rt') as f:
                    f.read(1000)  # Read first 1000 bytes
                file_valid = True
                print(f"Peak file appears valid: {peak_file}")
            except (EOFError, gzip.BadGzipFile, OSError) as e:
                print(f"Warning: Peak file appears corrupted: {e}")
                print(f"Removing corrupted file and re-downloading...")
                os.remove(peak_file)
                file_valid = False
    else:
        print(f"Peak file not found: {peak_file}")
    
    if not file_valid:
        print("Downloading Peak data (this may take a while for large file)...")
        download_file(
            "https://tcgaatacseq.s3.us-east-1.amazonaws.com/download/TCGA_ATAC_peak_Log2Counts_dedup_sample.gz",
            peak_file,
            timeout=1800  # 30 minutes timeout for large file
        )

    try:
        peak_all = pd.read_csv(peak_file, sep="\t", header=0, index_col=0, low_memory=False, compression='gzip')
    except (EOFError, gzip.BadGzipFile) as e:
        print(f"Error reading peak file: {e}")
        print("File may be corrupted. Please download manually or check network connection.")
        print(f"Expected file size: ~1.9 GB")
        print(f"Current file size: {os.path.getsize(peak_file) / (1024**3):.2f} GB")
        raise
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
        save_csv(tcga_luad_peak_samples, "Peak-LUAD.Rdata", path_save)

    print("Processing ATAC probe map...")
    probe_file = os.path.join(path_save, "TCGA_ATAC_peak.all.probeMap")
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
    save_csv(genes_adj_peak, "ATAC-LUAD.Rdata", path_save)

    print("Processing GTF file...")
    gtf_file = os.path.join(path_save, "hg38.ncbiRefSeq.gtf.gz")
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

    save_csv(geneinfo_df, "Geneinfo_df.Rdata", path_save)

    print("All data processing completed!")
    return 0


if __name__ == "__main__":
    sys.exit(main())
