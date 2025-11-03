import os
import sys
import subprocess
import urllib.request
import gzip
import shutil


def download_hg38(hg38_file):
    hg38_url = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz"
    hg38_gz = hg38_file + ".gz"
    
    if os.path.exists(hg38_file):
        print(f"{hg38_file} already exists...")
        return
    
    if os.path.exists(hg38_gz):
        print(f"Unzipping {hg38_gz}...")
        with gzip.open(hg38_gz, 'rb') as f_in:
            with open(hg38_file, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        return
    
    print(f"Downloading hg38 from {hg38_url}...")
    urllib.request.urlretrieve(hg38_url, hg38_gz)
    print(f"Unzipping {hg38_gz}...")
    with gzip.open(hg38_gz, 'rb') as f_in:
        with open(hg38_file, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    print(f"Downloaded {hg38_file}")


def check_bedtools():
    try:
        subprocess.run(['bedtools', '--version'], capture_output=True, check=True)
        return True
    except (subprocess.CalledProcessError, FileNotFoundError):
        return False


def get_sequences_bedtools(bed_file, hg38_file, results_file):
    if not check_bedtools():
        print("Error: bedtools is not installed. Please install bedtools first.")
        print("On macOS: brew install bedtools")
        print("On Linux: sudo apt-get install bedtools")
        return False
    
    print(f"Extracting sequences from {hg38_file} using {bed_file}...")
    with open(results_file, 'w') as f:
        subprocess.run(['bedtools', 'getfasta', '-fi', hg38_file, '-bed', bed_file],
                      stdout=f, check=True)
    return True


def get_sequences_pyfaidx(bed_file, hg38_file, results_file):
    try:
        from pyfaidx import Fasta
    except ImportError:
        print("Error: pyfaidx is not installed. Install with: pip install pyfaidx")
        return False
    
    print(f"Extracting sequences from {hg38_file} using {bed_file}...")
    fa = Fasta(hg38_file)
    
    with open(bed_file, 'r') as bed_f, open(results_file, 'w') as out_f:
        for line in bed_f:
            parts = line.strip().split('\t')
            if len(parts) >= 3:
                chrom = parts[0]
                start = int(parts[1])
                end = int(parts[2])
                
                if chrom in fa:
                    seq = fa[chrom][start:end].seq
                    name = parts[3] if len(parts) > 3 else f"{chrom}:{start}-{end}"
                    out_f.write(f">{name}\n{seq}\n")
    
    return True


def main():
    path_read = os.path.abspath(os.path.join(os.path.dirname(__file__), "data/"))
    path_save = os.path.abspath(os.path.join(os.path.dirname(__file__), "results/"))
    
    bed_file = os.path.join(path_save, "all_BED.bed")
    results_file = os.path.join(path_save, "DNA_sequences.txt")
    hg38_file = os.path.join(path_read, "hg38.fa")
    
    if not os.path.exists(bed_file):
        print(f"Error: BED file not found: {bed_file}")
        print("Please run Step-06_Peak.py first to generate the BED file.")
        return 1
    
    os.makedirs(os.path.dirname(hg38_file), exist_ok=True)
    download_hg38(hg38_file)
    
    if check_bedtools():
        success = get_sequences_bedtools(bed_file, hg38_file, results_file)
    else:
        print("bedtools not found, trying pyfaidx...")
        success = get_sequences_pyfaidx(bed_file, hg38_file, results_file)
    
    if success:
        print(f"Sequences extracted to {results_file}")
        return 0
    else:
        print("Failed to extract sequences")
        return 1


if __name__ == "__main__":
    sys.exit(main())

