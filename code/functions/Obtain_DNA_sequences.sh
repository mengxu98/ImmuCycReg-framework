#!/bin/bash
arg1=$1 # Input file
arg2=$2 # Result file
arg3=$3 # hg38 file

# Define variables
BED_FILE="$arg1"
RESULTS_FILE="$arg2"
HG38_FILE="$arg3"
# http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/
HG38_URL="https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz"
TF_URL="http://alggen.lsi.upc.es/cgi-bin/promo_v3/promo/promoinit.cgi?dirDB=TF_8.3"
GENECARDS_URL="https://www.genecards.org/"

# Exit immediately if any command fails
set -e

# Define functions
install_bedtools() {
    if command -v bedtools >/dev/null 2>&1 ; then
        echo "bedtools is already installed..."
    else
        echo "Installing bedtools..."
        sudo apt-get update -q
        sudo apt-get install bedtools -y -q
    fi
}

download_hg38() {
    if [ -e "$HG38_FILE" ]; then
        echo "$HG38_FILE already exists..."
    else
        if [ ! -e "$HG38_FILE.gz" ]; then
            echo "Downloading $HG38_FILE..."
            wget "$HG38_URL"
        fi
        echo "Running gunzip..."
        gunzip "$HG38_FILE.gz"
    fi
}

get_sequences() {
    if [ -e "$HG38_FILE" ] && [ -e "$BED_FILE" ]; then
        echo "Extracting sequences from $HG38_FILE using $BED_FILE..."
        bedtools getfasta -fi "$HG38_FILE" -bed "$BED_FILE" > "$RESULTS_FILE"
    else
        echo "BED files do not exist..."
        exit 1
    fi
}

get_tfs() {
    echo "Getting transcription factor list from $TF_URL..."
}

get_gene_name() {
    echo "Determining gene name using $GENECARDS_URL..."
}

# Main script
echo "Running..."
install_bedtools
download_hg38
get_sequences
get_tfs
get_gene_name
echo "End..."
