#!/bin/bash

# Run complete pipeline for R and Python scripts
# This script runs the data processing pipeline step by step

set -e  # Exit on error

BASE_DIR="$(cd "$(dirname "$0")/.." && pwd)"
CODE_DIR="$BASE_DIR/code"
DATA_DIR="$BASE_DIR/data"
RESULTS_DIR="$BASE_DIR/results"

echo "=========================================="
echo "Running ImmuCycReg Framework Pipeline"
echo "=========================================="
echo "Base directory: $BASE_DIR"
echo "Code directory: $CODE_DIR"
echo "Data directory: $DATA_DIR"
echo "Results directory: $RESULTS_DIR"
echo ""

# Create results directory if it doesn't exist
mkdir -p "$RESULTS_DIR"

# Function to run a script and check its status
run_script() {
    local script_name=$1
    local description=$2
    local script_path="$CODE_DIR/$script_name"
    
    echo "=========================================="
    echo "Running: $description"
    echo "Script: $script_name"
    echo "=========================================="
    
    if [ ! -f "$script_path" ]; then
        echo "⚠ SKIPPED: Script not found: $script_name"
        return 1
    fi
    
    # Check if it's R or Python script
    if [[ "$script_name" == *.R ]]; then
        echo "Executing R script..."
        if Rscript "$script_path"; then
            echo "✓ $description completed successfully"
            return 0
        else
            echo "✗ $description failed"
            return 1
        fi
    elif [[ "$script_name" == *.py ]]; then
        echo "Executing Python script..."
        if python "$script_path"; then
            echo "✓ $description completed successfully"
            return 0
        else
            echo "✗ $description failed"
            return 1
        fi
    else
        echo "⚠ Unknown script type: $script_name"
        return 1
    fi
}

# Step 1: Data Processing
echo ""
echo "=========================================="
echo "Step 1: Data Processing"
echo "=========================================="

# Try R first, then Python
if [ -f "$CODE_DIR/Step-01_Processdata.R" ]; then
    run_script "Step-01_Processdata.R" "Step 01: Data Processing (R)"
fi

if [ -f "$CODE_DIR/Step-01_Processdata.py" ]; then
    run_script "Step-01_Processdata.py" "Step 01: Data Processing (Python)"
fi

# Step 2: NMF Consensus Clustering
echo ""
echo "=========================================="
echo "Step 2: NMF Consensus Clustering"
echo "=========================================="

if [ -f "$CODE_DIR/Step-02_NMF-ConsensuClustering.R" ]; then
    run_script "Step-02_NMF-ConsensuClustering.R" "Step 02: NMF Consensus Clustering (R)"
fi

if [ -f "$CODE_DIR/Step-02_NMF-ConsensuClustering.py" ]; then
    run_script "Step-02_NMF-ConsensuClustering.py" "Step 02: NMF Consensus Clustering (Python)"
fi

# Step 3: CIBERSORT Data Preparation
echo ""
echo "=========================================="
echo "Step 3: CIBERSORT Data Preparation"
echo "=========================================="

if [ -f "$CODE_DIR/Step-03_Cibersort-data.R" ]; then
    run_script "Step-03_Cibersort-data.R" "Step 03: CIBERSORT Data Preparation (R)"
fi

if [ -f "$CODE_DIR/Step-03_Cibersort-data.py" ]; then
    run_script "Step-03_Cibersort-data.py" "Step 03: CIBERSORT Data Preparation (Python)"
fi

# Step 4: DESeq2 Analysis
echo ""
echo "=========================================="
echo "Step 4: DESeq2 Analysis"
echo "=========================================="

if [ -f "$CODE_DIR/Step-04_DESeq2.R" ]; then
    run_script "Step-04_DESeq2.R" "Step 04: DESeq2 Analysis (R)"
fi

if [ -f "$CODE_DIR/Step-04_DESeq2.py" ]; then
    run_script "Step-04_DESeq2.py" "Step 04: DESeq2 Analysis (Python)"
fi

# Step 5: CNV Correlation
echo ""
echo "=========================================="
echo "Step 5: CNV Correlation"
echo "=========================================="

if [ -f "$CODE_DIR/Step-05_CNV-Correlation.R" ]; then
    run_script "Step-05_CNV-Correlation.R" "Step 05: CNV Correlation (R)"
fi

if [ -f "$CODE_DIR/Step-05_CNV-Correlation.py" ]; then
    run_script "Step-05_CNV-Correlation.py" "Step 05: CNV Correlation (Python)"
fi

# Step 6: Peak Analysis
echo ""
echo "=========================================="
echo "Step 6: Peak Analysis"
echo "=========================================="

if [ -f "$CODE_DIR/Step-06_Peak.R" ]; then
    run_script "Step-06_Peak.R" "Step 06: Peak Analysis (R)"
fi

if [ -f "$CODE_DIR/Step-06_Peak.py" ]; then
    run_script "Step-06_Peak.py" "Step 06: Peak Analysis (Python)"
fi

# Step 7: DNA Sequence Extraction
echo ""
echo "=========================================="
echo "Step 7: DNA Sequence Extraction"
echo "=========================================="

if [ -f "$CODE_DIR/Step-07_Obtain_DNA_sequences.R" ]; then
    run_script "Step-07_Obtain_DNA_sequences.R" "Step 07: DNA Sequence Extraction (R)"
fi

if [ -f "$CODE_DIR/Step-07_Obtain_DNA_sequences.py" ]; then
    run_script "Step-07_Obtain_DNA_sequences.py" "Step 07: DNA Sequence Extraction (Python)"
fi

# Step 8: ImmuCycReg Framework
echo ""
echo "=========================================="
echo "Step 8: ImmuCycReg Framework"
echo "=========================================="

if [ -f "$CODE_DIR/Step-08_ImmuCycReg-framework.R" ]; then
    run_script "Step-08_ImmuCycReg-framework.R" "Step 08: ImmuCycReg Framework (R)"
fi

if [ -f "$CODE_DIR/Step-08_ImmuCycReg-framework.py" ]; then
    run_script "Step-08_ImmuCycReg-framework.py" "Step 08: ImmuCycReg Framework (Python)"
fi

# Summary
echo ""
echo "=========================================="
echo "Pipeline Execution Complete"
echo "=========================================="
echo ""
echo "Results directory: $RESULTS_DIR"
echo "Generated files:"
ls -lh "$RESULTS_DIR" | head -20

echo ""
echo "=========================================="
echo "Next Steps:"
echo "=========================================="
echo "1. Check results in: $RESULTS_DIR"
echo "2. Compare R and Python outputs:"
echo "   python code/verify_results.py"
echo "3. Convert RData to CSV if needed:"
echo "   Rscript code/convert_rdata_to_csv.R"
echo ""

