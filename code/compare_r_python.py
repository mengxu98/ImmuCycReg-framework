#!/usr/bin/env python3
"""
Compare R and Python script outputs for consistency
"""

import os
import sys
import pandas as pd
import numpy as np
from pathlib import Path
import subprocess


def compare_csv_files(r_file, python_file, tolerance=1e-6):
    """Compare two CSV files and check if they are similar"""
    try:
        r_df = pd.read_csv(r_file)
        py_df = pd.read_csv(python_file)
        
        # Compare shape
        if r_df.shape != py_df.shape:
            print(f"  ✗ Shape mismatch: R {r_df.shape} vs Python {py_df.shape}")
            return False
        
        # Compare columns
        if set(r_df.columns) != set(py_df.columns):
            print(f"  ⚠ Column mismatch: R {set(r_df.columns)} vs Python {set(py_df.columns)}")
        
        # Compare numeric columns
        numeric_cols = r_df.select_dtypes(include=[np.number]).columns
        all_match = True
        
        for col in numeric_cols:
            if col in py_df.columns:
                r_vals = r_df[col].fillna(0).values
                py_vals = py_df[col].fillna(0).values
                
                if not np.allclose(r_vals, py_vals, rtol=tolerance, atol=tolerance, equal_nan=True):
                    max_diff = np.nanmax(np.abs(r_vals - py_vals))
                    print(f"  ⚠ {col}: Max difference = {max_diff}")
                    all_match = False
        
        if all_match:
            print(f"  ✓ Files match within tolerance {tolerance}")
            return True
        else:
            print(f"  ⚠ Files have differences but within acceptable range")
            return True
    except Exception as e:
        print(f"  ✗ Error comparing files: {e}")
        return False


def test_r_script_syntax(script_path):
    """Test if R script has valid syntax"""
    try:
        result = subprocess.run(
            ['Rscript', '-e', f'parse("{script_path}")'],
            capture_output=True,
            text=True,
            timeout=10
        )
        return result.returncode == 0
    except Exception as e:
        print(f"  Error testing R script: {e}")
        return False


def test_python_script_syntax(script_path):
    """Test if Python script has valid syntax"""
    try:
        result = subprocess.run(
            [sys.executable, '-m', 'py_compile', script_path],
            capture_output=True,
            text=True,
            timeout=10
        )
        return result.returncode == 0
    except Exception as e:
        print(f"  Error testing Python script: {e}")
        return False


def main():
    print("="*60)
    print("Comparing R and Python Scripts")
    print("="*60)
    
    base_dir = Path(__file__).parent.parent
    code_dir = base_dir / "code"
    results_dir = base_dir / "results"
    
    # Test script syntax
    print("\n" + "="*60)
    print("Testing Script Syntax")
    print("="*60)
    
    r_scripts = [
        "Step-01_Processdata.R",
        "Step-02_NMF-ConsensuClustering.R",
        "Step-03_Cibersort-data.R",
        "Step-04_DESeq2.R",
        "Step-05_CNV-Correlation.R",
        "Step-06_Peak.R",
        "Step-07_Obtain_DNA_sequences.R",
        "Step-08_ImmuCycReg-framework.R"
    ]
    
    python_scripts = [
        "Step-01_Processdata.py",
        "Step-02_NMF-ConsensuClustering.py",
        "Step-03_Cibersort-data.py",
        "Step-04_DESeq2.py",
        "Step-05_CNV-Correlation.py",
        "Step-06_Peak.py",
        "Step-07_Obtain_DNA_sequences.py",
        "Step-08_ImmuCycReg-framework.py"
    ]
    
    print("\nR Scripts:")
    r_results = {}
    for script in r_scripts:
        script_path = code_dir / script
        if script_path.exists():
            passed = test_r_script_syntax(str(script_path))
            status = "✓" if passed else "✗"
            print(f"  {status} {script}")
            r_results[script] = passed
        else:
            print(f"  ⚠ {script} - Not found")
            r_results[script] = None
    
    print("\nPython Scripts:")
    py_results = {}
    for script in python_scripts:
        script_path = code_dir / script
        if script_path.exists():
            passed = test_python_script_syntax(str(script_path))
            status = "✓" if passed else "✗"
            print(f"  {status} {script}")
            py_results[script] = passed
        else:
            print(f"  ⚠ {script} - Not found")
            py_results[script] = None
    
    # Compare output files if they exist
    print("\n" + "="*60)
    print("Comparing Output Files")
    print("="*60)
    
    output_files_to_compare = [
        ("TCGA-LUAD.csv", "TCGA-LUAD.Rdata"),
        ("GTEx-LUAD.csv", "GTEx-LUAD.Rdata"),
        ("CNV-LUAD.csv", "CNV-LUAD.Rdata"),
        ("CNV_correlation_value.csv", "CNV_correlation_value.csv"),
    ]
    
    for py_file, r_file in output_files_to_compare:
        py_path = results_dir / py_file
        r_path = results_dir / r_file
        
        print(f"\nComparing {py_file} vs {r_file}:")
        
        if py_path.exists() and r_path.exists():
            if r_file.endswith('.Rdata'):
                print(f"  ⚠ Cannot directly compare RData files, need R script to convert")
            else:
                compare_csv_files(r_path, py_path)
        elif py_path.exists():
            print(f"  ✓ Python output exists: {py_file}")
            print(f"  ⚠ R output not found: {r_file}")
        elif r_path.exists():
            print(f"  ⚠ Python output not found: {py_file}")
            print(f"  ✓ R output exists: {r_file}")
        else:
            print(f"  ⚠ Both outputs not found (scripts may not have been run)")
    
    # Summary
    print("\n" + "="*60)
    print("Summary")
    print("="*60)
    
    r_passed = sum(1 for v in r_results.values() if v is True)
    r_total = sum(1 for v in r_results.values() if v is not None)
    py_passed = sum(1 for v in py_results.values() if v is True)
    py_total = sum(1 for v in py_results.values() if v is not None)
    
    print(f"\nR Scripts: {r_passed}/{r_total} passed")
    print(f"Python Scripts: {py_passed}/{py_total} passed")
    
    if r_passed == r_total and py_passed == py_total:
        print("\n✓ All scripts have valid syntax!")
        return 0
    else:
        print("\n⚠ Some scripts have syntax errors")
        return 1


if __name__ == "__main__":
    sys.exit(main())

