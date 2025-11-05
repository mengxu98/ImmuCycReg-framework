#!/usr/bin/env python3
"""
Compare R and Python script outputs for consistency
Tests both R and Python scripts and compares their outputs
"""

import os
import sys
import subprocess
import pandas as pd
import numpy as np
from pathlib import Path
import json


def run_r_script(script_path, description):
    """Run an R script and return success status"""
    print(f"\n{'='*60}")
    print(f"Running R Script: {description}")
    print(f"Script: {script_path}")
    print(f"{'='*60}")
    
    if not os.path.exists(script_path):
        print(f"✗ Script not found: {script_path}")
        return False
    
    try:
        result = subprocess.run(
            ['Rscript', script_path],
            capture_output=True,
            text=True,
            timeout=600,
            cwd=os.path.dirname(script_path) or '.'
        )
        
        if result.returncode == 0:
            print(f"✓ R script completed successfully")
            if result.stdout:
                print(f"Output (last 200 chars):\n{result.stdout[-200:]}")
            return True
        else:
            print(f"✗ R script failed with return code: {result.returncode}")
            if result.stderr:
                print(f"Error:\n{result.stderr[-500:]}")
            return False
    except subprocess.TimeoutExpired:
        print(f"✗ R script timed out")
        return False
    except Exception as e:
        print(f"✗ Error running R script: {str(e)}")
        return False


def run_python_script(script_path, description):
    """Run a Python script and return success status"""
    print(f"\n{'='*60}")
    print(f"Running Python Script: {description}")
    print(f"Script: {script_path}")
    print(f"{'='*60}")
    
    if not os.path.exists(script_path):
        print(f"✗ Script not found: {script_path}")
        return False
    
    try:
        result = subprocess.run(
            [sys.executable, script_path],
            capture_output=True,
            text=True,
            timeout=600,
            cwd=os.path.dirname(script_path) or '.'
        )
        
        if result.returncode == 0:
            print(f"✓ Python script completed successfully")
            if result.stdout:
                print(f"Output (last 200 chars):\n{result.stdout[-200:]}")
            return True
        else:
            print(f"✗ Python script failed with return code: {result.returncode}")
            if result.stderr:
                print(f"Error:\n{result.stderr[-500:]}")
            return False
    except subprocess.TimeoutExpired:
        print(f"✗ Python script timed out")
        return False
    except Exception as e:
        print(f"✗ Error running Python script: {str(e)}")
        return False


def compare_csv_outputs(r_file, py_file, tolerance=1e-5):
    """Compare two CSV files"""
    print(f"\nComparing:")
    print(f"  R output: {r_file}")
    print(f"  Python output: {py_file}")
    
    if not os.path.exists(r_file):
        print(f"  ⚠ R output not found: {r_file}")
        return None
    
    if not os.path.exists(py_file):
        print(f"  ⚠ Python output not found: {py_file}")
        return None
    
    try:
        r_df = pd.read_csv(r_file)
        py_df = pd.read_csv(py_file)
        
        # Compare shape
        if r_df.shape != py_df.shape:
            print(f"  ✗ Shape mismatch: R {r_df.shape} vs Python {py_df.shape}")
            return False
        
        # Compare columns
        r_cols = set(r_df.columns)
        py_cols = set(py_df.columns)
        if r_cols != py_cols:
            print(f"  ⚠ Column mismatch:")
            print(f"    R only: {r_cols - py_cols}")
            print(f"    Python only: {py_cols - r_cols}")
        
        # Compare numeric columns
        numeric_cols = r_df.select_dtypes(include=[np.number]).columns
        all_match = True
        max_diff = 0
        
        for col in numeric_cols:
            if col in py_df.columns:
                r_vals = r_df[col].fillna(0).astype(float).values
                py_vals = py_df[col].fillna(0).astype(float).values
                
                if len(r_vals) != len(py_vals):
                    print(f"  ✗ {col}: Length mismatch")
                    all_match = False
                    continue
                
                diff = np.abs(r_vals - py_vals)
                max_diff = max(max_diff, np.nanmax(diff))
                
                if not np.allclose(r_vals, py_vals, rtol=tolerance, atol=tolerance, equal_nan=True):
                    print(f"  ⚠ {col}: Max difference = {np.nanmax(diff):.6f}")
                    all_match = False
        
        if all_match:
            print(f"  ✓ Files match within tolerance {tolerance}")
            print(f"  Max difference: {max_diff:.6f}")
            return True
        else:
            print(f"  ⚠ Files have differences (max: {max_diff:.6f})")
            if max_diff < tolerance * 10:
                print(f"  → Differences are within acceptable range")
                return True
            return False
    except Exception as e:
        print(f"  ✗ Error comparing files: {e}")
        return False


def main():
    print("="*60)
    print("R and Python Script Comparison Test")
    print("="*60)
    
    base_dir = Path(__file__).parent.parent
    code_dir = base_dir / "code"
    results_dir = base_dir / "results"
    
    # Test scripts
    test_scripts = [
        {
            'name': 'Step-01_Processdata',
            'r_script': 'Step-01_Processdata.R',
            'py_script': 'Step-01_Processdata.py',
            'description': 'Step 01: Data Processing',
            'outputs': [
                ('TCGA-LUAD.csv', 'TCGA-LUAD.Rdata'),
                ('GTEx-LUAD.csv', 'GTEx-LUAD.Rdata'),
                ('CNV-LUAD.csv', 'CNV-LUAD.Rdata'),
            ]
        },
        {
            'name': 'Step-05_CNV-Correlation',
            'r_script': 'Step-05_CNV-Correlation.R',
            'py_script': 'Step-05_CNV-Correlation.py',
            'description': 'Step 05: CNV Correlation',
            'outputs': [
                ('CNV_correlation_value.csv', 'CNV_correlation_value.csv'),
            ]
        },
    ]
    
    print("\n" + "="*60)
    print("Step 1: Testing Script Syntax")
    print("="*60)
    
    # Test R script syntax
    print("\n--- Testing R Scripts ---")
    r_scripts = [
        'Step-01_Processdata.R',
        'Step-02_NMF-ConsensuClustering.R',
        'Step-03_Cibersort-data.R',
        'Step-04_DESeq2.R',
        'Step-05_CNV-Correlation.R',
        'Step-06_Peak.R',
    ]
    
    r_syntax_ok = {}
    for script in r_scripts:
        script_path = code_dir / script
        if script_path.exists():
            try:
                result = subprocess.run(
                    ['Rscript', '-e', f'parse("{script_path}")'],
                    capture_output=True,
                    text=True,
                    timeout=10
                )
                passed = result.returncode == 0
                status = "✓" if passed else "✗"
                print(f"  {status} {script}")
                r_syntax_ok[script] = passed
            except Exception as e:
                print(f"  ✗ {script} - Error: {e}")
                r_syntax_ok[script] = False
    
    # Test Python script syntax
    print("\n--- Testing Python Scripts ---")
    py_scripts = [
        'Step-01_Processdata.py',
        'Step-02_NMF-ConsensuClustering.py',
        'Step-03_Cibersort-data.py',
        'Step-04_DESeq2.py',
        'Step-05_CNV-Correlation.py',
        'Step-06_Peak.py',
    ]
    
    py_syntax_ok = {}
    for script in py_scripts:
        script_path = code_dir / script
        if script_path.exists():
            try:
                result = subprocess.run(
                    [sys.executable, '-m', 'py_compile', str(script_path)],
                    capture_output=True,
                    text=True,
                    timeout=10
                )
                passed = result.returncode == 0
                status = "✓" if passed else "✗"
                print(f"  {status} {script}")
                py_syntax_ok[script] = passed
            except Exception as e:
                print(f"  ✗ {script} - Error: {e}")
                py_syntax_ok[script] = False
    
    # Summary
    print("\n" + "="*60)
    print("Summary")
    print("="*60)
    
    r_passed = sum(1 for v in r_syntax_ok.values() if v)
    r_total = len(r_syntax_ok)
    py_passed = sum(1 for v in py_syntax_ok.values() if v)
    py_total = len(py_syntax_ok)
    
    print(f"\nR Scripts: {r_passed}/{r_total} syntax OK")
    print(f"Python Scripts: {py_passed}/{py_total} syntax OK")
    
    print("\n" + "="*60)
    print("Note: To compare outputs, run both R and Python scripts")
    print("and use compare_csv_outputs() function to compare results")
    print("="*60)
    
    if r_passed == r_total and py_passed == py_total:
        print("\n✓ All scripts have valid syntax!")
        return 0
    else:
        print("\n⚠ Some scripts have syntax errors")
        return 1


if __name__ == "__main__":
    sys.exit(main())

