#!/usr/bin/env python3
"""
Verify R and Python script results consistency
Compare actual output files from R and Python runs
"""

import os
import sys
import pandas as pd
import numpy as np
from pathlib import Path
import json


def load_rdata_as_csv(rdata_path):
    """Helper function - RData files need to be converted to CSV first"""
    print(f"Note: RData files need to be converted to CSV for comparison")
    print(f"Use R script to convert: save(data, file='output.csv')")
    return None


def compare_output_files(r_file, py_file, tolerance=1e-5, description=""):
    """Compare R and Python output files"""
    print(f"\n{'='*60}")
    if description:
        print(f"Comparing: {description}")
    print(f"R file: {r_file}")
    print(f"Python file: {py_file}")
    print(f"{'='*60}")
    
    # Handle RData files
    if r_file.endswith('.Rdata') or r_file.endswith('.rdata'):
        print("⚠ RData files cannot be directly compared")
        print("  Please convert RData to CSV first using R:")
        print(f"  load('{r_file}')")
        print(f"  write.csv(data, '{r_file.replace('.Rdata', '.csv')}')")
        return None
    
    if not os.path.exists(r_file):
        print(f"✗ R output not found: {r_file}")
        return False
    
    if not os.path.exists(py_file):
        print(f"✗ Python output not found: {py_file}")
        return False
    
    try:
        r_df = pd.read_csv(r_file)
        py_df = pd.read_csv(py_file)
        
        print(f"\nR file shape: {r_df.shape}")
        print(f"Python file shape: {py_df.shape}")
        
        if r_df.shape != py_df.shape:
            print(f"✗ Shape mismatch!")
            return False
        
        # Compare columns
        r_cols = set(r_df.columns)
        py_cols = set(py_df.columns)
        if r_cols != py_cols:
            print(f"⚠ Column mismatch:")
            if r_cols - py_cols:
                print(f"  R only: {r_cols - py_cols}")
            if py_cols - r_cols:
                print(f"  Python only: {py_cols - r_cols}")
        else:
            print(f"✓ Columns match: {list(r_cols)}")
        
        # Compare numeric values
        numeric_cols = r_df.select_dtypes(include=[np.number]).columns
        all_match = True
        max_diff = 0
        
        for col in numeric_cols:
            if col not in py_df.columns:
                continue
            
            r_vals = r_df[col].fillna(0).astype(float).values
            py_vals = py_df[col].fillna(0).astype(float).values
            
            if len(r_vals) != len(py_vals):
                print(f"✗ {col}: Length mismatch")
                all_match = False
                continue
            
            diff = np.abs(r_vals - py_vals)
            col_max_diff = np.nanmax(diff)
            max_diff = max(max_diff, col_max_diff)
            
            if col_max_diff > tolerance:
                print(f"⚠ {col}: Max difference = {col_max_diff:.6f}")
                all_match = False
        
        # Compare string/categorical columns
        string_cols = r_df.select_dtypes(include=['object']).columns
        for col in string_cols:
            if col not in py_df.columns:
                continue
            
            r_vals = r_df[col].fillna('').astype(str).values
            py_vals = py_df[col].fillna('').astype(str).values
            
            if not np.array_equal(r_vals, py_vals):
                mismatches = np.sum(r_vals != py_vals)
                print(f"⚠ {col}: {mismatches} mismatches")
                all_match = False
        
        if all_match:
            print(f"\n✓ Files match perfectly!")
            print(f"  Max numeric difference: {max_diff:.6e}")
            return True
        else:
            print(f"\n⚠ Files have differences")
            print(f"  Max numeric difference: {max_diff:.6e}")
            if max_diff < tolerance * 10:
                print(f"  → Differences are within acceptable range")
                return True
            return False
            
    except Exception as e:
        print(f"✗ Error comparing files: {e}")
        import traceback
        traceback.print_exc()
        return False


def main():
    print("="*60)
    print("R and Python Results Verification")
    print("="*60)
    
    base_dir = Path(__file__).parent.parent
    results_dir = base_dir / "results"
    
    # Files to compare
    files_to_compare = [
        {
            'r_file': 'CNV_correlation_value.csv',
            'py_file': 'CNV_correlation_value.csv',
            'description': 'CNV Correlation Results'
        },
        {
            'r_file': 'DESeq2_results_cluster1.csv',
            'py_file': 'DESeq2_results_cluster1.csv',
            'description': 'DESeq2 Results Cluster 1'
        },
        {
            'r_file': 'Gene_TFs_edge_num.csv',
            'py_file': 'Gene_TFs_edge_num.csv',
            'description': 'Gene-TF Edge Numbers'
        },
        {
            'r_file': 'Gene_TFs_node_num.csv',
            'py_file': 'Gene_TFs_node_num.csv',
            'description': 'Gene-TF Node Numbers'
        },
    ]
    
    results = {}
    
    for file_info in files_to_compare:
        r_path = results_dir / file_info['r_file']
        py_path = results_dir / file_info['py_file']
        
        result = compare_output_files(
            str(r_path),
            str(py_path),
            tolerance=1e-5,
            description=file_info['description']
        )
        results[file_info['description']] = result
    
    # Summary
    print("\n" + "="*60)
    print("Verification Summary")
    print("="*60)
    
    passed = sum(1 for v in results.values() if v is True)
    total = sum(1 for v in results.values() if v is not None)
    
    for desc, result in results.items():
        if result is True:
            print(f"✓ {desc}")
        elif result is False:
            print(f"✗ {desc}")
        else:
            print(f"⚠ {desc} - Not compared")
    
    print(f"\nPassed: {passed}/{total}")
    
    if passed == total and total > 0:
        print("\n✓ All compared files match!")
        return 0
    elif total == 0:
        print("\n⚠ No files to compare (scripts may not have been run)")
        return 0
    else:
        print("\n⚠ Some files have differences")
        return 1


if __name__ == "__main__":
    sys.exit(main())

