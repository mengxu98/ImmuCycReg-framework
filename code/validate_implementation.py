#!/usr/bin/env python3
"""
Validate R and Python script implementations for consistency
Compare the logic and structure of R and Python scripts
"""

import os
import sys
import re
from pathlib import Path


def extract_r_logic(script_path):
    """Extract key logic from R script"""
    with open(script_path, 'r') as f:
        content = f.read()
    
    logic = {
        'correlation_method': None,
        'data_loading': [],
        'data_processing': [],
        'output_format': None
    }
    
    # Check correlation method
    if 'corr.test' in content:
        logic['correlation_method'] = 'corr.test'
    elif 'cor(' in content:
        logic['correlation_method'] = 'cor'
    
    # Check data loading
    if 'load(' in content:
        logic['data_loading'].append('RData')
    if 'read.csv' in content or 'read.table' in content:
        logic['data_loading'].append('CSV/TXT')
    
    # Check output
    if 'write.csv' in content:
        logic['output_format'] = 'CSV'
    elif 'save(' in content:
        logic['output_format'] = 'RData'
    
    return logic


def extract_python_logic(script_path):
    """Extract key logic from Python script"""
    with open(script_path, 'r') as f:
        content = f.read()
    
    logic = {
        'correlation_method': None,
        'data_loading': [],
        'data_processing': [],
        'output_format': None
    }
    
    # Check correlation method
    if 'spearmanr' in content:
        logic['correlation_method'] = 'spearmanr'
    elif 'pearsonr' in content:
        logic['correlation_method'] = 'pearsonr'
    elif 'corrcoef' in content:
        logic['correlation_method'] = 'corrcoef'
    
    # Check data loading
    if 'read_csv' in content:
        logic['data_loading'].append('CSV')
    if 'read_table' in content:
        logic['data_loading'].append('TXT')
    
    # Check output
    if 'to_csv' in content:
        logic['output_format'] = 'CSV'
    
    return logic


def compare_scripts(r_script, py_script, description):
    """Compare R and Python script implementations"""
    print(f"\n{'='*60}")
    print(f"Comparing: {description}")
    print(f"{'='*60}")
    
    if not os.path.exists(r_script):
        print(f"✗ R script not found: {r_script}")
        return False
    
    if not os.path.exists(py_script):
        print(f"✗ Python script not found: {py_script}")
        return False
    
    r_logic = extract_r_logic(r_script)
    py_logic = extract_python_logic(py_script)
    
    print(f"\nR Script Logic:")
    print(f"  Correlation method: {r_logic['correlation_method']}")
    print(f"  Data loading: {', '.join(r_logic['data_loading'])}")
    print(f"  Output format: {r_logic['output_format']}")
    
    print(f"\nPython Script Logic:")
    print(f"  Correlation method: {py_logic['correlation_method']}")
    print(f"  Data loading: {', '.join(py_logic['data_loading'])}")
    print(f"  Output format: {py_logic['output_format']}")
    
    # Check consistency
    consistent = True
    
    # Correlation method mapping
    r_to_py = {
        'corr.test': 'spearmanr',
        'cor': 'spearmanr'
    }
    
    if r_logic['correlation_method']:
        expected_py = r_to_py.get(r_logic['correlation_method'], 'spearmanr')
        if py_logic['correlation_method'] != expected_py:
            print(f"\n⚠ Correlation method mismatch:")
            print(f"  R uses: {r_logic['correlation_method']}")
            print(f"  Python uses: {py_logic['correlation_method']}")
            print(f"  Expected Python: {expected_py}")
            consistent = False
    
    # Output format
    if r_logic['output_format'] == 'RData' and py_logic['output_format'] != 'CSV':
        print(f"\n⚠ Output format: R uses RData, Python should use CSV")
        consistent = False
    elif r_logic['output_format'] == 'CSV' and py_logic['output_format'] != 'CSV':
        print(f"\n⚠ Output format mismatch")
        consistent = False
    
    if consistent:
        print(f"\n✓ Scripts are logically consistent")
    else:
        print(f"\n⚠ Some inconsistencies detected")
    
    return consistent


def main():
    print("="*60)
    print("R and Python Script Implementation Validation")
    print("="*60)
    
    base_dir = Path(__file__).parent.parent
    code_dir = base_dir / "code"
    
    scripts_to_compare = [
        {
            'r_script': 'Step-05_CNV-Correlation.R',
            'py_script': 'Step-05_CNV-Correlation.py',
            'description': 'Step 05: CNV Correlation'
        },
    ]
    
    results = {}
    for script_info in scripts_to_compare:
        r_path = code_dir / script_info['r_script']
        py_path = code_dir / script_info['py_script']
        
        result = compare_scripts(
            str(r_path),
            str(py_path),
            script_info['description']
        )
        results[script_info['description']] = result
    
    # Summary
    print("\n" + "="*60)
    print("Validation Summary")
    print("="*60)
    
    passed = sum(1 for v in results.values() if v)
    total = len(results)
    
    for desc, result in results.items():
        status = "✓" if result else "✗"
        print(f"{status} {desc}")
    
    print(f"\nPassed: {passed}/{total}")
    
    if passed == total:
        print("\n✓ All scripts are logically consistent!")
        return 0
    else:
        print("\n⚠ Some inconsistencies detected")
        return 1


if __name__ == "__main__":
    sys.exit(main())

