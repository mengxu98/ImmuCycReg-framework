import os
import sys
import subprocess
import pandas as pd
import numpy as np
from pathlib import Path


def test_python_script(script_path, description):
    """Test if a Python script can run without errors"""
    print(f"\n{'='*60}")
    print(f"Testing: {description}")
    print(f"Script: {script_path}")
    print(f"{'='*60}")
    
    if not os.path.exists(script_path):
        print(f"ERROR: Script not found: {script_path}")
        return False
    
    try:
        result = subprocess.run(
            [sys.executable, script_path],
            capture_output=True,
            text=True,
            timeout=300,
            cwd=os.path.dirname(script_path) or '.'
        )
        
        if result.returncode == 0:
            print(f"✓ PASSED: {description}")
            if result.stdout:
                print(f"Output: {result.stdout[-500:]}")
            return True
        else:
            print(f"✗ FAILED: {description}")
            print(f"Return code: {result.returncode}")
            if result.stderr:
                print(f"Error: {result.stderr[-500:]}")
            return False
    except subprocess.TimeoutExpired:
        print(f"✗ TIMEOUT: {description}")
        return False
    except Exception as e:
        print(f"✗ ERROR: {description} - {str(e)}")
        return False


def test_imports():
    """Test if all required Python packages can be imported"""
    print("\n" + "="*60)
    print("Testing Python package imports")
    print("="*60)
    
    required_packages = [
        'pandas', 'numpy', 'scipy', 'matplotlib', 'seaborn',
        'sklearn', 'pydeseq2', 'l0learn'
    ]
    
    failed = []
    for package in required_packages:
        try:
            if package == 'sklearn':
                __import__('sklearn')
            elif package == 'pydeseq2':
                __import__('pydeseq2')
            elif package == 'l0learn':
                __import__('l0learn')
            else:
                __import__(package)
            print(f"✓ {package}")
        except ImportError as e:
            print(f"✗ {package} - {e}")
            failed.append(package)
    
    if failed:
        print(f"\nMissing packages: {', '.join(failed)}")
        print("Please install missing packages:")
        if 'pydeseq2' in failed:
            print("  pip install pydeseq2")
        if 'l0learn' in failed:
            print("  pip install l0learn")
        return False
    return True


def check_data_files():
    """Check if required data files exist"""
    print("\n" + "="*60)
    print("Checking required data files")
    print("="*60)
    
    base_dir = Path(__file__).parent.parent
    data_dir = base_dir / "data"
    results_dir = base_dir / "results"
    
    required_files = [
        data_dir / "sample_cluster_4.csv",
        data_dir / "Genes_17.txt",
        data_dir / "Genes_249.txt",
    ]
    
    missing = []
    for file_path in required_files:
        if file_path.exists():
            print(f"✓ {file_path.name}")
        else:
            print(f"✗ {file_path.name} - NOT FOUND")
            missing.append(str(file_path))
    
    if missing:
        print(f"\nMissing data files: {len(missing)}")
        print("Some scripts may fail without these files.")
        return False
    return True


def main():
    print("="*60)
    print("Testing Python Scripts")
    print("="*60)
    
    base_dir = Path(__file__).parent
    code_dir = base_dir
    
    all_passed = True
    
    if not test_imports():
        all_passed = False
    
    if not check_data_files():
        all_passed = False
    
    scripts_to_test = [
        ("Step-01_Processdata.py", "Step 01: Data Processing"),
        ("Step-02_NMF-ConsensuClustering.py", "Step 02: NMF Consensus Clustering"),
        ("Step-03_Cibersort-data.py", "Step 03: CIBERSORT Data Preparation"),
        ("Step-04_DESeq2.py", "Step 04: DESeq2 Analysis"),
        ("Step-05_CNV-Correlation.py", "Step 05: CNV Correlation"),
        ("Step-06_Peak.py", "Step 06: Peak Analysis"),
        ("Step-07_Obtain_DNA_sequences.py", "Step 07: DNA Sequence Extraction"),
    ]
    
    results = {}
    for script_name, description in scripts_to_test:
        script_path = code_dir / script_name
        if script_path.exists():
            passed = test_python_script(str(script_path), description)
            results[script_name] = passed
            if not passed:
                all_passed = False
        else:
            print(f"\n⚠ SKIPPED: {script_name} - File not found")
            results[script_name] = None
    
    print("\n" + "="*60)
    print("Test Summary")
    print("="*60)
    for script_name, result in results.items():
        if result is True:
            print(f"✓ {script_name}")
        elif result is False:
            print(f"✗ {script_name}")
        else:
            print(f"⚠ {script_name} - Not tested")
    
    if all_passed:
        print("\n✓ All tests passed!")
        return 0
    else:
        print("\n✗ Some tests failed. Please check the errors above.")
        return 1


if __name__ == "__main__":
    sys.exit(main())

