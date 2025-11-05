#!/usr/bin/env python3
"""
Run complete workflow for ImmuCycReg framework
This script orchestrates the execution of all steps in order
"""

import os
import sys
import subprocess
from pathlib import Path


def check_file_exists(filepath, description):
    """Check if a file exists and print status"""
    if os.path.exists(filepath):
        print(f"✓ {description}: {filepath}")
        return True
    else:
        print(f"✗ {description}: {filepath} - NOT FOUND")
        return False


def run_step(step_name, description, script_path, required_files=None):
    """Run a single step and check prerequisites"""
    print(f"\n{'='*60}")
    print(f"Step: {step_name}")
    print(f"Description: {description}")
    print(f"Script: {script_path}")
    print(f"{'='*60}")
    
    if not os.path.exists(script_path):
        print(f"✗ Script not found: {script_path}")
        return False
    
    # Check required files
    if required_files:
        print("\nChecking prerequisites:")
        all_exist = True
        for file_info in required_files:
            filepath = file_info['path']
            desc = file_info['description']
            if not check_file_exists(filepath, desc):
                all_exist = False
        
        if not all_exist:
            print(f"\n⚠ Some required files are missing. Step may fail.")
            response = input("Continue anyway? (y/n): ").strip().lower()
            if response != 'y':
                return False
    
    # Run the script
    print(f"\nRunning {script_path}...")
    try:
        if script_path.endswith('.py'):
            result = subprocess.run(
                [sys.executable, script_path],
                cwd=os.path.dirname(script_path) or '.',
                check=False,
                capture_output=False
            )
        elif script_path.endswith('.R'):
            result = subprocess.run(
                ['Rscript', script_path],
                cwd=os.path.dirname(script_path) or '.',
                check=False,
                capture_output=False
            )
        else:
            print(f"⚠ Unknown script type: {script_path}")
            return False
        
        if result.returncode == 0:
            print(f"✓ {step_name} completed successfully")
            return True
        else:
            print(f"✗ {step_name} failed with return code: {result.returncode}")
            return False
    except Exception as e:
        print(f"✗ Error running {step_name}: {e}")
        return False


def main():
    print("="*60)
    print("ImmuCycReg Framework - Complete Workflow")
    print("="*60)
    
    base_dir = Path(__file__).parent.parent
    code_dir = base_dir / "code"
    data_dir = base_dir / "data"
    results_dir = base_dir / "results"
    
    # Create directories if they don't exist
    os.makedirs(data_dir, exist_ok=True)
    os.makedirs(results_dir, exist_ok=True)
    
    print(f"\nBase directory: {base_dir}")
    print(f"Code directory: {code_dir}")
    print(f"Data directory: {data_dir}")
    print(f"Results directory: {results_dir}")
    
    # Check initial data files
    print("\n" + "="*60)
    print("Checking Initial Data Files")
    print("="*60)
    
    initial_files = [
        {'path': data_dir / "sample_cluster_4.csv", 'description': 'Sample cluster data'},
        {'path': data_dir / "Genes_17.txt", 'description': 'Genes list'},
        {'path': data_dir / "Genes_249.txt", 'description': 'Genes list (249)'},
    ]
    
    all_exist = True
    for file_info in initial_files:
        if not check_file_exists(str(file_info['path']), file_info['description']):
            all_exist = False
    
    if not all_exist:
        print("\n⚠ Some initial data files are missing.")
        print("Please ensure all required data files are in the data/ directory.")
    
    # Define workflow steps
    workflow_steps = [
        {
            'step': 'Step-01',
            'description': 'Data Processing',
            'script': 'Step-01_Processdata.py',
            'required_files': []
        },
        {
            'step': 'Step-02',
            'description': 'NMF Consensus Clustering',
            'script': 'Step-02_NMF-ConsensuClustering.py',
            'required_files': [
                {'path': results_dir / "TCGA-LUAD.csv", 'description': 'TCGA-LUAD expression data'}
            ]
        },
        {
            'step': 'Step-05',
            'description': 'CNV Correlation',
            'script': 'Step-05_CNV-Correlation.py',
            'required_files': [
                {'path': results_dir / "TCGA-LUAD.csv", 'description': 'TCGA-LUAD expression data'},
                {'path': results_dir / "CNV-LUAD.csv", 'description': 'CNV-LUAD data'}
            ]
        },
    ]
    
    # Run workflow
    print("\n" + "="*60)
    print("Running Workflow")
    print("="*60)
    
    results = {}
    for step_info in workflow_steps:
        script_path = code_dir / step_info['script']
        
        if script_path.exists():
            result = run_step(
                step_info['step'],
                step_info['description'],
                str(script_path),
                step_info.get('required_files')
            )
            results[step_info['step']] = result
        else:
            print(f"\n⚠ Script not found: {step_info['script']}")
            results[step_info['step']] = None
    
    # Summary
    print("\n" + "="*60)
    print("Workflow Summary")
    print("="*60)
    
    passed = sum(1 for v in results.values() if v is True)
    failed = sum(1 for v in results.values() if v is False)
    skipped = sum(1 for v in results.values() if v is None)
    
    for step, result in results.items():
        if result is True:
            print(f"✓ {step}: Success")
        elif result is False:
            print(f"✗ {step}: Failed")
        else:
            print(f"⚠ {step}: Skipped")
    
    print(f"\nPassed: {passed}, Failed: {failed}, Skipped: {skipped}")
    
    if failed == 0:
        print("\n✓ All executed steps completed successfully!")
        return 0
    else:
        print("\n⚠ Some steps failed. Please check the errors above.")
        return 1


if __name__ == "__main__":
    sys.exit(main())

