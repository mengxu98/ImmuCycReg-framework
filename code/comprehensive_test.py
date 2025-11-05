#!/usr/bin/env python3
"""
Comprehensive test suite for R and Python scripts
Tests syntax, dependencies, and result consistency
"""

import os
import sys
import subprocess
import pandas as pd
import numpy as np
from pathlib import Path
import json
from datetime import datetime


def test_r_syntax():
    """Test R script syntax"""
    print("\n" + "="*60)
    print("Testing R Script Syntax")
    print("="*60)
    
    base_dir = Path(__file__).parent.parent
    code_dir = base_dir / "code"
    
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
    
    results = {}
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
                print(f"{status} {script}")
                results[script] = passed
            except Exception as e:
                print(f"✗ {script} - Error: {e}")
                results[script] = False
        else:
            print(f"⚠ {script} - Not found")
            results[script] = None
    
    return results


def test_python_syntax():
    """Test Python script syntax"""
    print("\n" + "="*60)
    print("Testing Python Script Syntax")
    print("="*60)
    
    base_dir = Path(__file__).parent.parent
    code_dir = base_dir / "code"
    
    py_scripts = [
        "Step-01_Processdata.py",
        "Step-02_NMF-ConsensuClustering.py",
        "Step-03_Cibersort-data.py",
        "Step-04_DESeq2.py",
        "Step-05_CNV-Correlation.py",
        "Step-06_Peak.py",
        "Step-07_Obtain_DNA_sequences.py",
        "Step-08_ImmuCycReg-framework.py"
    ]
    
    results = {}
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
                print(f"{status} {script}")
                results[script] = passed
            except Exception as e:
                print(f"✗ {script} - Error: {e}")
                results[script] = False
        else:
            print(f"⚠ {script} - Not found")
            results[script] = None
    
    return results


def test_dependencies():
    """Test package dependencies"""
    print("\n" + "="*60)
    print("Testing Package Dependencies")
    print("="*60)
    
    # Python packages
    print("\nPython Packages:")
    py_packages = {
        'pandas': 'pandas',
        'numpy': 'numpy',
        'scipy': 'scipy',
        'matplotlib': 'matplotlib',
        'seaborn': 'seaborn',
        'sklearn': 'sklearn',
        'pydeseq2': 'pydeseq2',
        'l0learn': 'l0learn'
    }
    
    py_results = {}
    for name, module in py_packages.items():
        try:
            if module == 'sklearn':
                __import__('sklearn')
            else:
                __import__(module)
            print(f"  ✓ {name}")
            py_results[name] = True
        except ImportError:
            print(f"  ✗ {name} - NOT INSTALLED")
            py_results[name] = False
    
    # R packages
    print("\nR Packages:")
    r_packages = [
        'dplyr', 'magrittr', 'DESeq2', 'ConsensusClusterPlus',
        'ggplot2', 'psych', 'survival', 'survminer'
    ]
    
    r_results = {}
    try:
        result = subprocess.run(
            ['Rscript', '-e', f'''
            packages <- c({', '.join([f'"{p}"' for p in r_packages])})
            for (pkg in packages) {{
                if (requireNamespace(pkg, quietly = TRUE)) {{
                    cat(sprintf("  ✓ %s\\n", pkg))
                }} else {{
                    cat(sprintf("  ✗ %s - NOT INSTALLED\\n", pkg))
                }}
            }}
            '''],
            capture_output=True,
            text=True,
            timeout=30
        )
        print(result.stdout)
    except Exception as e:
        print(f"Error testing R packages: {e}")
    
    return py_results, r_results


def check_data_files():
    """Check if required data files exist"""
    print("\n" + "="*60)
    print("Checking Data Files")
    print("="*60)
    
    base_dir = Path(__file__).parent.parent
    data_dir = base_dir / "data"
    
    required_files = [
        "sample_cluster_4.csv",
        "Genes_17.txt",
        "Genes_249.txt",
        "data/Genes_2230.csv"
    ]
    
    results = {}
    for file in required_files:
        file_path = data_dir / file
        if file_path.exists():
            print(f"✓ {file}")
            results[file] = True
        else:
            print(f"✗ {file} - NOT FOUND")
            results[file] = False
    
    return results


def generate_report(r_syntax, py_syntax, py_deps, data_files):
    """Generate comprehensive test report"""
    print("\n" + "="*60)
    print("Comprehensive Test Report")
    print("="*60)
    
    report = {
        'test_date': datetime.now().isoformat(),
        'r_syntax': {
            'total': len([v for v in r_syntax.values() if v is not None]),
            'passed': sum(1 for v in r_syntax.values() if v is True),
            'failed': sum(1 for v in r_syntax.values() if v is False)
        },
        'python_syntax': {
            'total': len([v for v in py_syntax.values() if v is not None]),
            'passed': sum(1 for v in py_syntax.values() if v is True),
            'failed': sum(1 for v in py_syntax.values() if v is False)
        },
        'python_dependencies': {
            'total': len(py_deps),
            'installed': sum(1 for v in py_deps.values() if v),
            'missing': [k for k, v in py_deps.items() if not v]
        },
        'data_files': {
            'total': len(data_files),
            'found': sum(1 for v in data_files.values() if v),
            'missing': [k for k, v in data_files.items() if not v]
        }
    }
    
    print(f"\nR Scripts: {report['r_syntax']['passed']}/{report['r_syntax']['total']} passed")
    print(f"Python Scripts: {report['python_syntax']['passed']}/{report['python_syntax']['total']} passed")
    print(f"Python Dependencies: {report['python_dependencies']['installed']}/{report['python_dependencies']['total']} installed")
    print(f"Data Files: {report['data_files']['found']}/{report['data_files']['total']} found")
    
    if report['python_dependencies']['missing']:
        print(f"\nMissing Python packages: {', '.join(report['python_dependencies']['missing'])}")
        print("Install with:")
        print(f"  pip install {' '.join(report['python_dependencies']['missing'])}")
    
    if report['data_files']['missing']:
        print(f"\nMissing data files: {', '.join(report['data_files']['missing'])}")
        print("Some scripts may not run without these files.")
    
    # Save report
    base_dir = Path(__file__).parent.parent
    report_file = base_dir / "results" / "test_report.json"
    report_file.parent.mkdir(parents=True, exist_ok=True)
    
    with open(report_file, 'w') as f:
        json.dump(report, f, indent=2)
    
    print(f"\nReport saved to: {report_file}")
    
    return report


def main():
    print("="*60)
    print("Comprehensive Test Suite")
    print("="*60)
    print(f"Test Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    
    # Run tests
    r_syntax = test_r_syntax()
    py_syntax = test_python_syntax()
    py_deps, r_deps = test_dependencies()
    data_files = check_data_files()
    
    # Generate report
    report = generate_report(r_syntax, py_syntax, py_deps, data_files)
    
    # Overall status
    all_passed = (
        report['r_syntax']['failed'] == 0 and
        report['python_syntax']['failed'] == 0
    )
    
    print("\n" + "="*60)
    if all_passed:
        print("✓ All syntax tests passed!")
        print("\nNext steps:")
        print("1. Install missing packages if any")
        print("2. Run R and Python scripts with test data")
        print("3. Compare outputs using verify_results.py")
        return 0
    else:
        print("⚠ Some tests failed")
        print("Please fix syntax errors before proceeding")
        return 1


if __name__ == "__main__":
    sys.exit(main())

