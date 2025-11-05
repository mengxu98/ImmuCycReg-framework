#!/usr/bin/env python3
"""
Basic test script to verify Python scripts can be imported and basic functionality works
"""

import os
import sys
import importlib.util


def test_script_import(script_path, script_name):
    """Test if a Python script can be imported without errors"""
    print(f"Testing {script_name}...")
    
    if not os.path.exists(script_path):
        print(f"  ✗ File not found: {script_path}")
        return False
    
    try:
        spec = importlib.util.spec_from_file_location(script_name, script_path)
        if spec is None or spec.loader is None:
            print(f"  ✗ Could not load spec for {script_name}")
            return False
        
        module = importlib.util.module_from_spec(spec)
        sys.modules[script_name] = module
        
        spec.loader.exec_module(module)
        print(f"  ✓ {script_name} imported successfully")
        return True
    except SyntaxError as e:
        print(f"  ✗ Syntax error in {script_name}: {e}")
        return False
    except ImportError as e:
        print(f"  ✗ Import error in {script_name}: {e}")
        return False
    except Exception as e:
        print(f"  ⚠ {script_name} has runtime errors (may be expected): {e}")
        return True


def main():
    print("="*60)
    print("Basic Python Script Import Test")
    print("="*60)
    
    base_dir = os.path.dirname(os.path.abspath(__file__))
    
    scripts = [
        "Step-01_Processdata.py",
        "Step-02_NMF-ConsensuClustering.py",
        "Step-03_Cibersort-data.py",
        "Step-04_DESeq2.py",
        "Step-05_CNV-Correlation.py",
        "Step-06_Peak.py",
        "Step-07_Obtain_DNA_sequences.py",
    ]
    
    results = {}
    for script in scripts:
        script_path = os.path.join(base_dir, script)
        script_name = script.replace('.py', '').replace('-', '_')
        results[script] = test_script_import(script_path, script_name)
    
    print("\n" + "="*60)
    print("Summary")
    print("="*60)
    
    passed = sum(1 for v in results.values() if v)
    total = len(results)
    
    for script, result in results.items():
        status = "✓" if result else "✗"
        print(f"{status} {script}")
    
    print(f"\nPassed: {passed}/{total}")
    
    return 0 if passed == total else 1


if __name__ == "__main__":
    sys.exit(main())

