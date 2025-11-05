#!/bin/bash

echo "=========================================="
echo "Testing R and Python Scripts"
echo "=========================================="

cd "$(dirname "$0")"
BASE_DIR="$(cd .. && pwd)"

echo ""
echo "1. Testing R script syntax..."
Rscript test_r_scripts.R

echo ""
echo "2. Testing Python script syntax..."
python test_basic.py

echo ""
echo "3. Comparing R and Python scripts..."
python compare_r_python.py

echo ""
echo "4. Testing Python package imports..."
python -c "
import sys
packages = ['pandas', 'numpy', 'scipy', 'matplotlib', 'seaborn', 'sklearn']
missing = []
for pkg in packages:
    try:
        __import__(pkg)
        print(f'  ✓ {pkg}')
    except ImportError:
        print(f'  ✗ {pkg} - NOT INSTALLED')
        missing.append(pkg)

if missing:
    print(f'\nMissing packages: {', '.join(missing)}')
    print('Install with: pip install ' + ' '.join(missing))
    sys.exit(1)
else:
    print('\n✓ All core packages installed')
    sys.exit(0)
"

echo ""
echo "5. Testing R package imports..."
Rscript -e "
required <- c('dplyr', 'magrittr', 'DESeq2', 'ConsensusClusterPlus', 'ggplot2')
missing <- c()
for (pkg in required) {
    if (requireNamespace(pkg, quietly = TRUE)) {
        cat(sprintf('  ✓ %s\n', pkg))
    } else {
        cat(sprintf('  ✗ %s - NOT INSTALLED\n', pkg))
        missing <- c(missing, pkg)
    }
}
if (length(missing) > 0) {
    cat(sprintf('\nMissing packages: %s\n', paste(missing, collapse=', ')))
    quit(status = 1)
} else {
    cat('\n✓ All R packages installed\n')
    quit(status = 0)
}
"

echo ""
echo "=========================================="
echo "Test Complete"
echo "=========================================="
