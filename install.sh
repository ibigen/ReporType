#!/bin/bash
echo "Starting detect_type installation!"

## remove environment if exist
if conda env list | grep "^detect_type" >/dev/null 2>&1; then
    conda remove -n detect_type --all
fi

## create a new environment
conda create -n detect_type python=3.8

echo "Install softwares in conda..."
source activate detect_type && pip install -r requirements.txt && conda install nanofilt --yes && conda install spades --yes && conda install abricate=1.0.1 --yes && conda install raven-assembler trimmomatic emboss --yes;

echo "Checking instalation..."

# Testing Python instalation
if source activate detect_type && python3 --version | grep -q "Python"; then
    echo "Python correctly"
else
    echo "Error installing Python"
fi

# Testing Pandas instalation
if source activate detect_type && python3 -c "import pandas; print('pandas', pandas.__version__)" | grep -q "pandas"; then
    echo "Pandas installed correctly"
else
    echo "Error installing Pandas"
fi

# Testing NumPy instalation
if source activate detect_type &&  python3 -c "import numpy; print('numpy',numpy.__version__)" | grep -q "numpy"; then
    echo "NumPy installed correctly"
else
    echo "Error installing NumPy"
fi

# Testing Biopython instalation
if source activate detect_type && python3 -c "import Bio; print('biopython',Bio.__version__)" | grep -q "biopython"; then
    echo "Biopython installed correctly"
else
    echo "Error installing Biopython"
fi

# Testing Snakemake instalation
if source activate detect_type && command -v snakemake >/dev/null 2>&1 ; then
    echo "Snakemake installed correctly"
else
    echo "Erro installing o Snakemake."
fi

# Testing Abricate instalation
if source activate detect_type && abricate --version | grep -q "abricate"; then
    echo "Abricate installed correctly"
else
    echo "Error installing Abricate"
fi

# Testing Raven instalation
if source activate detect_type && command -v raven >/dev/null 2>&1 ; then
    echo "Raven installed correctly"
else
    echo "Erro installing Raven"
fi

# Testing SPAdes instalation
if source activate detect_type && spades.py --version | grep -q "SPAdes genome assembler"; then
    echo "SPAdes installed correctly"
else
    echo "Error installing SPAdes"
fi

# Testing trimmomatic instalation
if source activate detect_type && command -v trimmomatic >/dev/null 2>&1 ; then
    echo "Trimmomatic installed correctly"
else
    echo "Erro installing Trimmomatic."
fi

# Testing EMBOSS instalation
if source activate detect_type && which abiview | grep -q "abiview"; then
    echo "EMBOSS installed correctly"
else
    echo "Error installing EMBOSS"
fi

# Testing NanoFilt instalation
if source activate detect_type && NanoFilt --version | grep -q "NanoFilt"; then
    echo "NanoFilt installed correctly"
else
    echo "Error installing NanoFilt"
fi

echo ""
echo "To activate the enviroment $ source activate detect_type"
echo ""
