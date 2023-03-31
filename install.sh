#!/bin/bash
echo "Starting detect_type installation!"

conda create -n detect_type python=3.8

echo "Start softwares in conda..."
source activate detect_type

pip install -r requirements.txt && conda install -c bioconda -c conda-forge nanofilt --yes && conda install -c bioconda -c conda-forge spades --yes && conda install -c bioconda -c conda-forge abricate=1.0.1 --yes && conda install -c bioconda -c conda-forge raven-assembler trimmomatic emboss --yes

echo "Checking instalation..."

# Testing Python instalation
if python3 --version | grep -q "Python"; then
    echo "Python correctly"
else
    echo "Error installing Python"
fi

# Testing Pandas instalation
if python3 -c "import pandas; print('pandas', pandas.__version__)" | grep -q "pandas"; then
    echo "Pandas installed correctly"
else
    echo "Error installing Pandas"
fi

# Testing NumPy instalation
if python3 -c "import numpy; print('numpy',numpy.__version__)" | grep -q "numpy"; then
    echo "NumPy installed correctly"
else
    echo "Error installing NumPy"
fi

# Testing Biopython instalation
if python3 -c "import Bio; print('biopython',Bio.__version__)" | grep -q "biopython"; then
    echo "Biopython installed correctly"
else
    echo "Error installing Biopython"
fi

# Testing Snakemake instalation
if command -v snakemake >/dev/null 2>&1 ; then
    echo "Snakemake installed correctly"
else
    echo "Erro installing o Snakemake."
fi

# Testing Abricate instalation
if abricate --version | grep -q "abricate"; then
    echo "Abricate installed correctly"
else
    echo "Error installing Abricate"
fi

# Testing Raven instalation
if command -v raven >/dev/null 2>&1 ; then
    echo "Raven installed correctly"
else
    echo "Erro installing Raven"
fi

# Testing SPAdes instalation
if spades.py --version | grep -q "SPAdes genome assembler"; then
    echo "SPAdes installed correctly"
else
    echo "Error installing SPAdes"
fi

# Testing trimmomatic instalation
if command -v trimmomatic >/dev/null 2>&1 ; then
    echo "Trimmomatic installed correctly"
else
    echo "Erro installing Trimmomatic."
fi

# Testing EMBOSS instalation
if which abiview | grep -q "abiview"; then
    echo "EMBOSS installed correctly"
else
    echo "Error installing EMBOSS"
fi

# Testing NanoFilt instalation
if NanoFilt --version | grep -q "NanoFilt"; then
    echo "NanoFilt installed correctly"
else
    echo "Error installing NanoFilt"
fi
