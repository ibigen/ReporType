#!/bin/bash
echo "Starting detect_type installation!"

conda create -n detect_type python=3.8

source activate detect_type

pip install -r requirements.txt

conda env update -f environment.yml

echo "Checking dependencies instalation..."

##Testar instalação de dependências
# Testa a versão do Python
if python3 --version | grep -q "Python 3.8"; then
    echo "Python 3.8 installed correctly"
else
    echo "Error installing Python 3.8"
fi

# Testa a versão do Pandas
if python3 -c "import pandas; print(pandas.__version__)" | grep -q "1.4.2"; then
    echo "Pandas 1.5.3 installed correctly"
else
    echo "Error installing Pandas 1.4.2"
fi

# Testa a versão do NumPy
if python3 -c "import numpy; print(numpy.__version__)" | grep -q "1.22.3"; then
    echo "NumPy 1.22.3 installed correctly"
else
    echo "Error installing NumPy 1.22.3"
fi

# Testa a versão do Biopython
if python3 -c "import Bio; print(Bio.__version__)" | grep -q "1.81"; then
    echo "Biopython 1.81 installed correctly"
else
    echo "Error installing Biopython 1.81"
fi

# Testa a versão do Snakemake
if snakemake --version | grep -q "7.22.0"; then
    echo "Snakemake 7.22.0 installed correctly"
else
    echo "Error installing Snakemake 7.22.0"
fi

# Testa a versão do Abricate
if abricate --version | grep -q "1.0.1"; then
    echo "Abricate 1.0.1 installed correctly"
else
    echo "Error installing Abricate 1.0.1"
fi

# Testa a versão do Raven
if raven --version | grep -q "1.8.1"; then
    echo "Raven 1.8.1 installed correctly"
else
    echo "Error installing Raven 1.8.1"
fi

# Testa a versão do SPAdes
if spades.py --test | grep -q "test passed"; then
    echo "SPAdes 3.13.1 installed correctly"
else
    echo "Error installing SPAdes 3.13.1"
fi

# Testa a versão do Trimmomatic
if trimmomatic -version | grep -q "0.39"; then
    echo "Trimmomatic 0.39 installed correctly"
else
    echo "Error installing Trimmomatic 0.39"
fi

# Testa a versão do EMBOSS
if abiview --version | grep -q "EMBOSS:6.6.0.0"; then
    echo "EMBOSS 6.6.0.0 installed correctly"
else
    echo "Error installing EMBOSS 6.6.0.0"
fi

# Testa a versão do NanoFilt
if NanoFilt --version | grep -q "2.8.0"; then
    echo "NanoFilt 2.8.0 installed correctly"
else
    echo "Error installing NanoFilt 2.8.0"
fi

### Transformar comando "snakemake" em "detect_type"
alias detect_type="conda activate detect_type && snakemake"