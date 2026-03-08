#!/bin/bash
# setup.sh – installs all dependencies for the pipeline on a fresh Ubuntu system

set -e  # exit on error

echo "========================================="
echo " Setting up environment for RNA-seq/WGS pipeline"
echo "========================================="

# 1. Install Miniconda if not present
if ! command -v conda &> /dev/null; then
    echo "Installing Miniconda..."
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
    bash miniconda.sh -b -p $HOME/miniconda
    rm miniconda.sh
    eval "$($HOME/miniconda/bin/conda shell.bash hook)"
    conda init
else
    echo "Conda already installed."
fi

# 2. Create and activate the pipeline environment
ENV_NAME="rnaseq_wgs_pipeline"
if conda env list | grep -q $ENV_NAME; then
    echo "Environment $ENV_NAME already exists. Updating..."
else
    echo "Creating conda environment: $ENV_NAME"
    conda create -y -n $ENV_NAME python=3.9
fi

# Activate environment
source $(conda info --base)/etc/profile.d/conda.sh
conda activate $ENV_NAME

# 3. Install bioinformatics tools via conda
echo "Installing bcftools, samtools, htslib, bedtools..."
conda install -y -c bioconda bcftools samtools htslib bedtools

# 4. Install Python packages
echo "Installing Python packages..."
pip install pyyaml pandas numpy matplotlib seaborn scikit-learn scipy

# 5. Verify installations
echo "Verifying installations:"
bcftools --version | head -n1
bedtools --version
python -c "import sklearn; print('scikit-learn', sklearn.__version__)"

echo "========================================="
echo "Setup complete. To activate the environment, run:"
echo "  conda activate $ENV_NAME"
echo "========================================="