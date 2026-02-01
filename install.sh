#!/bin/bash
#
# Clinical NGS Pipeline - Installation Script
# Maintained by: Prabir
#
# Conda-first, reproducible installation

set -e

ENV_NAME="clinical_NGS_nf"
INSTALL_USER="${SUDO_USER:-$USER}"
CONDA_BIN="/home/${INSTALL_USER}/miniconda3/condabin/conda"

echo "╔════════════════════════════════════════════════════════╗"
echo "║                                                        ║"
echo "║      Clinical NGS Pipeline - Installation              ║"
echo "║               Maintained by: Prabir                    ║"
echo "║                                                        ║"
echo "╚════════════════════════════════════════════════════════╝"
echo ""

############################################
# Root check (system packages only)
############################################
if [[ $EUID -ne 0 ]]; then
   echo "❌ Run this script with sudo:"
   echo "   sudo ./install.sh"
   exit 1
fi

############################################
# System dependencies
############################################
echo "📦 Installing system dependencies..."
apt-get update
apt-get install -y \
    openjdk-11-jdk \
    curl \
    wget \
    git \
    unzip \
    python3 \
    python3-venv

############################################
# Nextflow
############################################
echo ""
echo "🔧 Installing Nextflow..."
if ! command -v nextflow &> /dev/null; then
    curl -s https://get.nextflow.io | bash
    mv nextflow /usr/local/bin/
    chmod +x /usr/local/bin/nextflow
    echo "✅ Nextflow installed"
else
    echo "✅ Nextflow already installed"
fi

############################################
# Conda validation (USER scope)
############################################
echo ""
echo "🐍 Checking Conda for user: ${INSTALL_USER}"

if [[ ! -x "${CONDA_BIN}" ]]; then
    echo "❌ Conda not found at expected location:"
    echo "   ${CONDA_BIN}"
    echo "Fix conda installation before proceeding."
    exit 1
fi

############################################
# Conda environment + Python deps (NON-ROOT)
############################################
echo ""
echo "📦 Creating Conda environment: ${ENV_NAME}"

sudo -u "${INSTALL_USER}" bash << EOF
set -e

export PATH="/home/${INSTALL_USER}/miniconda3/bin:\$PATH"

CONDA_BASE=\$(${CONDA_BIN} info --base)

source "\$CONDA_BASE/etc/profile.d/conda.sh"

if conda env list | awk '{print \$1}' | grep -qx "${ENV_NAME}"; then
    echo "⚠️  Conda environment '${ENV_NAME}' already exists"
else
    conda create -y -n ${ENV_NAME} python=3.10
    echo "✅ Conda environment '${ENV_NAME}' created"
fi

conda activate ${ENV_NAME}

echo "📦 Installing Python dependencies inside ${ENV_NAME}..."
pip install --upgrade pip

pip install \
    flask \
    flask-cors \
    pandas \
    numpy \
    matplotlib \
    seaborn \
    pyyaml

echo "✅ Python dependencies installed"
EOF

############################################
# Docker
############################################
echo ""
echo "🐳 Checking Docker..."
if command -v docker &> /dev/null; then
    echo "✅ Docker is installed"
    echo "Pulling GATK Docker image..."
    docker pull broadinstitute/gatk:latest
else
    echo "⚠️  Docker not installed."
    echo "Install manually: https://docs.docker.com/get-docker/"
fi

############################################
# Project structure
############################################
echo ""
echo "🧪 Creating project directories..."
mkdir -p test_runs results work logs

############################################
# Done
############################################
echo ""
echo "✅ Installation completed successfully!"
echo ""
echo "════════════════════════════════════════════════════════"
echo "Next Steps:"
echo "════════════════════════════════════════════════════════"
echo ""
echo "1. Activate environment:"
echo "   conda activate ${ENV_NAME}"
echo ""
echo "2. Verify pipeline:"
echo "   nextflow run main.nf --help"
echo ""
echo "3. Run test pipeline:"
echo "   nextflow run main.nf \\"
echo "       --analysis_type germline \\"
echo "       --sample_sheet test_data/samples.csv \\"
echo "       --reference_genome /path/to/reference.fa \\"
echo "       --target_bed test_data/target_regions.bed"
echo ""
echo "4. Start web API:"
echo "   cd web && python api_server.py"
echo ""
echo "════════════════════════════════════════════════════════"
echo "Maintained by: Prabir"
echo "════════════════════════════════════════════════════════"
