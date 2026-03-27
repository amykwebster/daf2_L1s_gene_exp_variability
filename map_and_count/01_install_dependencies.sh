#!/bin/bash

# CEL-seq2 Pipeline: Install Dependencies
# This script installs the necessary tools for the CEL-seq2 pipeline

echo "=== Installing CEL-seq2 Pipeline Dependencies ==="
echo "================================================"
echo ""

# Check if running on macOS
if [[ "$OSTYPE" != "darwin"* ]]; then
    echo "❌ This script is designed for macOS. Please adapt for your system."
    exit 1
fi

# Function to check if a command exists
command_exists() {
    command -v "$1" &> /dev/null
}

# 1. Check/Install Python packages
echo "📦 Installing Python packages..."
pip3 install --user umi_tools pandas numpy

# 2. Install Subread if not present
if ! command_exists subread-align; then
    echo ""
    echo "📦 Installing Subread..."
    mkdir -p tools
    cd tools
    
    # Download Subread for macOS
    curl -L https://sourceforge.net/projects/subread/files/subread-2.0.6/subread-2.0.6-macOS-x86_64.tar.gz/download -o subread-2.0.6-macOS-x86_64.tar.gz
    tar -xzf subread-2.0.6-macOS-x86_64.tar.gz
    
    # Create local bin directory
    cd ..
    mkdir -p bin
    
    # Create symlinks
    for tool in tools/subread-2.0.6-macOS-x86_64/bin/*; do
        ln -sf "../$tool" bin/
    done
    
    echo "✅ Subread installed in tools/"
    echo ""
    echo "⚠️  Add this to your ~/.zshrc or ~/.bashrc:"
    echo "export PATH=\"$PWD/bin:\$PATH\""
else
    echo "✅ Subread already installed"
fi

# 3. Check for samtools
if ! command_exists samtools; then
    echo ""
    echo "⚠️  samtools not found. Please install with:"
    echo "    brew install samtools"
    echo "    or"
    echo "    conda install -c bioconda samtools"
fi

# 4. Check umi_tools installation
echo ""
echo "🔍 Checking umi_tools installation..."
if python3 -c "import umi_tools" 2>/dev/null; then
    echo "✅ umi_tools is installed"
    UMI_PATH=$(python3 -c "import os; import umi_tools; print(os.path.dirname(umi_tools.__file__))")
    echo "   Location: $UMI_PATH"
else
    echo "❌ umi_tools import failed. Trying to find executable..."
    if command_exists umi_tools; then
        echo "✅ umi_tools executable found at: $(which umi_tools)"
    else
        # Try common locations
        if [ -f "$HOME/Library/Python/3.12/bin/umi_tools" ]; then
            echo "✅ umi_tools found at: $HOME/Library/Python/3.12/bin/umi_tools"
            echo "⚠️  Add to PATH: export PATH=\"\$HOME/Library/Python/3.12/bin:\$PATH\""
        else
            echo "❌ umi_tools not found. Please ensure it's installed and in PATH"
        fi
    fi
fi

echo ""
echo "✅ Dependency check complete!"
echo ""
echo "📋 Summary:"
echo "   - Python packages: umi_tools, pandas, numpy"
echo "   - Subread: mapping and counting"
echo "   - samtools: BAM file manipulation"
echo ""
echo "🔧 Before running the pipeline, ensure your PATH includes:"
echo "   - $PWD/bin (for Subread)"
echo "   - Python bin directory (for umi_tools)"
