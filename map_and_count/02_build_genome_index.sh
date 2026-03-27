#!/bin/bash

# CEL-seq2 Pipeline: Build Genome Index
# This script builds the Subread index for C. elegans genome

echo "=== Building Subread Genome Index ==="
echo "===================================="
echo ""

# Check if subread-buildindex is available
if ! command -v subread-buildindex &> /dev/null; then
    echo "❌ Error: subread-buildindex not found"
    echo "   Please run 01_install_dependencies.sh first"
    exit 1
fi

# Set up genome directory
GENOME_DIR="genome_files"
mkdir -p "$GENOME_DIR"

# File paths
FASTA_FILE="$GENOME_DIR/Caenorhabditis_elegans.WBcel235.dna.toplevel.fa"
GTF_FILE="$GENOME_DIR/Caenorhabditis_elegans.WBcel235.104.gtf"
INDEX_PREFIX="$GENOME_DIR/c_elegans_index"

# Check if genome files exist
if [ ! -f "$FASTA_FILE" ]; then
    echo "❌ Error: Genome FASTA not found at $FASTA_FILE"
    echo ""
    echo "📥 Please download the C. elegans genome:"
    echo "   1. Go to: https://www.ensembl.org/Caenorhabditis_elegans/Info/Index"
    echo "   2. Download: Caenorhabditis_elegans.WBcel235.dna.toplevel.fa.gz"
    echo "   3. Extract and place in: $GENOME_DIR/"
    exit 1
fi

if [ ! -f "$GTF_FILE" ]; then
    echo "❌ Error: GTF file not found at $GTF_FILE"
    echo ""
    echo "📥 Please download the C. elegans GTF:"
    echo "   1. Go to: https://www.ensembl.org/Caenorhabditis_elegans/Info/Index"
    echo "   2. Download: Caenorhabditis_elegans.WBcel235.104.gtf.gz"
    echo "   3. Extract and place in: $GENOME_DIR/"
    exit 1
fi

# Check if index already exists
if [ -f "${INDEX_PREFIX}.00.b.array" ]; then
    echo "✅ Index already exists at: $INDEX_PREFIX"
    echo "   To rebuild, delete existing index files first"
    exit 0
fi

# Build index
echo "📊 Building index for: $FASTA_FILE"
echo "📁 Output prefix: $INDEX_PREFIX"
echo ""
echo "🔨 Building index (this may take 5-10 minutes)..."

subread-buildindex -o "$INDEX_PREFIX" "$FASTA_FILE"

# Check if successful
if [ -f "${INDEX_PREFIX}.00.b.array" ]; then
    echo ""
    echo "✅ Index built successfully!"
    echo ""
    echo "Index files created:"
    ls -lh "${INDEX_PREFIX}".*
else
    echo ""
    echo "❌ Index building failed"
    exit 1
fi
