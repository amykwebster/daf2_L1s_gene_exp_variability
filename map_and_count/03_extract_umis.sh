#!/bin/bash

# CEL-seq2 Pipeline: Extract UMIs
# This script extracts UMIs from Read 1 and prepares reads for demultiplexing

echo "=== CEL-seq2 UMI Extraction ==="
echo "==============================="
echo ""

# Check command line arguments
if [ $# -lt 2 ]; then
    echo "Usage: $0 <R1_fastq.gz> <R2_fastq.gz> [output_prefix]"
    echo ""
    echo "Example:"
    echo "  $0 sample_R1.fastq.gz sample_R2.fastq.gz sample_name"
    exit 1
fi

R1_FILE="$1"
R2_FILE="$2"
OUTPUT_PREFIX="${3:-extracted}"

# Check if input files exist
if [ ! -f "$R1_FILE" ]; then
    echo "❌ Error: R1 file not found: $R1_FILE"
    exit 1
fi

if [ ! -f "$R2_FILE" ]; then
    echo "❌ Error: R2 file not found: $R2_FILE"
    exit 1
fi

# Find umi_tools
UMI_TOOLS=""
if command -v umi_tools &> /dev/null; then
    UMI_TOOLS="umi_tools"
elif [ -f "$HOME/Library/Python/3.12/bin/umi_tools" ]; then
    UMI_TOOLS="$HOME/Library/Python/3.12/bin/umi_tools"
else
    echo "❌ Error: umi_tools not found"
    echo "   Please run 01_install_dependencies.sh first"
    exit 1
fi

echo "📁 Input files:"
echo "   R1: $R1_FILE"
echo "   R2: $R2_FILE"
echo ""
echo "🔍 Extracting UMIs and barcodes..."
echo "   Pattern: NNNNNNNCCCCCC (7bp UMI + 6bp barcode)"
echo ""

# Extract UMIs from R1
# Pattern: NNNNNNNCCCCCC means 7 Ns (UMI) followed by 6 Cs (cell barcode)
$UMI_TOOLS extract \
    --bc-pattern=NNNNNNNCCCCCC \
    --stdin="$R1_FILE" \
    --stdout="${OUTPUT_PREFIX}_R1.fastq.gz" \
    --read2-in="$R2_FILE" \
    --read2-out="${OUTPUT_PREFIX}_R2.fastq.gz" \
    --log="${OUTPUT_PREFIX}_umi_extract.log"

# Check if extraction was successful
if [ -f "${OUTPUT_PREFIX}_R1.fastq.gz" ] && [ -f "${OUTPUT_PREFIX}_R2.fastq.gz" ]; then
    echo ""
    echo "✅ UMI extraction complete!"
    echo ""
    echo "📊 Statistics:"
    grep -E "(Input Reads:|Reads output:)" "${OUTPUT_PREFIX}_umi_extract.log"
    echo ""
    echo "📁 Output files:"
    echo "   R1: ${OUTPUT_PREFIX}_R1.fastq.gz"
    echo "   R2: ${OUTPUT_PREFIX}_R2.fastq.gz"
    echo "   Log: ${OUTPUT_PREFIX}_umi_extract.log"
    echo ""
    echo "ℹ️  Next step: Run demultiplexing on the extracted files"
else
    echo ""
    echo "❌ UMI extraction failed. Check the log file for details."
    exit 1
fi
