#!/bin/bash

# CEL-seq2 Pipeline: Merge Multiple Lanes
# This script merges L001, L002, etc. files from the same pool

echo "=== CEL-seq2 Lane Merger ==="
echo "============================"
echo ""

# Check arguments
if [ $# -lt 1 ]; then
    echo "Usage: $0 <pool_name> [input_dir] [output_dir]"
    echo ""
    echo "Example:"
    echo "  $0 pool1                    # Looks in current directory"
    echo "  $0 pool1 ../raw_data/pool1  # Specify input directory"
    echo ""
    echo "This will merge:"
    echo "  pool1_L001_R1.fastq.gz + pool1_L002_R1.fastq.gz → pool1_merged_R1.fastq.gz"
    echo "  pool1_L001_R2.fastq.gz + pool1_L002_R2.fastq.gz → pool1_merged_R2.fastq.gz"
    exit 1
fi

POOL_NAME="$1"
INPUT_DIR="${2:-.}"
OUTPUT_DIR="${3:-$INPUT_DIR}"

# Find all R1 files for this pool
R1_FILES=($(ls "$INPUT_DIR"/${POOL_NAME}_L*_R1*.fastq.gz 2>/dev/null | sort))
R2_FILES=($(ls "$INPUT_DIR"/${POOL_NAME}_L*_R2*.fastq.gz 2>/dev/null | sort))

if [ ${#R1_FILES[@]} -eq 0 ]; then
    echo "❌ Error: No files found matching ${POOL_NAME}_L*_R1*.fastq.gz in $INPUT_DIR"
    exit 1
fi

echo "📊 Found ${#R1_FILES[@]} lanes for $POOL_NAME"
echo ""

# If only one lane, just create symlinks
if [ ${#R1_FILES[@]} -eq 1 ]; then
    echo "ℹ️  Only one lane found, creating symlinks..."
    ln -sf "$(basename ${R1_FILES[0]})" "$OUTPUT_DIR/${POOL_NAME}_merged_R1.fastq.gz"
    ln -sf "$(basename ${R2_FILES[0]})" "$OUTPUT_DIR/${POOL_NAME}_merged_R2.fastq.gz"
    echo "✅ Created symlinks for single lane"
    exit 0
fi

# Show files to be merged
echo "Files to merge:"
echo "R1 files:"
for f in "${R1_FILES[@]}"; do
    echo "  - $(basename $f)"
done
echo ""
echo "R2 files:"
for f in "${R2_FILES[@]}"; do
    echo "  - $(basename $f)"
done
echo ""

# Count reads in each lane
echo "📊 Counting reads per lane..."
for i in "${!R1_FILES[@]}"; do
    r1_file="${R1_FILES[$i]}"
    lane=$(basename "$r1_file" | grep -o 'L[0-9]\+')
    read_count=$(zcat "$r1_file" | echo $(wc -l) / 4 | bc)
    echo "  $lane: $read_count reads"
done
echo ""

# Merge files
OUTPUT_R1="$OUTPUT_DIR/${POOL_NAME}_merged_R1.fastq.gz"
OUTPUT_R2="$OUTPUT_DIR/${POOL_NAME}_merged_R2.fastq.gz"

echo "🔄 Merging R1 files..."
cat "${R1_FILES[@]}" > "$OUTPUT_R1"

echo "🔄 Merging R2 files..."
cat "${R2_FILES[@]}" > "$OUTPUT_R2"

# Verify merge
echo ""
echo "📊 Verifying merged files..."
merged_r1_count=$(zcat "$OUTPUT_R1" | echo $(wc -l) / 4 | bc)
merged_r2_count=$(zcat "$OUTPUT_R2" | echo $(wc -l) / 4 | bc)

echo "  Merged R1: $merged_r1_count reads"
echo "  Merged R2: $merged_r2_count reads"

if [ "$merged_r1_count" -eq "$merged_r2_count" ]; then
    echo "  ✅ Read counts match"
else
    echo "  ⚠️  Warning: R1 and R2 have different read counts!"
fi

echo ""
echo "✅ Lane merging complete!"
echo ""
echo "📁 Merged files:"
echo "  R1: $OUTPUT_R1"
echo "  R2: $OUTPUT_R2"
echo ""
echo "🚀 You can now run the pipeline with:"
echo "  ./run_celseq2_pipeline.sh \\"
echo "    -r1 $OUTPUT_R1 \\"
echo "    -r2 $OUTPUT_R2 \\"
echo "    -b $INPUT_DIR/barcodes.csv \\"
echo "    -p $POOL_NAME"
