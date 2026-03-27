#!/usr/bin/env python3
"""
CEL-seq2 Pipeline: Demultiplex by Cell Barcode
This script demultiplexes reads based on the 6bp cell barcode extracted by umi_tools
"""

import gzip
import sys
import os
from collections import defaultdict
import argparse

def read_barcodes(barcode_file):
    """Read barcode file and return dict of barcode -> sample name"""
    barcodes = {}
    with open(barcode_file, 'r') as f:
        # Try to detect header
        first_line = f.readline().strip()
        if 'sample' in first_line.lower() or 'barcode' in first_line.lower():
            # Has header, skip it
            pass
        else:
            # No header, process first line
            f.seek(0)
        
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
                
            parts = line.split(',')
            if len(parts) >= 2:
                sample_name = parts[0].strip()
                barcode = parts[1].strip()
                if barcode and sample_name:
                    barcodes[barcode] = sample_name
    
    return barcodes

def parse_umi_header(header):
    """Extract barcode from UMI-tools modified header
    UMI-tools format: @original_name_BARCODE_UMI additional_info
    """
    # Split by whitespace first to separate read name from other info
    parts = header.split()
    if not parts:
        return None
    
    # The read name with barcode/UMI is the first part
    read_name = parts[0]
    
    # Split by underscore to get components
    name_parts = read_name.split('_')
    
    # UMI-tools adds barcode and UMI as last two components
    if len(name_parts) >= 3:
        # Second to last should be the barcode
        barcode = name_parts[-2]
        return barcode
    
    return None

def hamming_distance(s1, s2):
    """Calculate Hamming distance between two strings"""
    if len(s1) != len(s2):
        return float('inf')
    return sum(c1 != c2 for c1, c2 in zip(s1, s2))

def find_matching_barcode(barcode, barcode_dict, max_mismatches=1):
    """Find matching barcode allowing for mismatches"""
    # First try exact match
    if barcode in barcode_dict:
        return barcode, 0
    
    # If no exact match and mismatches allowed, find closest
    if max_mismatches > 0:
        best_match = None
        best_distance = float('inf')
        
        for ref_barcode in barcode_dict:
            dist = hamming_distance(barcode, ref_barcode)
            if dist <= max_mismatches and dist < best_distance:
                best_match = ref_barcode
                best_distance = dist
        
        if best_match is not None:
            return best_match, best_distance
    
    return None, None

def demultiplex_reads(r1_file, r2_file, barcode_dict, output_dir, max_mismatches=1):
    """Demultiplex paired-end reads based on barcodes"""
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Open file handles for each sample
    r1_handles = {}
    r2_handles = {}
    
    for barcode, sample in barcode_dict.items():
        r1_out = os.path.join(output_dir, f"{sample}_R1.fastq.gz")
        r2_out = os.path.join(output_dir, f"{sample}_R2.fastq.gz")
        r1_handles[barcode] = gzip.open(r1_out, 'wt')
        r2_handles[barcode] = gzip.open(r2_out, 'wt')
    
    # Files for unmatched reads
    r1_unmatched = gzip.open(os.path.join(output_dir, "unmatched_R1.fastq.gz"), 'wt')
    r2_unmatched = gzip.open(os.path.join(output_dir, "unmatched_R2.fastq.gz"), 'wt')
    
    # Statistics
    stats = defaultdict(int)
    total_reads = 0
    
    # Process reads
    print("Processing reads...")
    with gzip.open(r1_file, 'rt') as f1, gzip.open(r2_file, 'rt') as f2:
        while True:
            # Read FASTQ entries
            r1_header = f1.readline().strip()
            r1_seq = f1.readline().strip()
            r1_plus = f1.readline().strip()
            r1_qual = f1.readline().strip()
            
            r2_header = f2.readline().strip()
            r2_seq = f2.readline().strip()
            r2_plus = f2.readline().strip()
            r2_qual = f2.readline().strip()
            
            if not r1_header:
                break
            
            total_reads += 1
            
            # Extract barcode from header
            barcode = parse_umi_header(r1_header)
            
            if not barcode:
                # If we can't parse the barcode, put in unmatched
                r1_unmatched.write(f"{r1_header}\n{r1_seq}\n{r1_plus}\n{r1_qual}\n")
                r2_unmatched.write(f"{r2_header}\n{r2_seq}\n{r2_plus}\n{r2_qual}\n")
                stats['unparseable'] += 1
                continue
            
            # Find matching barcode
            matched_barcode, distance = find_matching_barcode(barcode, barcode_dict, max_mismatches)
            
            if matched_barcode is not None:
                sample = barcode_dict[matched_barcode]
                r1_handles[matched_barcode].write(f"{r1_header}\n{r1_seq}\n{r1_plus}\n{r1_qual}\n")
                r2_handles[matched_barcode].write(f"{r2_header}\n{r2_seq}\n{r2_plus}\n{r2_qual}\n")
                stats[sample] += 1
                if distance > 0:
                    stats['corrected'] += 1
            else:
                r1_unmatched.write(f"{r1_header}\n{r1_seq}\n{r1_plus}\n{r1_qual}\n")
                r2_unmatched.write(f"{r2_header}\n{r2_seq}\n{r2_plus}\n{r2_qual}\n")
                stats['unmatched'] += 1
            
            if total_reads % 100000 == 0:
                print(f"  Processed {total_reads:,} reads...")
    
    # Close all file handles
    for handle in r1_handles.values():
        handle.close()
    for handle in r2_handles.values():
        handle.close()
    r1_unmatched.close()
    r2_unmatched.close()
    
    # Print statistics
    print(f"\n✅ Demultiplexing complete!")
    print(f"Total reads processed: {total_reads:,}")
    
    print(f"\n📊 Reads per sample:")
    for sample in sorted(barcode_dict.values()):
        count = stats[sample]
        percentage = (count / total_reads * 100) if total_reads > 0 else 0
        print(f"  {sample}: {count:,} ({percentage:.2f}%)")
    
    if stats['corrected'] > 0:
        print(f"\n🔧 Reads with corrected barcodes: {stats['corrected']:,}")
    
    if stats['unparseable'] > 0:
        print(f"\n⚠️  Reads with unparseable headers: {stats['unparseable']:,}")
    
    unmatched = stats['unmatched']
    unmatched_percentage = (unmatched / total_reads * 100) if total_reads > 0 else 0
    print(f"\n❌ Unmatched reads: {unmatched:,} ({unmatched_percentage:.2f}%)")
    
    # Return stats for use by other scripts
    return stats, total_reads

def main():
    parser = argparse.ArgumentParser(
        description='Demultiplex CEL-seq2 reads based on cell barcodes',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s --r1 extracted_R1.fastq.gz --r2 extracted_R2.fastq.gz --barcodes barcodes.csv
  
Barcode file format (CSV):
  sample_name,barcode
  Sample1,AAAAAA
  Sample2,CCCCCC
  ...
        """
    )
    
    parser.add_argument('--r1', required=True, help='R1 FASTQ file (with UMIs extracted)')
    parser.add_argument('--r2', required=True, help='R2 FASTQ file (with UMIs extracted)')
    parser.add_argument('--barcodes', required=True, help='CSV file with sample names and barcodes')
    parser.add_argument('--output-dir', default='demultiplexed', help='Output directory (default: demultiplexed)')
    parser.add_argument('--max-mismatches', type=int, default=1, 
                       help='Maximum allowed mismatches in barcode (default: 1)')
    
    args = parser.parse_args()
    
    # Check input files
    for f in [args.r1, args.r2, args.barcodes]:
        if not os.path.exists(f):
            print(f"❌ Error: File not found: {f}")
            sys.exit(1)
    
    # Read barcodes
    print(f"📖 Reading barcodes from {args.barcodes}...")
    barcode_dict = read_barcodes(args.barcodes)
    print(f"   Found {len(barcode_dict)} barcodes")
    
    if not barcode_dict:
        print("❌ Error: No barcodes found in barcode file")
        sys.exit(1)
    
    # Demultiplex
    print(f"\n🔄 Demultiplexing reads...")
    print(f"   Max mismatches allowed: {args.max_mismatches}")
    demultiplex_reads(args.r1, args.r2, barcode_dict, args.output_dir, args.max_mismatches)

if __name__ == '__main__':
    main()
