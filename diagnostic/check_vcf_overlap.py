#!/usr/bin/env python3
"""
check_vcf_overlap.py

A standalone script to compare positions between two VCF files (e.g., RNA and WGS).
It extracts all chromosome and position pairs from each VCF and computes the overlap,
helping to identify issues such as chromosome naming mismatches or different genome builds.

Usage:
    python check_vcf_overlap.py --rna RNA.vcf.gz --wgs WGS.vcf.gz [--chroms CHR1,CHR2] [--max POS] --output DIR

Options:
    --rna FILE      RNA VCF file (bgzipped)
    --wgs FILE      WGS VCF file (bgzipped)
    --chroms LIST   Comma-separated list of chromosomes to analyze (e.g., chr1,chr2). If omitted, uses all common contigs.
    --max N         Maximum number of positions to load per VCF (for very large files, set to e.g. 10_000_000)
    --output DIR    Output directory (will be created)
"""

import os
import sys
import argparse
import subprocess
import pandas as pd
import tempfile
import shutil
from collections import Counter

def run_cmd(cmd, check=True):
    """Run a shell command and return stdout."""
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if check and result.returncode != 0:
        raise RuntimeError(f"Command failed: {cmd}\n{result.stderr}")
    return result.stdout.strip()

def get_vcf_contigs(vcf):
    """Return sorted list of contigs from VCF index."""
    cmd = f"bcftools index -s {vcf} | cut -f1 | sort -V"
    result = run_cmd(cmd, check=False)
    if result:
        return result.split('\n')
    else:
        # fallback: extract from header
        cmd = f"bcftools view -h {vcf} | grep '^##contig=<ID=' | sed 's/^##contig=<ID=\\([^,]*\\),.*$/\\1/'"
        result = run_cmd(cmd)
        return result.split('\n')

def ensure_vcf_index(vcf):
    """Create index if missing."""
    if not os.path.exists(vcf + '.tbi'):
        print(f"Indexing {vcf}...")
        run_cmd(f"bcftools index {vcf}")

def load_positions(vcf, chroms, max_pos=None):
    """
    Load all (chrom, pos) pairs from VCF for given chromosomes.
    Returns a set of tuples.
    """
    positions = set()
    total = 0
    for chrom in chroms:
        if max_pos and total >= max_pos:
            break
        cmd = f"bcftools query -f '%CHROM\t%POS\n' -r {chrom} {vcf}"
        proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        for line in proc.stdout:
            if max_pos and total >= max_pos:
                proc.terminate()
                break
            line = line.strip()
            if line:
                parts = line.split('\t')
                if len(parts) >= 2:
                    positions.add((parts[0], int(parts[1])))
                    total += 1
        proc.wait()
    return positions

def main():
    parser = argparse.ArgumentParser(description="Compare positions between two VCFs")
    parser.add_argument("--rna", required=True, help="RNA VCF file")
    parser.add_argument("--wgs", required=True, help="WGS VCF file")
    parser.add_argument("--chroms", help="Comma-separated list of chromosomes (e.g., chr1,chr2). Default: common contigs")
    parser.add_argument("--max", type=int, default=0, help="Maximum number of positions to load per VCF (0 = all)")
    parser.add_argument("--output", required=True, help="Output directory")
    args = parser.parse_args()

    out_dir = args.output
    os.makedirs(out_dir, exist_ok=True)

    print("="*60)
    print("VCF Position Overlap Checker")
    print("="*60)

    # Ensure VCFs are indexed
    ensure_vcf_index(args.rna)
    ensure_vcf_index(args.wgs)

    # Get contigs
    rna_contigs_all = get_vcf_contigs(args.rna)
    wgs_contigs_all = get_vcf_contigs(args.wgs)

    if args.chroms:
        desired = [c.strip() for c in args.chroms.split(',')]
    else:
        # Use common contigs between the two VCFs
        desired = sorted(set(rna_contigs_all) & set(wgs_contigs_all))
        if not desired:
            print("No common contigs found. Trying to infer naming...")
            # Check if one has "chr" and other doesn't
            rna_has_chr = any(c.startswith('chr') for c in rna_contigs_all)
            wgs_has_chr = any(c.startswith('chr') for c in wgs_contigs_all)
            if rna_has_chr != wgs_has_chr:
                print("Chromosome naming mismatch detected.")
                # Propose conversion
                if rna_has_chr:
                    # RNA has "chr", WGS does not
                    rna_stripped = {c.replace('chr', '') for c in rna_contigs_all if c.startswith('chr')}
                    common = rna_stripped & set(wgs_contigs_all)
                    if common:
                        print(f"After stripping 'chr' from RNA, {len(common)} common contigs found.")
                        print("Example: RNA 'chr1' would match WGS '1'.")
                        # Use stripped names for RNA? We'll handle in load with conversion.
                        # For simplicity, we'll load RNA with original names and then convert during intersection.
                    else:
                        print("No matches even after stripping.")
                else:
                    # WGS has "chr", RNA does not
                    wgs_stripped = {c.replace('chr', '') for c in wgs_contigs_all if c.startswith('chr')}
                    common = wgs_stripped & set(rna_contigs_all)
                    if common:
                        print(f"After stripping 'chr' from WGS, {len(common)} common contigs found.")
                        print("Example: WGS 'chr1' would match RNA '1'.")
            # Fallback to using all contigs from RNA? But we'll just exit.
            sys.exit("Unable to determine chromosomes to compare. Please specify --chroms manually.")

    print(f"Chromosomes to analyze: {', '.join(desired)}")

    # Load positions
    print(f"Loading positions from RNA VCF...")
    rna_pos = load_positions(args.rna, desired, args.max)
    print(f"  Loaded {len(rna_pos)} RNA positions.")

    print(f"Loading positions from WGS VCF...")
    wgs_pos = load_positions(args.wgs, desired, args.max)
    print(f"  Loaded {len(wgs_pos)} WGS positions.")

    # Original intersection
    intersect = rna_pos & wgs_pos
    print(f"Overlap (original names): {len(intersect)} positions")

    # Try stripping "chr" from RNA if no overlap
    if len(intersect) == 0 and len(rna_pos) > 0 and len(wgs_pos) > 0:
        print("No overlap with original names. Attempting to fix chromosome naming...")
        # Check if RNA has "chr" and WGS doesn't
        rna_sample = next(iter(rna_pos))[0]
        wgs_sample = next(iter(wgs_pos))[0]
        if rna_sample.startswith('chr') and not wgs_sample.startswith('chr'):
            # Strip "chr" from RNA
            rna_stripped = {(c.replace('chr', ''), p) for c, p in rna_pos if c.startswith('chr')}
            intersect_stripped = rna_stripped & wgs_pos
            print(f"  After stripping 'chr' from RNA: {len(intersect_stripped)} positions")
        elif not rna_sample.startswith('chr') and wgs_sample.startswith('chr'):
            # Add "chr" to RNA
            rna_added = {('chr' + c, p) for c, p in rna_pos}
            intersect_added = rna_added & wgs_pos
            print(f"  After adding 'chr' to RNA: {len(intersect_added)} positions")
        else:
            print("  No simple naming mismatch detected.")

    # Save summary
    with open(os.path.join(out_dir, 'overlap_summary.txt'), 'w') as f:
        f.write(f"RNA positions: {len(rna_pos)}\n")
        f.write(f"WGS positions: {len(wgs_pos)}\n")
        f.write(f"Overlap (original): {len(intersect)}\n")
        f.write(f"Overlap fraction (RNA): {len(intersect)/len(rna_pos) if rna_pos else 0:.6f}\n")
        f.write(f"Overlap fraction (WGS): {len(intersect)/len(wgs_pos) if wgs_pos else 0:.6f}\n")

    # Save sample of overlapping positions (if any)
    if len(intersect) > 0:
        sample = list(intersect)[:1000]
        with open(os.path.join(out_dir, 'overlap_sample.tsv'), 'w') as f:
            f.write("chrom\tpos\n")
            for chrom, pos in sorted(sample):
                f.write(f"{chrom}\t{pos}\n")
        print(f"Sample of overlapping positions saved to {out_dir}/overlap_sample.tsv")
    else:
        print("No overlapping positions found.")

    print(f"Summary saved to {out_dir}/overlap_summary.txt")
    print("Done.")

if __name__ == "__main__":
    main()