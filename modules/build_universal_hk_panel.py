#!/usr/bin/env python3
"""
build_universal_hk_panel.py

A standalone script to create a universal SNP panel from housekeeping gene exons
that are present in both RNA‑seq and WGS VCFs. The output TSV (chrom, pos, ref, alt)
can be used in the pipeline's study_snp_analysis module.

Requirements:
  - bcftools (in PATH)
  - bedtools (in PATH)
  - Python packages: pandas, requests (optional, for download)

Usage:
  python build_universal_hk_panel.py \
      --rna_vcf RNA.vcf.gz \
      --wgs_vcf WGS.vcf.gz \
      --gtf gencode.v38.annotation.gtf \
      --output universal_hk_snps.tsv
"""

import os
import sys
import subprocess
import tempfile
import argparse
import pandas as pd
import gzip
from collections import defaultdict

# ------------------------------------------------------------------------------
# Utility functions
# ------------------------------------------------------------------------------
def run_cmd(cmd, description="", check=True):
    """Run shell command and return stdout."""
    print(f"  [CMD] {cmd}")
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if check and result.returncode != 0:
        raise RuntimeError(f"{description} failed:\n{result.stderr}")
    return result.stdout.strip()

def ensure_vcf_index(vcf):
    """Index VCF if missing."""
    if not os.path.exists(vcf + '.tbi'):
        print(f"  Indexing {vcf}...")
        run_cmd(f"bcftools index {vcf}", "Indexing VCF")

def download_hk_genes(output_file):
    """
    Download housekeeping gene list from HRT Atlas (via MSigDB).
    Fallback to a curated list if download fails.
    """
    print("Downloading housekeeping gene list...")
    url = "https://www.gsea-msigdb.org/gsea/msigdb/download_geneset.jsp?geneSetName=HOUNKPE_HOUSEKEEPING_GENES&fileType=txt"
    try:
        import requests
        response = requests.get(url, timeout=30)
        response.raise_for_status()
        genes = [line.strip() for line in response.text.splitlines() if line and not line.startswith('#')]
        with open(output_file, 'w') as f:
            f.write('\n'.join(genes))
        print(f"  Downloaded {len(genes)} housekeeping genes.")
        return genes
    except Exception as e:
        print(f"  Download failed: {e}")
        print("  Using fallback list of common housekeeping genes.")
        return fallback_hk_genes(output_file)

def fallback_hk_genes(output_file):
    """Curated list of well‑known housekeeping genes."""
    genes = [
        "ACTB", "GAPDH", "RPLP0", "B2M", "HPRT1", "TBP", "UBC",
        "YWHAZ", "PPIA", "PGK1", "ALDOA", "RPS27A", "RPL19",
        "RPL11", "NONO", "ARHGDIA", "RPL32", "RPS18", "HSP90AB1",
        "ATP5F1", "CYC1", "SDHA", "GUSB", "HMBS", "IPO8", "POLR2A"
    ]
    with open(output_file, 'w') as f:
        f.write('\n'.join(genes))
    print(f"  Using fallback list of {len(genes)} housekeeping genes.")
    return genes

# ------------------------------------------------------------------------------
# Extract exon coordinates from GTF for given genes (pure Python)
# ------------------------------------------------------------------------------
def extract_exons_python(gtf_file, gene_list_file, output_bed):
    """
    Read GTF line by line, extract exons for genes in gene_list_file.
    Output BED6: chrom, start, end, gene, score, strand.
    """
    print("Extracting exons for housekeeping genes from GTF using Python...")
    
    # Read gene list into a set
    with open(gene_list_file) as f:
        genes = {line.strip() for line in f if line.strip()}
    
    exon_count = 0
    opener = gzip.open if gtf_file.endswith('.gz') else open
    mode = 'rt' if gtf_file.endswith('.gz') else 'r'
    
    with opener(gtf_file, mode) as gtf, open(output_bed, 'w') as out:
        for line in gtf:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue
            chrom, source, feature, start, end, score, strand, frame, attributes = parts
            if feature != 'exon':
                continue
            
            # Parse attributes to find gene_name
            # Attributes are key-value pairs separated by ';' and often spaces
            attr_dict = {}
            for attr in attributes.strip().split(';'):
                if not attr:
                    continue
                # Split on first space or '='
                if ' ' in attr:
                    key, val = attr.strip().split(' ', 1)
                elif '=' in attr:
                    key, val = attr.strip().split('=', 1)
                else:
                    continue
                # Remove quotes from value
                val = val.strip().strip('"')
                attr_dict[key] = val
            
            gene_name = attr_dict.get('gene_name')
            if gene_name and gene_name in genes:
                # Write BED6: chrom, start-1, end, gene, score, strand
                out.write(f"{chrom}\t{int(start)-1}\t{end}\t{gene_name}\t.\t{strand}\n")
                exon_count += 1
    
    print(f"  Found {exon_count} exons for housekeeping genes.")
    return exon_count

# ------------------------------------------------------------------------------
# Extract positions from a VCF (using bcftools query)
# ------------------------------------------------------------------------------
def extract_positions(vcf, chroms=None):
    """
    Extract all (chrom, pos) pairs from VCF, optionally restricted to chromosomes.
    Returns a set of tuples.
    """
    cmd = "bcftools query -f '%CHROM\t%POS\n'"
    if chroms:
        chrom_list = ','.join(chroms)
        cmd += f" -r {chrom_list}"
    cmd += f" {vcf}"
    result = run_cmd(cmd, f"Extracting positions from {vcf}")
    positions = set()
    for line in result.split('\n'):
        if not line:
            continue
        chrom, pos = line.split()
        positions.add((chrom, int(pos)))
    return positions

# ------------------------------------------------------------------------------
# Get alleles for positions from VCF
# ------------------------------------------------------------------------------
def get_alleles(vcf, positions, output_tsv):
    """
    For a given set of (chrom,pos), retrieve REF and ALT from VCF.
    Writes a TSV with columns: chrom, pos, ref, alt.
    """
    # Create a temporary BED file
    with tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False) as f:
        bed_file = f.name
        for chrom, pos in sorted(positions):
            f.write(f"{chrom}\t{pos-1}\t{pos}\n")
    # Use bcftools view to extract those positions and query alleles
    query_file = tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False).name
    cmd = f"bcftools view -R {bed_file} {vcf} | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' > {query_file}"
    run_cmd(cmd, "Extracting alleles")
    # Read alleles and remove duplicates
    df = pd.read_csv(query_file, sep='\t', header=None, names=['chrom','pos','ref','alt'])
    df.drop_duplicates(inplace=True)
    df.to_csv(output_tsv, sep='\t', index=False, header=False)
    # Cleanup
    os.unlink(bed_file)
    os.unlink(query_file)
    return len(df)

# ------------------------------------------------------------------------------
# Main
# ------------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(description="Build universal housekeeping SNP panel from RNA and WGS VCFs")
    parser.add_argument("--rna_vcf", required=True, help="RNA-seq VCF (bgzipped)")
    parser.add_argument("--wgs_vcf", required=True, help="WGS VCF (bgzipped)")
    parser.add_argument("--gtf", required=True, help="GENCODE GTF file (can be .gtf or .gtf.gz)")
    parser.add_argument("--output", default="universal_hk_snps.tsv", help="Output TSV file")
    parser.add_argument("--tmp_dir", default="/tmp", help="Temporary directory")
    parser.add_argument("--keep_tmp", action="store_true", help="Keep temporary files")
    args = parser.parse_args()

    # Check input files
    for f in [args.rna_vcf, args.wgs_vcf, args.gtf]:
        if not os.path.exists(f):
            print(f"ERROR: File not found: {f}")
            sys.exit(1)

    # Ensure VCFs are indexed
    ensure_vcf_index(args.rna_vcf)
    ensure_vcf_index(args.wgs_vcf)

    # Create temporary directory
    tmp_dir = tempfile.mkdtemp(prefix="hk_panel_", dir=args.tmp_dir)
    print(f"Temporary directory: {tmp_dir}")

    try:
        # Step 1: Get housekeeping gene list
        gene_file = os.path.join(tmp_dir, "hk_genes.txt")
        genes = download_hk_genes(gene_file)

        # Step 2: Extract exons for those genes (using Python)
        exons_bed = os.path.join(tmp_dir, "hk_exons.bed")
        exon_count = extract_exons_python(args.gtf, gene_file, exons_bed)
        if exon_count == 0:
            print("ERROR: No exons found. Check GTF and gene list.")
            return

        # Step 3: Get common chromosomes between the two VCFs (optional but recommended)
        rna_contigs = set(run_cmd(f"bcftools index -s {args.rna_vcf} | cut -f1", "Getting RNA contigs").split())
        wgs_contigs = set(run_cmd(f"bcftools index -s {args.wgs_vcf} | cut -f1", "Getting WGS contigs").split())
        common_contigs = sorted(rna_contigs & wgs_contigs)
        print(f"Common chromosomes: {len(common_contigs)}")

        # Step 4: Extract positions from both VCFs (only on common chromosomes)
        print("Extracting RNA positions...")
        rna_pos = extract_positions(args.rna_vcf, common_contigs)
        print(f"  {len(rna_pos)} RNA positions")

        print("Extracting WGS positions...")
        wgs_pos = extract_positions(args.wgs_vcf, common_contigs)
        print(f"  {len(wgs_pos)} WGS positions")

        # Step 5: Compute intersection
        common_pos = rna_pos & wgs_pos
        print(f"Common positions between RNA and WGS: {len(common_pos)}")

        if len(common_pos) == 0:
            print("ERROR: No common positions found. Cannot build panel.")
            return

        # Step 6: Intersect common positions with exons
        # Convert common_pos to BED
        common_bed = os.path.join(tmp_dir, "common_pos.bed")
        with open(common_bed, 'w') as f:
            for chrom, pos in sorted(common_pos):
                f.write(f"{chrom}\t{pos-1}\t{pos}\n")

        # Run bedtools intersect
        intersect_bed = os.path.join(tmp_dir, "intersect.bed")
        cmd = f"bedtools intersect -a {exons_bed} -b {common_bed} -wa | sort -u > {intersect_bed}"
        run_cmd(cmd, "Intersecting with exons")
        overlap_count = int(run_cmd(f"wc -l < {intersect_bed}", check=False))
        print(f"Positions in exons: {overlap_count}")

        if overlap_count == 0:
            print("ERROR: No common positions fall within housekeeping exons.")
            return

        # Step 7: Get alleles from WGS VCF (could also use RNA, but WGS is more reliable)
        print("Retrieving REF/ALT alleles...")
        # We need to convert intersect_bed back to a set of positions
        intersect_pos = set()
        with open(intersect_bed) as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 3:
                    chrom = parts[0]
                    pos = int(parts[2])  # end coordinate (1-based)
                    intersect_pos.add((chrom, pos))

        final_count = get_alleles(args.wgs_vcf, intersect_pos, args.output)

        # Summary
        print("\n" + "="*60)
        print("Panel built successfully!")
        print(f"Housekeeping genes used: {len(genes)}")
        print(f"Exons extracted: {exon_count}")
        print(f"Common RNA-WGS positions: {len(common_pos)}")
        print(f"Overlapping exons: {overlap_count}")
        print(f"Final panel SNPs: {final_count}")
        print(f"Output file: {args.output}")
        print("="*60)

    finally:
        if not args.keep_tmp:
            import shutil
            shutil.rmtree(tmp_dir)
            print(f"Temporary directory {tmp_dir} removed.")

if __name__ == "__main__":
    main()


# python build_universal_hk_panel.py \
#     --rna_vcf /path/to/RNA.vcf.gz \
#     --wgs_vcf /path/to/WGS.vcf.gz \
#     --gtf /path/to/gencode.v38.annotation.gtf \
#     --output universal_hk_snps.tsv