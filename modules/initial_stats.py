"""
Module: initial_stats
----------------------
Performs initial diagnostic statistics on input VCFs before any processing.
Runs first in the pipeline to help identify issues like chromosome naming mismatches,
build differences, or empty VCFs.
Now includes full position‑level overlap analysis to determine exactly how many
positions are common between RNA and WGS VCFs.
Produces:
  - variant_counts_per_chrom.tsv – per‑chromosome variant counts for RNA and WGS
  - vcf_summary.tsv – total variants, samples, Ti/Tv ratio
  - position_overlap_summary.tsv – detailed overlap statistics (original and stripped names)
  - overlapping_positions_sample.tsv – a sample of up to 1000 overlapping positions
  - chromosome_barplot.png/pdf – bar chart of variant counts per chromosome
  - variant_density_scatter.png/pdf – scatter plot of RNA vs WGS counts per chromosome
  - titv_comparison.png/pdf – bar chart comparing Ti/Tv ratios (if computable)
All outputs are saved in the 'initial_stats' subdirectory.
"""

import os
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from . import utils
import tempfile
import shutil
import sys
import subprocess

def run(rna_vcf, wgs_vcf, output_dir, params):
    """
    Main entry point for initial_stats module.
    
    Args:
        rna_vcf (str): Path to RNA‑seq VCF.
        wgs_vcf (str): Path to WGS VCF.
        output_dir (str): Main pipeline output directory.
        params (dict): Module parameters from config:
            - chromosomes (str or list): Which chromosomes to analyze. "all" or list.
            - test_region (str): Region for quick overlap test (e.g., "chr1:1-1000000").
            - plot_formats (list): List of formats to save plots.
            - max_positions_for_overlap (int): Maximum number of positions to load for full overlap
                                               (if VCFs are huge, set to e.g. 10_000_000 to sample).
                                               Default: 0 (meaning load all).
    """
    logger = utils.get_logger()
    module_out = os.path.join(output_dir, "initial_stats")
    os.makedirs(module_out, exist_ok=True)

    logger.info("Starting initial statistics module...")

    # -------------------------------------------------------------------------
    # 1. Basic VCF information
    # -------------------------------------------------------------------------
    rna_samples = utils.get_vcf_samples(rna_vcf)
    wgs_samples = utils.get_vcf_samples(wgs_vcf)
    logger.info(f"  RNA VCF contains {len(rna_samples)} samples.")
    logger.info(f"  WGS VCF contains {len(wgs_samples)} samples.")

    # Ensure VCFs are indexed for later operations
    utils.ensure_vcf_index(rna_vcf)
    utils.ensure_vcf_index(wgs_vcf)

    # -------------------------------------------------------------------------
    # 2. Get chromosome lists and determine which to analyze
    # -------------------------------------------------------------------------
    rna_contigs_all = set(utils.get_vcf_contigs(rna_vcf))
    wgs_contigs_all = set(utils.get_vcf_contigs(wgs_vcf))

    chrom_param = params.get('chromosomes', 'all')
    if chrom_param == 'all':
        # Default set: chr1-22, chrX, chrY (with and without 'chr')
        desired_chroms = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
        desired_chroms += [str(i) for i in range(1, 23)] + ["X", "Y"]
    elif isinstance(chrom_param, list):
        desired_chroms = chrom_param
    else:
        logger.error("Invalid value for 'chromosomes' parameter. Must be 'all' or a list.")
        raise ValueError("Invalid chromosomes parameter.")

    # Intersect with actual contigs
    rna_contigs = sorted(rna_contigs_all & set(desired_chroms))
    wgs_contigs = sorted(wgs_contigs_all & set(desired_chroms))
    common_contigs = sorted(set(rna_contigs) & set(wgs_contigs))

    logger.info(f"  RNA VCF has {len(rna_contigs_all)} total contigs, after filtering to desired set: {len(rna_contigs)}")
    logger.info(f"  WGS VCF has {len(wgs_contigs_all)} total contigs, after filtering to desired set: {len(wgs_contigs)}")
    logger.info(f"  Common contigs in desired set: {len(common_contigs)}")

    # -------------------------------------------------------------------------
    # 3. Count variants per chromosome
    # -------------------------------------------------------------------------
    def get_variant_counts(vcf, contigs):
        counts = {}
        for chrom in contigs:
            cmd = f"bcftools view -H -r {chrom} {vcf} | wc -l"
            try:
                result = utils.run_cmd(cmd, check=False)
                counts[chrom] = int(result.strip())
            except:
                counts[chrom] = 0
        return counts

    rna_counts = get_variant_counts(rna_vcf, rna_contigs)
    wgs_counts = get_variant_counts(wgs_vcf, wgs_contigs)

    # Save per‑chromosome counts (only common contigs for direct comparison)
    common_df = pd.DataFrame({
        'chromosome': common_contigs,
        'RNA': [rna_counts.get(c, 0) for c in common_contigs],
        'WGS': [wgs_counts.get(c, 0) for c in common_contigs]
    })
    common_df.to_csv(os.path.join(module_out, "variant_counts_per_chrom.tsv"), sep='\t', index=False)

    # Also save all contigs (for reference)
    all_rna_df = pd.DataFrame(list(rna_counts.items()), columns=['chromosome', 'RNA'])
    all_wgs_df = pd.DataFrame(list(wgs_counts.items()), columns=['chromosome', 'WGS'])
    all_df = pd.merge(all_rna_df, all_wgs_df, on='chromosome', how='outer').fillna(0)
    all_df.to_csv(os.path.join(module_out, "all_contig_counts.tsv"), sep='\t', index=False)

    total_rna = sum(rna_counts.values())
    total_wgs = sum(wgs_counts.values())
    logger.info(f"  Total variants in RNA VCF (on all contigs): {total_rna}")
    logger.info(f"  Total variants in WGS VCF (on all contigs): {total_wgs}")

    # -------------------------------------------------------------------------
    # 4. Full position‑level overlap analysis
    # -------------------------------------------------------------------------
    logger.info("  Performing full position‑level overlap analysis...")
    max_positions = params.get('max_positions_for_overlap', 0)  # 0 means load all

    def load_positions(vcf, contigs, max_pos=None):
        """Load all (chrom, pos) from VCF, optionally limiting to max_pos."""
        positions = set()
        total = 0
        for chrom in contigs:
            if max_pos and total >= max_pos:
                break
            cmd = f"bcftools query -f '%CHROM\t%POS\n' -r {chrom} {vcf}"
            # We'll read line by line to avoid memory blow for huge VCFs
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
        logger.debug(f"Loaded {len(positions)} positions from {vcf}")
        return positions

    # Load positions from both VCFs (only on common contigs to save time)
    logger.info(f"    Loading positions from RNA VCF (on {len(common_contigs)} common contigs)...")
    rna_positions = load_positions(rna_vcf, common_contigs, max_positions)
    logger.info(f"      Loaded {len(rna_positions)} RNA positions.")

    logger.info(f"    Loading positions from WGS VCF (on {len(common_contigs)} common contigs)...")
    wgs_positions = load_positions(wgs_vcf, common_contigs, max_positions)
    logger.info(f"      Loaded {len(wgs_positions)} WGS positions.")

    # Intersection with original names
    intersect_original = rna_positions & wgs_positions
    logger.info(f"      Overlap (original names): {len(intersect_original)} positions.")

    # Check for naming mismatch: if RNA has "chr" and WGS doesn't, try stripping
    if len(intersect_original) == 0 and len(rna_positions) > 0 and len(wgs_positions) > 0:
        logger.info("    No overlap with original names. Trying to strip 'chr' from RNA positions...")
        rna_stripped = {(c.replace('chr', ''), p) for c, p in rna_positions if c.startswith('chr')}
        intersect_stripped = rna_stripped & wgs_positions
        logger.info(f"      Overlap after stripping 'chr': {len(intersect_stripped)} positions.")
        if len(intersect_stripped) > 0:
            logger.info("    Chromosome naming mismatch detected! RNA has 'chr' prefix, WGS does not.")
            # Also try adding 'chr' to WGS
            wgs_with_chr = {('chr' + c, p) for c, p in wgs_positions}
            intersect_added = rna_positions & wgs_with_chr
            logger.info(f"      Overlap after adding 'chr' to WGS: {len(intersect_added)} positions.")
    else:
        intersect_stripped = set()
        intersect_added = set()

    # Save overlap statistics
    with open(os.path.join(module_out, 'position_overlap_summary.tsv'), 'w') as f:
        f.write("Metric\tValue\n")
        f.write(f"RNA_positions_loaded\t{len(rna_positions)}\n")
        f.write(f"WGS_positions_loaded\t{len(wgs_positions)}\n")
        f.write(f"Overlap_original_names\t{len(intersect_original)}\n")
        f.write(f"Overlap_after_stripping_chr_from_RNA\t{len(intersect_stripped)}\n")
        f.write(f"Overlap_after_adding_chr_to_WGS\t{len(intersect_added)}\n")
        f.write(f"Fraction_of_RNA_overlapping_original\t{len(intersect_original)/len(rna_positions) if rna_positions else 0:.6f}\n")
        f.write(f"Fraction_of_WGS_overlapping_original\t{len(intersect_original)/len(wgs_positions) if wgs_positions else 0:.6f}\n")

    # Save a sample of overlapping positions (if any)
    if len(intersect_original) > 0:
        sample = list(intersect_original)[:1000]
        with open(os.path.join(module_out, 'overlapping_positions_sample.tsv'), 'w') as f:
            f.write("chrom\tpos\n")
            for chrom, pos in sorted(sample):
                f.write(f"{chrom}\t{pos}\n")
    elif len(intersect_stripped) > 0:
        sample = list(intersect_stripped)[:1000]
        with open(os.path.join(module_out, 'overlapping_positions_sample_stripped.tsv'), 'w') as f:
            f.write("chrom\tpos\n")
            for chrom, pos in sorted(sample):
                f.write(f"{chrom}\t{pos}\n")
    else:
        logger.warning("    No overlapping positions found even after name adjustments.")

    # -------------------------------------------------------------------------
    # 5. Compute Ti/Tv ratio for each VCF (if possible)
    # -------------------------------------------------------------------------
    def compute_titv(vcf, chroms=None, max_sites=1000000):
        """
        Compute Ti/Tv ratio by sampling up to max_sites variants.
        Returns (ti, tv, ratio) or (None, None, None) if not computable.
        """
        if chroms:
            chrom_filter = " -r " + ",".join(chroms[:5])  # limit to first 5 to avoid huge queries
        else:
            chrom_filter = ""
        cmd = f"bcftools view {vcf} {chrom_filter} -H | head -n {max_sites} | cut -f4,5"
        result = utils.run_cmd(cmd, check=False)
        if not result:
            return None, None, None
        ti = 0
        tv = 0
        for line in result.split('\n'):
            if not line:
                continue
            ref, alt = line.split('\t')
            alt = alt.split(',')[0]  # take first alternate
            if len(ref) == 1 and len(alt) == 1:
                if {ref, alt} in [{'A','G'}, {'C','T'}]:
                    ti += 1
                else:
                    tv += 1
        if tv == 0:
            return ti, tv, None
        return ti, tv, ti/tv

    rna_ti, rna_tv, rna_titv = compute_titv(rna_vcf, common_contigs)
    wgs_ti, wgs_tv, wgs_titv = compute_titv(wgs_vcf, common_contigs)

    # -------------------------------------------------------------------------
    # 6. Quick overlap test on a small region (original, for consistency)
    # -------------------------------------------------------------------------
    test_region = params.get('test_region', 'chr1:1-1000000')
    logger.info(f"  Performing quick overlap test on {test_region}...")
    tmp_dir = tempfile.mkdtemp(dir=module_out, prefix="test_isec_")
    cmd = f"bcftools isec -r {test_region} -c all -n=2 {rna_vcf} {wgs_vcf} -p {tmp_dir}"
    try:
        utils.run_cmd(cmd, "Quick overlap test", check=True)
        intersect_file = os.path.join(tmp_dir, "0002.vcf")
        if os.path.exists(intersect_file) and os.path.getsize(intersect_file) > 0:
            overlap_count = int(utils.run_cmd(f"grep -v '^#' {intersect_file} | wc -l"))
            logger.info(f"    Found {overlap_count} common sites in {test_region}.")
        else:
            overlap_count = 0
            logger.info("    No common sites found in this region.")
    except Exception as e:
        logger.error(f"    Quick overlap test failed: {e}")
        overlap_count = -1
    finally:
        shutil.rmtree(tmp_dir, ignore_errors=True)

    # -------------------------------------------------------------------------
    # 7. Write summary file
    # -------------------------------------------------------------------------
    with open(os.path.join(module_out, 'vcf_summary.tsv'), 'w') as f:
        f.write("Metric\tRNA\tWGS\n")
        f.write(f"Total_variants\t{total_rna}\t{total_wgs}\n")
        f.write(f"Number_of_samples\t{len(rna_samples)}\t{len(wgs_samples)}\n")
        f.write(f"Ti_count\t{rna_ti if rna_ti is not None else 'NA'}\t{wgs_ti if wgs_ti is not None else 'NA'}\n")
        f.write(f"Tv_count\t{rna_tv if rna_tv is not None else 'NA'}\t{wgs_tv if wgs_tv is not None else 'NA'}\n")
        f.write(f"Ti/Tv_ratio\t{rna_titv if rna_titv is not None else 'NA'}\t{wgs_titv if wgs_titv is not None else 'NA'}\n")

    # -------------------------------------------------------------------------
    # 8. Generate plots
    # -------------------------------------------------------------------------
    utils.set_publication_style()
    plot_formats = params.get('plot_formats', ['png', 'pdf'])

    # Bar plot of variant counts per chromosome (top 24)
    if not common_df.empty:
        top_chroms = common_df.head(24)
        plt.figure(figsize=(12,6))
        x = np.arange(len(top_chroms))
        width = 0.35
        plt.bar(x - width/2, top_chroms['RNA'], width, label='RNA', color='steelblue')
        plt.bar(x + width/2, top_chroms['WGS'], width, label='WGS', color='darkorange')
        plt.xlabel('Chromosome')
        plt.ylabel('Number of variants')
        plt.title('Variant counts per chromosome (common contigs)')
        plt.xticks(x, top_chroms['chromosome'], rotation=45)
        plt.legend()
        plt.tight_layout()
        for fmt in plot_formats:
            plt.savefig(os.path.join(module_out, f'chromosome_barplot.{fmt}'), dpi=300)
        plt.close()
    else:
        logger.warning("No common chromosomes to plot bar chart.")

    # Scatter plot of RNA vs WGS counts per chromosome
    if not common_df.empty:
        plt.figure(figsize=(8,8))
        plt.scatter(common_df['RNA'], common_df['WGS'], alpha=0.7, s=50, color='steelblue')
        max_val = max(common_df['RNA'].max(), common_df['WGS'].max())
        plt.plot([0, max_val], [0, max_val], 'r--', alpha=0.5, label='y=x')
        plt.xlabel('RNA variant count')
        plt.ylabel('WGS variant count')
        plt.title('Per‑chromosome variant counts (common contigs)')
        for _, row in common_df.iterrows():
            plt.annotate(row['chromosome'], (row['RNA'], row['WGS']), fontsize=8, alpha=0.7)
        plt.legend()
        plt.tight_layout()
        for fmt in plot_formats:
            plt.savefig(os.path.join(module_out, f'variant_density_scatter.{fmt}'), dpi=300)
        plt.close()
    else:
        logger.warning("No common chromosomes to plot scatter.")

    # Ti/Tv comparison bar plot (if available)
    if rna_titv is not None and wgs_titv is not None:
        plt.figure(figsize=(5,6))
        plt.bar(['RNA', 'WGS'], [rna_titv, wgs_titv], color=['steelblue', 'darkorange'])
        plt.ylabel('Ti/Tv ratio')
        plt.title('Ti/Tv comparison')
        plt.axhline(y=2.0, color='gray', linestyle='--', label='Expected DNA Ti/Tv ~2.0')
        plt.axhline(y=2.8, color='lightgray', linestyle='--', label='Expected RNA Ti/Tv ~2.5‑3.0')
        plt.legend()
        plt.tight_layout()
        for fmt in plot_formats:
            plt.savefig(os.path.join(module_out, f'titv_comparison.{fmt}'), dpi=300)
        plt.close()
    else:
        logger.warning("Ti/Tv could not be computed for one or both VCFs; skipping plot.")

    logger.info(f"Initial statistics completed. Results in {module_out}")