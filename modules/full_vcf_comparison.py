"""
Module: full_vcf_comparison
----------------------------
Compares the full RNA‑seq VCF with the full WGS VCF at the site level.
Provides a global overview of variant overlap, independent of sample mapping.
Useful for quick QC to ensure VCFs are from the same build and have expected overlap.
User can specify:
  - which RNA VCF to use (original or filtered)
  - which chromosomes to include
  - minimum variants per chromosome for plotting
  - plot formats
Produces:
  - variant_counts.tsv – total variants and sample counts
  - chromosome_comparison.tsv – per‑chromosome variant counts
  - common_sites_summary.tsv – number of positions shared
  - venn_diagram.{png,pdf} – Venn diagram of site overlap
  - chromosome_barplot.{png,pdf} – bar chart of variant counts per chromosome
  - variant_density_scatter.{png,pdf} – scatter plot of RNA vs WGS counts per chromosome
"""

import os
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
import seaborn as sns
from . import utils
import tempfile
import shutil

def run(rna_vcf, wgs_vcf, sample_map, output_dir, params):
    """
    Main entry point for full_vcf_comparison module.
    
    Args:
        rna_vcf (str): Path to RNA‑seq VCF (original or filtered).
        wgs_vcf (str): Path to WGS VCF.
        sample_map (str): Path to sample mapping file (optional).
        output_dir (str): Main pipeline output directory.
        params (dict): Module parameters from config:
            - use_filtered_rna (bool): Whether the RNA VCF is filtered (for info only).
            - chromosomes (str or list): Which chromosomes to analyze. Can be "all" or a list.
            - min_variants_per_chrom (int): Minimum variants to include a chromosome in plots.
            - plot_formats (list): List of formats to save plots, e.g. ["png","pdf"].
    """
    logger = utils.get_logger()
    module_out = os.path.join(output_dir, "full_vcf_comparison")
    os.makedirs(module_out, exist_ok=True)

    logger.info("Starting full VCF comparison module...")

    # -------------------------------------------------------------------------
    # 1. Basic VCF information
    # -------------------------------------------------------------------------
    rna_samples = utils.get_vcf_samples(rna_vcf)
    wgs_samples = utils.get_vcf_samples(wgs_vcf)
    logger.info(f"  RNA VCF contains {len(rna_samples)} samples.")
    logger.info(f"  WGS VCF contains {len(wgs_samples)} samples.")

    if sample_map and os.path.exists(sample_map):
        # We don't need column mapping here; just count rows
        mapping_df = pd.read_csv(sample_map, sep='\t', header=None)
        logger.info(f"  Sample mapping file provided with {len(mapping_df)} rows.")
    else:
        logger.info("  No sample mapping file provided or file not found; skipping sample‑level comparison.")

    # -------------------------------------------------------------------------
    # 2. Get chromosome lists and check naming consistency
    # -------------------------------------------------------------------------
    rna_contigs = utils.get_vcf_contigs(rna_vcf)
    wgs_contigs = utils.get_vcf_contigs(wgs_vcf)
    common_contigs = sorted(set(rna_contigs) & set(wgs_contigs))

    logger.info(f"  RNA VCF has {len(rna_contigs)} contigs, WGS VCF has {len(wgs_contigs)} contigs.")
    logger.info(f"  Common contigs: {len(common_contigs)}")

    # Check chromosome naming
    rna_has_chr = any(c.startswith('chr') for c in rna_contigs)
    wgs_has_chr = any(c.startswith('chr') for c in wgs_contigs)
    if rna_has_chr != wgs_has_chr:
        logger.warning("    Chromosome naming mismatch: RNA has 'chr' prefix? %s, WGS has 'chr' prefix? %s",
                       rna_has_chr, wgs_has_chr)
    else:
        logger.info("    Chromosome naming consistent (both %s 'chr' prefix).",
                    "have" if rna_has_chr else "lack")

    # -------------------------------------------------------------------------
    # 3. Determine which chromosomes to analyze
    # -------------------------------------------------------------------------
    chrom_param = params.get('chromosomes', 'all')
    if chrom_param == 'all':
        chroms_to_analyze = common_contigs
    elif isinstance(chrom_param, list):
        # Keep only those that are in common
        chroms_to_analyze = [c for c in chrom_param if c in common_contigs]
        if not chroms_to_analyze:
            logger.error("None of the specified chromosomes are common between VCFs.")
            raise ValueError("No common chromosomes in specified list.")
    else:
        logger.error("Invalid value for 'chromosomes' parameter. Must be 'all' or a list.")
        raise ValueError("Invalid chromosomes parameter.")

    logger.info(f"  Analyzing {len(chroms_to_analyze)} chromosomes.")

    # -------------------------------------------------------------------------
    # 4. Count variants per chromosome
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

    rna_counts = get_variant_counts(rna_vcf, chroms_to_analyze)
    wgs_counts = get_variant_counts(wgs_vcf, chroms_to_analyze)

    # Save per‑chromosome counts
    chrom_df = pd.DataFrame({
        'chromosome': chroms_to_analyze,
        'RNA': [rna_counts.get(c, 0) for c in chroms_to_analyze],
        'WGS': [wgs_counts.get(c, 0) for c in chroms_to_analyze]
    })
    chrom_df.to_csv(os.path.join(module_out, "chromosome_comparison.tsv"), sep='\t', index=False)

    # Total variant counts
    total_rna = chrom_df['RNA'].sum()
    total_wgs = chrom_df['WGS'].sum()
    logger.info(f"  Total variants in RNA VCF (on analyzed chromosomes): {total_rna}")
    logger.info(f"  Total variants in WGS VCF (on analyzed chromosomes): {total_wgs}")

    # -------------------------------------------------------------------------
    # 5. Compute site overlap (by position, ignoring alleles)
    # -------------------------------------------------------------------------
    logger.info("  Computing site overlap between VCFs...")

    tmp_dir = tempfile.mkdtemp(dir=module_out, prefix="isec_")
    all_intersect = []
    total_union = 0

    for chrom in chroms_to_analyze:
        chrom_dir = os.path.join(tmp_dir, chrom)
        os.makedirs(chrom_dir, exist_ok=True)
        cmd = f"bcftools isec -r {chrom} -c all -n=2 {rna_vcf} {wgs_vcf} -p {chrom_dir}"
        try:
            utils.run_cmd(cmd, f"Intersecting {chrom}", check=True)
            intersect_file = os.path.join(chrom_dir, "0002.vcf")
            if os.path.exists(intersect_file) and os.path.getsize(intersect_file) > 0:
                # Count lines (skip header)
                count = int(utils.run_cmd(f"grep -v '^#' {intersect_file} | wc -l"))
                all_intersect.append(count)
            else:
                all_intersect.append(0)
            # Also get union: number of sites in either VCF on this chromosome
            rna_file = os.path.join(chrom_dir, "0000.vcf")
            wgs_file = os.path.join(chrom_dir, "0001.vcf")
            rna_count = 0
            wgs_count = 0
            if os.path.exists(rna_file) and os.path.getsize(rna_file) > 0:
                rna_count = int(utils.run_cmd(f"grep -v '^#' {rna_file} | wc -l"))
            if os.path.exists(wgs_file) and os.path.getsize(wgs_file) > 0:
                wgs_count = int(utils.run_cmd(f"grep -v '^#' {wgs_file} | wc -l"))
            union = rna_count + wgs_count + all_intersect[-1]  # isec already counted intersection in both
            total_union += union
        except Exception as e:
            logger.error(f"    Intersection failed for {chrom}: {e}")
            all_intersect.append(0)

    total_intersect = sum(all_intersect)

    # Clean up
    shutil.rmtree(tmp_dir, ignore_errors=True)

    logger.info(f"  Total positions shared (by position): {total_intersect}")
    logger.info(f"  Total union of positions (on analyzed chromosomes): {total_union}")

    # -------------------------------------------------------------------------
    # 6. Write summary files
    # -------------------------------------------------------------------------
    with open(os.path.join(module_out, 'variant_counts.tsv'), 'w') as f:
        f.write("VCF\tTotal_variants\tSamples\n")
        f.write(f"RNA\t{total_rna}\t{len(rna_samples)}\n")
        f.write(f"WGS\t{total_wgs}\t{len(wgs_samples)}\n")

    with open(os.path.join(module_out, 'common_sites_summary.tsv'), 'w') as f:
        f.write("Metric\tValue\n")
        f.write(f"Analyzed_chromosomes\t{len(chroms_to_analyze)}\n")
        f.write(f"Total_RNA_variants\t{total_rna}\n")
        f.write(f"Total_WGS_variants\t{total_wgs}\n")
        f.write(f"Shared_positions\t{total_intersect}\n")
        f.write(f"Union_positions\t{total_union}\n")
        f.write(f"Fraction_of_RNA_in_intersection\t{total_intersect/total_rna if total_rna>0 else 0:.4f}\n")
        f.write(f"Fraction_of_WGS_in_intersection\t{total_intersect/total_wgs if total_wgs>0 else 0:.4f}\n")

    # -------------------------------------------------------------------------
    # 7. Generate plots
    # -------------------------------------------------------------------------
    utils.set_publication_style()
    plot_formats = params.get('plot_formats', ['png', 'pdf'])
    min_variants = params.get('min_variants_per_chrom', 0)

    # Filter chromosomes for plotting (optional)
    plot_chrom_df = chrom_df[chrom_df['RNA'] + chrom_df['WGS'] >= min_variants]
    if len(plot_chrom_df) == 0:
        logger.warning("  No chromosomes meet min_variants_per_chrom; using all chromosomes.")
        plot_chrom_df = chrom_df

    # Venn diagram of total sites
    plt.figure(figsize=(6,6))
    venn2(subsets=(total_rna - total_intersect, total_wgs - total_intersect, total_intersect),
          set_labels=('RNA', 'WGS'))
    plt.title(f'Overlap of variant sites\n(total RNA: {total_rna}, total WGS: {total_wgs})')
    for fmt in plot_formats:
        plt.savefig(os.path.join(module_out, f'venn_diagram.{fmt}'), dpi=300, bbox_inches='tight')
    plt.close()

    # Bar plot of variant counts per chromosome (top 24)
    top_chroms = plot_chrom_df.head(24)  # limit to 24 for readability
    plt.figure(figsize=(12,6))
    x = np.arange(len(top_chroms))
    width = 0.35
    plt.bar(x - width/2, top_chroms['RNA'], width, label='RNA', color='steelblue')
    plt.bar(x + width/2, top_chroms['WGS'], width, label='WGS', color='darkorange')
    plt.xlabel('Chromosome')
    plt.ylabel('Number of variants')
    plt.title('Variant counts per chromosome')
    plt.xticks(x, top_chroms['chromosome'], rotation=45)
    plt.legend()
    plt.tight_layout()
    for fmt in plot_formats:
        plt.savefig(os.path.join(module_out, f'chromosome_barplot.{fmt}'), dpi=300)
    plt.close()

    # Scatter plot of RNA vs WGS counts per chromosome
    plt.figure(figsize=(8,8))
    plt.scatter(plot_chrom_df['RNA'], plot_chrom_df['WGS'], alpha=0.7, s=50, color='steelblue')
    max_val = max(plot_chrom_df['RNA'].max(), plot_chrom_df['WGS'].max())
    plt.plot([0, max_val], [0, max_val], 'r--', alpha=0.5, label='y=x')
    plt.xlabel('RNA variant count')
    plt.ylabel('WGS variant count')
    plt.title('Per‑chromosome variant counts: RNA vs WGS')
    for _, row in plot_chrom_df.iterrows():
        plt.annotate(row['chromosome'], (row['RNA'], row['WGS']), fontsize=8, alpha=0.7)
    plt.legend()
    plt.tight_layout()
    for fmt in plot_formats:
        plt.savefig(os.path.join(module_out, f'variant_density_scatter.{fmt}'), dpi=300)
    plt.close()

    logger.info(f"Full VCF comparison completed. Results in {module_out}")