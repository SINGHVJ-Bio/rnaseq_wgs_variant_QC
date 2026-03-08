"""
Module 1: variant_qc
----------------------
Compares RNA‑seq VCF against VerifyBAMID2 SNP panel.
Now also optionally performs QC on a WGS VCF and generates extensive QC plots.
Produces:
  - verifybamid_snps.tsv – list of resource SNPs
  - rna_filtered.vcf.gz – filtered RNA VCF (if apply_filters true)
  - rna_metrics.tsv – RNA‑specific INFO fields
  - rna_verifybamid_sites.vcf.gz – intersection with resource SNPs
  - mapping_stats.txt – counts of resource SNPs present/missing
  - snp_details.tsv – per‑SNP details
  - ti_tv_stats.txt – Ti/Tv for intersecting SNPs
  - rna_plots/ – extensive QC plots for RNA (publication‑quality)
  - wgs_plots/ – extensive QC plots for WGS (if run_wgs_qc true)
  - comparison_plots/ – comparison plots if both VCFs provided
"""

import os
import subprocess
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from . import utils

# -----------------------------------------------------------------------------
# Main entry point
# -----------------------------------------------------------------------------
def run(resource_bed, rna_vcf, wgs_vcf, output_dir, stats_path, params):
    """
    Main entry point for variant_qc module.
    
    Args:
        resource_bed (str): Path to VerifyBAMID2 resource BED file (5 columns).
        rna_vcf (str): Path to RNA‑seq joint‑call VCF (bgzipped).
        wgs_vcf (str): Optional path to WGS VCF (from top‑level config).
        output_dir (str): Main pipeline output directory.
        stats_path (str): Path to pipeline statistics JSON file.
        params (dict): Module parameters from config:
            - min_qual (int): Minimum QUAL score (used only if apply_filters True).
            - max_missing (float): Max missing fraction (used only if apply_filters True).
            - strip_chr (bool): If True, remove "chr" prefix from resource chromosomes.
            - apply_filters (bool): If False, skip QUAL/missingness filtering.
            - run_wgs_qc (bool): If True and wgs_vcf provided, run WGS QC.
            - extra_plots (bool): If True, generate additional QC plots (default True).
    """
    logger = utils.get_logger()
    module_out = os.path.join(output_dir, "variant_qc")
    os.makedirs(module_out, exist_ok=True)

    # Initialize statistics for this module
    stats = {
        'input_rna_vcf': rna_vcf,
        'resource_snps_total': 0,
        'filtered_vcf_variants': 0,
        'resource_snps_present': 0,
        'resource_snps_missing': 0,
        'fraction_present': 0.0,
        'ti_count': 0,
        'tv_count': 0,
        'titv_ratio': 0.0
    }

    # -------------------------------------------------------------------------
    # Step 1: Convert resource BED to list of 1‑based positions with ref/alt
    # -------------------------------------------------------------------------
    snp_tsv = os.path.join(module_out, "verifybamid_snps.tsv")
    _convert_bed_to_snp_list(resource_bed, snp_tsv, params.get('strip_chr', False))
    
    # Get total resource SNPs for stats
    df_snps = pd.read_csv(snp_tsv, sep='\t')
    stats['resource_snps_total'] = len(df_snps)

    # -------------------------------------------------------------------------
    # Step 2: Optionally filter RNA VCF by QUAL and missingness
    # -------------------------------------------------------------------------
    apply_filters = params.get('apply_filters', True)
    if apply_filters:
        filtered_vcf = os.path.join(module_out, "rna_filtered.vcf.gz")
        _filter_vcf(rna_vcf, filtered_vcf, params)
        logger.info(f"  Filtered VCF saved to {filtered_vcf}")
        # Count variants in filtered VCF
        cmd = f"bcftools view -H {filtered_vcf} | wc -l"
        result = utils.run_cmd(cmd, "Counting variants in filtered VCF")
        stats['filtered_vcf_variants'] = int(result.strip())
        rna_vcf_to_use = filtered_vcf
    else:
        logger.info("  apply_filters is false; using input RNA VCF directly.")
        rna_vcf_to_use = rna_vcf
        utils.ensure_vcf_index(rna_vcf_to_use)
        cmd = f"bcftools view -H {rna_vcf_to_use} | wc -l"
        result = utils.run_cmd(cmd, "Counting variants in input VCF")
        stats['filtered_vcf_variants'] = int(result.strip())

    # -------------------------------------------------------------------------
    # Step 3: Extract RNA‑specific metrics (only fields present in VCF header)
    # -------------------------------------------------------------------------
    metrics_tsv = os.path.join(module_out, "rna_metrics.tsv")
    _extract_rna_metrics(rna_vcf_to_use, metrics_tsv)

    # -------------------------------------------------------------------------
    # Step 4: Intersect filtered RNA VCF with resource SNP list
    # -------------------------------------------------------------------------
    targets_file = os.path.join(module_out, "targets.txt")
    with open(targets_file, 'w') as f:
        df_snps = pd.read_csv(snp_tsv, sep='\t')
        for _, row in df_snps.iterrows():
            f.write(f"{row['chrom']}\t{row['pos']}\n")

    subset_vcf = os.path.join(module_out, "rna_verifybamid_sites.vcf.gz")
    cmd = f"bcftools view -T {targets_file} {rna_vcf_to_use} -Oz -o {subset_vcf}"
    utils.run_cmd(cmd, "Intersecting with resource SNPs")
    utils.run_cmd(f"bcftools index {subset_vcf}", "Indexing subset VCF")

    # -------------------------------------------------------------------------
    # Step 5: Generate mapping statistics (counts present/missing)
    # -------------------------------------------------------------------------
    stats_file = os.path.join(module_out, "mapping_stats.txt")
    present_count, missing_count = _generate_stats(snp_tsv, rna_vcf_to_use, subset_vcf, stats_file)
    stats['resource_snps_present'] = present_count
    stats['resource_snps_missing'] = missing_count
    stats['fraction_present'] = present_count / stats['resource_snps_total'] if stats['resource_snps_total'] > 0 else 0.0

    # -------------------------------------------------------------------------
    # Step 6: Create detailed per‑SNP table and Ti/Tv stats
    # -------------------------------------------------------------------------
    detail_file = os.path.join(module_out, "snp_details.tsv")
    ti_tv_file = os.path.join(module_out, "ti_tv_stats.txt")
    ti_count, tv_count = _create_detail_table(snp_tsv, subset_vcf, detail_file, ti_tv_file)
    stats['ti_count'] = ti_count
    stats['tv_count'] = tv_count
    if tv_count > 0:
        stats['titv_ratio'] = ti_count / tv_count

    # Update global statistics
    utils.update_stats(stats_path, 'variant_qc', stats)

    # -------------------------------------------------------------------------
    # Step 7: Generate QC plots for RNA VCF (basic and extra)
    # -------------------------------------------------------------------------
    extra_plots = params.get('extra_plots', True)
    _generate_vcf_plots(rna_vcf_to_use, module_out, "rna", params, extra_plots)

    # -------------------------------------------------------------------------
    # Step 8: If WGS VCF provided and run_wgs_qc true, generate WGS QC plots
    # -------------------------------------------------------------------------
    run_wgs_qc = params.get('run_wgs_qc', False)
    if wgs_vcf and os.path.exists(wgs_vcf) and run_wgs_qc:
        logger.info("  Generating WGS QC plots...")
        _generate_vcf_plots(wgs_vcf, module_out, "wgs", params, extra_plots)

        # ---------------------------------------------------------------------
        # Step 9: Generate comparison plots if both VCFs available
        # ---------------------------------------------------------------------
        logger.info("  Generating RNA vs WGS comparison plots...")
        _generate_comparison_plots(rna_vcf_to_use, wgs_vcf, module_out, params)

    logger.info(f"variant_qc completed. Results in: {module_out}")
    return stats

# -----------------------------------------------------------------------------
# Helper functions (existing, unchanged)
# -----------------------------------------------------------------------------
def _convert_bed_to_snp_list(bed_file, out_tsv, strip_chr):
    df = pd.read_csv(bed_file, sep='\t', header=None,
                     names=['chrom', 'start0', 'end1', 'ref', 'alt'])
    df['pos'] = df['start0'] + 1
    if strip_chr:
        df['chrom'] = df['chrom'].str.replace('^chr', '', regex=True)
    df[['chrom', 'pos', 'ref', 'alt']].to_csv(out_tsv, sep='\t', index=False)
    logger = utils.get_logger()
    logger.info(f"  Converted {len(df)} SNPs to {out_tsv}")

def _filter_vcf(in_vcf, out_vcf, params):
    logger = utils.get_logger()
    min_qual = params.get('min_qual', 0)
    max_missing = params.get('max_missing', 1.0)
    filter_expr = f'QUAL>={min_qual} && F_MISSING<={max_missing}'
    cmd = f"bcftools view -i '{filter_expr}' {in_vcf} -Oz -o {out_vcf}"
    utils.run_cmd(cmd, "Filtering VCF")
    utils.run_cmd(f"bcftools index {out_vcf}", "Indexing filtered VCF")
    logger.info(f"  Filtered VCF saved to {out_vcf}")

def _extract_rna_metrics(vcf, out_tsv):
    logger = utils.get_logger()
    desired = ['VDB', 'RPB', 'MQB', 'BQB', 'MQ0F', 'DP', 'QD']
    present = utils.get_vcf_info_fields(vcf)
    available = [f for f in desired if f in present]
    if not available:
        logger.warning("  None of the RNA-specific INFO fields are present.")
        with open(out_tsv, 'w') as f:
            f.write("chrom\tpos\t" + "\t".join(desired) + "\n")
        return
    fields = "%CHROM\t%POS\t" + "\t".join([f"%{f}" for f in available])
    cmd = f"bcftools query -f '{fields}\n' {vcf}"
    result = utils.run_cmd(cmd, "Extracting RNA metrics")
    with open(out_tsv, 'w') as f:
        f.write("chrom\tpos\t" + "\t".join(available) + "\n")
        f.write(result)
    logger.info(f"  RNA metrics saved to {out_tsv} (fields: {', '.join(available)})")

def _generate_stats(snp_tsv, filtered_vcf, subset_vcf, stats_file):
    logger = utils.get_logger()
    df_snps = pd.read_csv(snp_tsv, sep='\t')
    total_snps = len(df_snps)
    cmd = f"bcftools query -f '%CHROM\t%POS\n' {subset_vcf}"
    result = utils.run_cmd(cmd, "Querying present SNPs")
    present_positions = set()
    for line in result.split('\n'):
        if line:
            chrom, pos = line.split()
            present_positions.add((chrom, int(pos)))
    present_count = len(present_positions)
    missing_count = total_snps - present_count
    with open(stats_file, 'w') as f:
        f.write(f"Total VerifyBAMID SNPs\t{total_snps}\n")
        f.write(f"SNPs present in filtered RNA VCF\t{present_count}\n")
        f.write(f"SNPs missing from filtered RNA VCF\t{missing_count}\n")
        f.write(f"Fraction present\t{present_count/total_snps:.3f}\n")
    logger.info(f"  Statistics written to {stats_file}")
    return present_count, missing_count

def _create_detail_table(snp_tsv, subset_vcf, detail_file, ti_tv_file):
    logger = utils.get_logger()
    df_snps = pd.read_csv(snp_tsv, sep='\t')
    df_snps['in_vcf'] = 'no'
    df_snps['vcf_filter'] = '.'
    df_snps['allele_match'] = 'NA'
    df_snps['vcf_ref'] = ''
    df_snps['vcf_alt'] = ''
    df_snps['variant_type'] = ''
    cmd = f"bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%FILTER\n' {subset_vcf}"
    result = utils.run_cmd(cmd, "Querying subset VCF")
    subset_data = {}
    for line in result.split('\n'):
        if line:
            parts = line.split('\t')
            if len(parts) >= 5:
                chrom, pos, ref, alt, filt = parts[:5]
                subset_data[(chrom, int(pos))] = (ref, alt, filt)
    def get_variant_type(ref, alt):
        alt = alt.split(',')[0]
        if len(ref) == 1 and len(alt) == 1:
            if {ref, alt} in [{'A','G'}, {'C','T'}]:
                return 'ts'
            else:
                return 'tv'
        else:
            return 'indel'
    ti_count, tv_count = 0, 0
    for idx, row in df_snps.iterrows():
        key = (row['chrom'], row['pos'])
        if key in subset_data:
            df_snps.at[idx, 'in_vcf'] = 'yes'
            vcf_ref, vcf_alt, vcf_filt = subset_data[key]
            df_snps.at[idx, 'vcf_ref'] = vcf_ref
            df_snps.at[idx, 'vcf_alt'] = vcf_alt
            df_snps.at[idx, 'vcf_filter'] = vcf_filt
            df_snps.at[idx, 'allele_match'] = 'yes' if (vcf_ref == row['ref'] and vcf_alt == row['alt']) else 'no'
            vtype = get_variant_type(vcf_ref, vcf_alt)
            df_snps.at[idx, 'variant_type'] = vtype
            if vtype == 'ts':
                ti_count += 1
            elif vtype == 'tv':
                tv_count += 1
    df_snps.to_csv(detail_file, sep='\t', index=False)
    logger.info(f"  Detailed SNP table saved to {detail_file}")
    with open(ti_tv_file, 'w') as f:
        f.write(f"Ti count (among present SNPs)\t{ti_count}\n")
        f.write(f"Tv count (among present SNPs)\t{tv_count}\n")
        if tv_count > 0:
            f.write(f"Ti/Tv ratio\t{ti_count/tv_count:.3f}\n")
        else:
            f.write("Ti/Tv ratio\tundefined (no Tv)\n")
    logger.info(f"  Ti/Tv stats written to {ti_tv_file}")
    return ti_count, tv_count

# -----------------------------------------------------------------------------
# Enhanced plotting functions (all with publication style)
# -----------------------------------------------------------------------------
def _generate_vcf_plots(vcf, out_dir, prefix, params, extra_plots=True):
    """
    Generate extensive QC plots for a single VCF, all with publication style.
    """
    # Apply publication style
    utils.set_publication_style()
    logger = utils.get_logger()
    plot_dir = os.path.join(out_dir, f"{prefix}_plots")
    os.makedirs(plot_dir, exist_ok=True)
    plot_formats = ['png', 'pdf']

    # Basic plots
    _basic_vcf_plots(vcf, plot_dir, prefix, plot_formats)

    if not extra_plots:
        return

    # Extra plots
    logger.info(f"    Generating extra QC plots for {prefix}...")
    _plot_overall_titv(vcf, plot_dir, prefix, plot_formats)
    _plot_variant_types(vcf, plot_dir, prefix, plot_formats)
    _plot_het_hom(vcf, plot_dir, prefix, plot_formats)
    _plot_missingness(vcf, plot_dir, prefix, plot_formats)
    _plot_allele_balance(vcf, plot_dir, prefix, plot_formats)

def _basic_vcf_plots(vcf, plot_dir, prefix, plot_formats):
    """Generate basic QC plots (QUAL, DP, Ti/Tv per chrom, variant counts per chrom)."""
    logger = utils.get_logger()

    # QUAL distribution
    cmd = f"bcftools query -f '%QUAL\n' {vcf} | head -n 1000000"
    result = utils.run_cmd(cmd, check=False)
    if result:
        quals = [float(x) for x in result.split('\n') if x]
        plt.figure(figsize=(8,4))
        plt.hist(quals, bins=50, alpha=0.7, color='steelblue', edgecolor='black')
        plt.xlabel('QUAL')
        plt.ylabel('Frequency')
        plt.title(f'QUAL distribution ({prefix})')
        plt.tight_layout()
        for fmt in plot_formats:
            plt.savefig(os.path.join(plot_dir, f'qual_dist.{fmt}'), dpi=300)
        plt.close()
        logger.info(f"    Saved QUAL distribution for {prefix}")

    # DP distribution (if FORMAT/DP present)
    cmd = f"bcftools query -f '[%DP\n]' {vcf} | head -n 1000000"
    result = utils.run_cmd(cmd, check=False)
    if result:
        dps = [int(x) for x in result.split('\n') if x and x != '.']
        if dps:
            plt.figure(figsize=(8,4))
            plt.hist(dps, bins=50, alpha=0.7, color='steelblue', edgecolor='black')
            plt.xlabel('Depth (DP)')
            plt.ylabel('Frequency')
            plt.title(f'Depth distribution ({prefix})')
            plt.tight_layout()
            for fmt in plot_formats:
                plt.savefig(os.path.join(plot_dir, f'dp_dist.{fmt}'), dpi=300)
            plt.close()
            logger.info(f"    Saved DP distribution for {prefix}")

    # Ti/Tv per chromosome
    chroms = utils.get_vcf_contigs(vcf)[:24]  # limit to main chromosomes
    titv_data = []
    for chrom in chroms:
        cmd = f"bcftools view -H -r {chrom} {vcf} | cut -f4,5 | head -n 50000"
        result = utils.run_cmd(cmd, check=False)
        ti, tv = 0, 0
        for line in result.split('\n'):
            if not line:
                continue
            ref, alt = line.split()
            alt = alt.split(',')[0]
            if len(ref)==1 and len(alt)==1:
                if {ref, alt} in [{'A','G'}, {'C','T'}]:
                    ti += 1
                else:
                    tv += 1
        if tv > 0:
            titv_data.append((chrom, ti/tv))
    if titv_data:
        df = pd.DataFrame(titv_data, columns=['chromosome', 'TiTv'])
        plt.figure(figsize=(12,6))
        sns.barplot(data=df, x='chromosome', y='TiTv', color='steelblue')
        plt.axhline(y=2.0, color='gray', linestyle='--', label='Expected DNA Ti/Tv ~2.0')
        plt.axhline(y=2.8, color='lightgray', linestyle='--', label='Expected RNA Ti/Tv ~2.5‑3.0')
        plt.xlabel('Chromosome')
        plt.ylabel('Ti/Tv ratio')
        plt.title(f'Ti/Tv per chromosome ({prefix})')
        plt.legend()
        plt.tight_layout()
        for fmt in plot_formats:
            plt.savefig(os.path.join(plot_dir, f'titv_per_chrom.{fmt}'), dpi=300)
        plt.close()
        logger.info(f"    Saved Ti/Tv per chromosome for {prefix}")

    # Variant counts per chromosome
    counts = []
    for chrom in chroms:
        cmd = f"bcftools view -H -r {chrom} {vcf} | wc -l"
        result = utils.run_cmd(cmd, check=False)
        if result:
            counts.append(int(result.strip()))
        else:
            counts.append(0)
    df_counts = pd.DataFrame({'chromosome': chroms, 'count': counts})
    plt.figure(figsize=(12,6))
    sns.barplot(data=df_counts, x='chromosome', y='count', color='steelblue')
    plt.xlabel('Chromosome')
    plt.ylabel('Variant count')
    plt.title(f'Variant counts per chromosome ({prefix})')
    plt.xticks(rotation=45)
    plt.tight_layout()
    for fmt in plot_formats:
        plt.savefig(os.path.join(plot_dir, f'variant_counts.{fmt}'), dpi=300)
    plt.close()
    logger.info(f"    Saved variant counts per chromosome for {prefix}")

def _plot_overall_titv(vcf, plot_dir, prefix, plot_formats):
    """Plot overall Ti/Tv ratio as a bar."""
    logger = utils.get_logger()
    cmd = f"bcftools view -H {vcf} | cut -f4,5 | head -n 1000000"
    result = utils.run_cmd(cmd, check=False)
    if not result:
        return
    ti, tv = 0, 0
    for line in result.split('\n'):
        if not line:
            continue
        ref, alt = line.split()
        alt = alt.split(',')[0]
        if len(ref) == 1 and len(alt) == 1:
            if {ref, alt} in [{'A','G'}, {'C','T'}]:
                ti += 1
            else:
                tv += 1
    if tv == 0:
        return
    titv = ti/tv
    plt.figure(figsize=(4,6))
    plt.bar([prefix], [titv], color='steelblue')
    plt.axhline(y=2.0, color='gray', linestyle='--', label='DNA expected ~2.0')
    plt.axhline(y=2.8, color='lightgray', linestyle='--', label='RNA expected ~2.5‑3.0')
    plt.ylabel('Ti/Tv ratio')
    plt.title('Overall Ti/Tv')
    plt.legend()
    plt.tight_layout()
    for fmt in plot_formats:
        plt.savefig(os.path.join(plot_dir, f'overall_titv.{fmt}'), dpi=300)
    plt.close()
    logger.info(f"    Saved overall Ti/Tv for {prefix}")

def _plot_variant_types(vcf, plot_dir, prefix, plot_formats):
    """Plot SNP vs indel counts as a pie chart."""
    logger = utils.get_logger()
    cmd = f"bcftools view -H {vcf} | head -n 1000000 | cut -f4,5"
    result = utils.run_cmd(cmd, check=False)
    if not result:
        return
    snp = 0
    indel = 0
    for line in result.split('\n'):
        if not line:
            continue
        ref, alt = line.split()
        alt = alt.split(',')[0]
        if len(ref) == 1 and len(alt) == 1:
            snp += 1
        else:
            indel += 1
    if snp + indel == 0:
        return
    plt.figure(figsize=(5,5))
    plt.pie([snp, indel], labels=['SNP', 'Indel'], autopct='%1.1f%%',
            colors=['steelblue', 'darkorange'])
    plt.title('Variant type distribution')
    plt.tight_layout()
    for fmt in plot_formats:
        plt.savefig(os.path.join(plot_dir, f'variant_types.{fmt}'), dpi=300)
    plt.close()
    logger.info(f"    Saved variant type distribution for {prefix}")

def _plot_het_hom(vcf, plot_dir, prefix, plot_formats):
    """Plot heterozygous/homozygous ratio per sample (first 50 samples)."""
    logger = utils.get_logger()
    samples = utils.get_vcf_samples(vcf)
    if not samples:
        return
    samples = samples[:50]  # limit for readability
    data = []
    for sample in samples:
        cmd = f"bcftools view -s {sample} -H {vcf} | cut -f10 | head -n 100000"
        result = utils.run_cmd(cmd, check=False)
        if not result:
            continue
        gt_lines = result.split('\n')
        het = 0
        hom = 0
        for gt in gt_lines:
            if not gt:
                continue
            gt_val = gt.split(':')[0]
            if gt_val in ('0/1', '1/0', '0|1', '1|0'):
                het += 1
            elif gt_val in ('0/0', '1/1', '0|0', '1|1'):
                hom += 1
        if het + hom > 0:
            ratio = het / (het + hom) if hom > 0 else 1.0
            data.append((sample, ratio))
    if not data:
        return
    df = pd.DataFrame(data, columns=['sample', 'het_ratio'])
    plt.figure(figsize=(max(8, 0.3*len(df)), 6))
    sns.barplot(data=df, x='sample', y='het_ratio', color='steelblue')
    plt.xticks(rotation=90)
    plt.ylabel('Heterozygous / (Heterozygous+Homozygous)')
    plt.title('Heterozygosity ratio per sample')
    plt.tight_layout()
    for fmt in plot_formats:
        plt.savefig(os.path.join(plot_dir, f'het_hom_ratio.{fmt}'), dpi=300)
    plt.close()
    logger.info(f"    Saved heterozygosity ratio per sample for {prefix}")

def _plot_missingness(vcf, plot_dir, prefix, plot_formats):
    """Plot missing genotype rate per sample (first 50 samples)."""
    logger = utils.get_logger()
    samples = utils.get_vcf_samples(vcf)
    if not samples:
        return
    samples = samples[:50]
    data = []
    for sample in samples:
        cmd = f"bcftools view -s {sample} -H {vcf} | cut -f10 | head -n 100000"
        result = utils.run_cmd(cmd, check=False)
        if not result:
            continue
        gt_lines = result.split('\n')
        missing = 0
        total = 0
        for gt in gt_lines:
            if not gt:
                continue
            total += 1
            gt_val = gt.split(':')[0]
            if gt_val in ('./.', '.|.'):
                missing += 1
        if total > 0:
            miss_rate = missing / total
            data.append((sample, miss_rate))
    if not data:
        return
    df = pd.DataFrame(data, columns=['sample', 'missing_rate'])
    plt.figure(figsize=(max(8, 0.3*len(df)), 6))
    sns.barplot(data=df, x='sample', y='missing_rate', color='steelblue')
    plt.xticks(rotation=90)
    plt.ylabel('Missing genotype rate')
    plt.title('Missingness per sample')
    plt.tight_layout()
    for fmt in plot_formats:
        plt.savefig(os.path.join(plot_dir, f'missingness.{fmt}'), dpi=300)
    plt.close()
    logger.info(f"    Saved missingness per sample for {prefix}")

def _plot_allele_balance(vcf, plot_dir, prefix, plot_formats):
    """Plot histogram of alternative allele ratio at heterozygous sites (if AD present)."""
    logger = utils.get_logger()
    # Check if AD field exists
    header = utils.run_cmd(f"bcftools view -h {vcf} | grep '##FORMAT=<ID=AD'", check=False)
    if not header:
        logger.debug(f"    No AD field found for {prefix}, skipping allele balance plot.")
        return
    # Sample heterozygous sites and get AD
    cmd = (f"bcftools view -H {vcf} | grep -E '0/1|1/0|0|1|1|0' | head -n 10000 "
           f"| cut -f10 | cut -d':' -f1,2")
    result = utils.run_cmd(cmd, check=False)
    if not result:
        return
    ratios = []
    for line in result.split('\n'):
        if not line:
            continue
        parts = line.split(':')
        if len(parts) < 2:
            continue
        gt, ad = parts[0], parts[1]
        if gt not in ('0/1', '1/0', '0|1', '1|0'):
            continue
        try:
            ref, alt = map(int, ad.split(','))
            if ref + alt > 0:
                ratio = alt / (ref + alt)
                ratios.append(ratio)
        except:
            continue
    if not ratios:
        return
    plt.figure(figsize=(8,4))
    plt.hist(ratios, bins=50, alpha=0.7, color='steelblue', edgecolor='black')
    plt.xlabel('Alternative allele ratio')
    plt.ylabel('Frequency')
    plt.title('Allele balance at heterozygous sites')
    plt.tight_layout()
    for fmt in plot_formats:
        plt.savefig(os.path.join(plot_dir, f'allele_balance.{fmt}'), dpi=300)
    plt.close()
    logger.info(f"    Saved allele balance histogram for {prefix}")

def _generate_comparison_plots(rna_vcf, wgs_vcf, out_dir, params):
    """Generate comparison plots for RNA and WGS VCFs (publication style)."""
    utils.set_publication_style()
    logger = utils.get_logger()
    comp_dir = os.path.join(out_dir, "comparison_plots")
    os.makedirs(comp_dir, exist_ok=True)
    plot_formats = ['png', 'pdf']

    # Get common chromosomes
    rna_contigs = set(utils.get_vcf_contigs(rna_vcf))
    wgs_contigs = set(utils.get_vcf_contigs(wgs_vcf))
    common_chroms = sorted(rna_contigs & wgs_contigs)

    # Get variant counts per chromosome
    rna_counts = {}
    wgs_counts = {}
    for chrom in common_chroms:
        cmd = f"bcftools view -H -r {chrom} {rna_vcf} | wc -l"
        rna_counts[chrom] = int(utils.run_cmd(cmd, check=False) or 0)
        cmd = f"bcftools view -H -r {chrom} {wgs_vcf} | wc -l"
        wgs_counts[chrom] = int(utils.run_cmd(cmd, check=False) or 0)

    # Scatter plot
    df = pd.DataFrame({
        'chromosome': common_chroms,
        'RNA': [rna_counts[c] for c in common_chroms],
        'WGS': [wgs_counts[c] for c in common_chroms]
    })
    plt.figure(figsize=(8,8))
    plt.scatter(df['RNA'], df['WGS'], alpha=0.7, s=50, color='steelblue')
    max_val = max(df['RNA'].max(), df['WGS'].max())
    plt.plot([0, max_val], [0, max_val], 'r--', alpha=0.5, label='y=x')
    plt.xlabel('RNA variant count')
    plt.ylabel('WGS variant count')
    plt.title('Per‑chromosome variant counts: RNA vs WGS')
    for _, row in df.iterrows():
        plt.annotate(row['chromosome'], (row['RNA'], row['WGS']), fontsize=8, alpha=0.7)
    plt.legend()
    plt.tight_layout()
    for fmt in plot_formats:
        plt.savefig(os.path.join(comp_dir, f'variant_counts_scatter.{fmt}'), dpi=300)
    plt.close()
    logger.info("    Saved variant counts scatter plot")

    # Total counts bar plot
    total_rna = sum(rna_counts.values())
    total_wgs = sum(wgs_counts.values())
    plt.figure(figsize=(5,6))
    plt.bar(['RNA', 'WGS'], [total_rna, total_wgs], color=['steelblue', 'darkorange'])
    plt.ylabel('Total variants')
    plt.title('Total variant counts')
    for fmt in plot_formats:
        plt.savefig(os.path.join(comp_dir, f'total_counts.{fmt}'), dpi=300)
    plt.close()
    logger.info("    Saved total counts bar plot")