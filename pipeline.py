#!/usr/bin/env python3
"""
RNA‑seq vs WGS Variant QC Pipeline
===================================
Run: python pipeline.py --config config.yaml

This pipeline performs quality control on RNA‑seq variant calls by comparing
them with matched WGS data. It consists of six modules:
  0. initial_stats       – quick diagnostic statistics of input VCFs (run first)
  1. variant_qc          – compare RNA VCF with a reference SNP panel (VerifyBAMID2)
  2. compare_wgs_rnaseq  – detailed comparison between RNA and WGS VCFs
  3. study_snp_analysis  – analyse user‑defined custom SNP panels (e.g., forensic panels)
  4. full_vcf_comparison – quick global comparison of full VCFs (site overlap, variant counts)
  5. exploratory_plotting – generate publication‑quality plots

All results are saved in the specified output directory. A comprehensive
report (output_report.tsv) lists every output file with descriptions and
interpretation guidance. A detailed log file and a statistics summary are also created.
"""

import os
import sys
import subprocess
import yaml
import argparse
from datetime import datetime
from modules import initial_stats, variant_qc, compare_wgs_rnaseq, exploratory_plotting, study_snp_analysis, full_vcf_comparison
from modules.utils import setup_logging, get_logger, init_stats, write_stats_summary

def check_dependencies():
    """Ensure required external tools (bcftools, bedtools) are available."""
    required_tools = ['bcftools', 'bedtools']
    missing = []
    for tool in required_tools:
        try:
            subprocess.run([tool, '--version'], capture_output=True, check=True)
        except (subprocess.CalledProcessError, FileNotFoundError):
            missing.append(tool)
    if missing:
        sys.exit(f"ERROR: Missing required tools: {', '.join(missing)}. Please run setup.sh first.")

def generate_report(output_dir, config):
    """
    Generate a comprehensive tabular report of all output files.
    The report is saved as 'output_report.tsv' in the main output directory.
    Columns: Category, File Name, Full Path, Description, Importance, Interpretation.
    """
    logger = get_logger()
    report_path = os.path.join(output_dir, "output_report.tsv")
    with open(report_path, 'w') as f:
        f.write("Category\tFile Name\tFull Path\tDescription\tImportance\tInterpretation\n")

        # Helper to add a row if the file exists
        def add_row(category, filename, description, importance, interpretation):
            full_path = os.path.join(output_dir, category, filename)
            if os.path.exists(full_path):
                f.write(f"{category}\t{filename}\t{full_path}\t{description}\t{importance}\t{interpretation}\n")

        # ---------------------------------------------------------------------
        # initial_stats outputs
        # ---------------------------------------------------------------------
        cat = "initial_stats"
        add_row(cat, "vcf_summary.tsv",
                "Summary statistics for each VCF: total variants, sample count, Ti/Tv ratio.",
                "Quick overview of input data quality and comparability.",
                "Check that both VCFs have variants and Ti/Tv ratios are reasonable.")
        add_row(cat, "variant_counts_per_chrom.tsv",
                "Variant counts per chromosome for common chromosomes.",
                "Reveals chromosome coverage and potential naming issues.",
                "If counts are zero for many chromosomes, check naming conventions.")
        add_row(cat, "all_contig_counts.tsv",
                "Variant counts for all contigs in each VCF.",
                "Useful for spotting unexpected contigs (e.g., unplaced scaffolds).",
                "Large numbers on non‑standard contigs may indicate alignment issues.")
        add_row(cat, "quick_overlap_report.tsv",
                f"Number of common sites in a small test region (default: chr1:1-1,000,000).",
                "Quick check to see if VCFs share any positions.",
                "If zero, likely a chromosome naming or build mismatch.")
        add_row(cat, "chromosome_barplot.png",
                "Bar chart of variant counts per chromosome for common contigs.",
                "Visual comparison of per‑chromosome variant distribution.",
                "RNA should roughly follow WGS pattern; large deviations may indicate issues.")
        add_row(cat, "chromosome_barplot.pdf",
                "Vector version of bar chart.",
                "Publication‑ready figure.",
                "Same interpretation as PNG.")
        add_row(cat, "variant_density_scatter.png",
                "Scatter plot of RNA vs WGS variant counts per chromosome.",
                "Each point is a chromosome; points should lie near the y=x line.",
                "Systematic deviation suggests bias (e.g., RNA missing many variants on large chromosomes).")
        add_row(cat, "variant_density_scatter.pdf",
                "Vector version of scatter plot.",
                "Publication‑ready figure.",
                "Same interpretation as PNG.")
        add_row(cat, "titv_comparison.png",
                "Bar chart comparing Ti/Tv ratios of RNA and WGS (if computable).",
                "Ti/Tv ratio is a classic QC metric; deviations suggest false positives.",
                "Expected RNA Ti/Tv ~2.5‑3.0, DNA ~2.0‑2.1.")
        add_row(cat, "titv_comparison.pdf",
                "Vector version of Ti/Tv bar chart.",
                "Publication‑ready figure.",
                "Same interpretation as PNG.")

        # ---------------------------------------------------------------------
        # variant_qc outputs
        # ---------------------------------------------------------------------
        cat = "variant_qc"
        add_row(cat, "verifybamid_snps.tsv",
                "List of resource SNPs from VerifyBAMID2 (chrom, pos, ref, alt).",
                "Defines the set of variants used for contamination estimation.",
                "Check that the resource matches your genome build (GRCh38).")
        add_row(cat, "rna_filtered.vcf.gz",
                "RNA‑seq VCF after applying QUAL and missingness filters (if apply_filters true).",
                "Serves as the high‑confidence RNA variant set for downstream comparison.",
                "Number of variants remaining indicates stringency of filters.")
        add_row(cat, "rna_metrics.tsv",
                "RNA‑specific INFO fields (VDB, RPB, MQB, BQB, MQ0F, DP, QD) for all sites.",
                "These metrics help identify technical artefacts (strand bias, read position bias).",
                "Plot distributions to set optimal thresholds for filtering.")
        add_row(cat, "rna_verifybamid_sites.vcf.gz",
                "Intersection of filtered RNA VCF with resource SNP positions.",
                "Shows which resource SNPs are captured in your RNA data.",
                "Low overlap may indicate expression issues or sample mismatch.")
        add_row(cat, "mapping_stats.txt",
                "Counts of resource SNPs present/missing in RNA VCF.",
                "Quick overview of how well the RNA data covers the reference panel.",
                "Fraction present should be high (>80%) for reliable contamination estimates.")
        add_row(cat, "snp_details.tsv",
                "Per‑SNP table with presence, allele match, variant type.",
                "Detailed view of each resource SNP's behaviour in RNA data.",
                "Allele mismatches could indicate sequencing errors or RNA editing.")
        add_row(cat, "ti_tv_stats.txt",
                "Ti/Tv ratio for the intersecting SNPs.",
                "Ti/Tv ratio is a classic QC metric; deviations suggest false positives.",
                "Expected RNA Ti/Tv ~2.5‑3.0 in coding regions.")

        # ---------------------------------------------------------------------
        # compare_wgs_rnaseq outputs
        # ---------------------------------------------------------------------
        cat = "compare_wgs_rnaseq"
        add_row(cat, "filtered_sample_map.tsv",
                "Sample mapping after filtering to those present in both VCFs (and expression).",
                "Ensures all downstream analyses use only common samples.",
                "Verify that all expected sample pairs are present.")
        add_row(cat, "filtered_expression_matrix.csv",
                "Expression matrix subset to common samples (if expression file was provided).",
                "Allows stratification of variant concordance by expression level.",
                "Check that sample IDs match and expression values are reasonable.")
        add_row(cat, "common_sites.vcf.gz",
                "Sites present in both RNA and WGS VCFs.",
                "The basis for all pairwise comparisons.",
                "Number of sites reflects overlap between call sets.")
        add_row(cat, "genotypes.tsv",
                "Genotype table for all sample pairs at common sites (GT, DP, GQ).",
                "Raw data for computing concordance, ASE, etc.",
                "Used as input for all subsequent QC metrics.")
        add_row(cat, "concordance_summary.tsv",
                "Per‑sample precision, recall, concordance (RNA vs WGS).",
                "Quantifies how well RNA‑seq variants match the DNA truth.",
                "High precision (>0.9) indicates few false positives; high recall (>0.8) indicates good sensitivity.")
        add_row(cat, "concordance_by_region.tsv",
                "Same metrics stratified by exonic / intronic regions (if exome BED provided).",
                "Reveals whether variant calling performance differs between regions.",
                "Exonic regions should show higher concordance due to better coverage.")
        add_row(cat, "titv_comparison.tsv",
                "Ti/Tv ratio for RNA and WGS separately.",
                "Compares mutation spectra between the two technologies.",
                "Similar Ti/Tv values suggest consistent variant quality.")
        add_row(cat, "ase_analysis.tsv",
                "Allelic balance at WGS‑heterozygous sites (WGS alt ratio vs RNA alt ratio).",
                "Detects allele‑specific expression or mapping bias.",
                "Deviation from diagonal may indicate biological ASE or technical artefacts.")
        add_row(cat, "depth_correlation.tsv",
                "Pearson correlation of log10(depth) between RNA and WGS.",
                "Measures whether coverage patterns are similar.",
                "High correlation (>0.5) suggests consistent coverage across technologies.")

        # ---------------------------------------------------------------------
        # full_vcf_comparison outputs
        # ---------------------------------------------------------------------
        cat = "full_vcf_comparison"
        add_row(cat, "variant_counts.tsv",
                "Total variant counts and sample numbers for each VCF.",
                "Quick sanity check: RNA should have fewer variants than WGS (only expressed regions).",
                "If RNA has more variants than WGS, it may indicate RNA editing or artefacts.")
        add_row(cat, "chromosome_comparison.tsv",
                "Per‑chromosome variant counts for RNA and WGS.",
                "Reveals whether coverage is uniform across chromosomes.",
                "Check for chromosomes with zero variants (e.g., Y chromosome in RNA).")
        add_row(cat, "common_sites_summary.tsv",
                "Number of positions shared between VCFs (by position only).",
                "Indicates how much of the RNA variant set overlaps with WGS.",
                "Low overlap (<50%) may suggest different genome builds or sample mismatches.")
        add_row(cat, "venn_diagram.png",
                "Venn diagram of total variant sites overlap.",
                "Visual summary of shared and unique variants.",
                "RNA‑specific sites could be RNA editing or false positives.")
        add_row(cat, "venn_diagram.pdf",
                "Vector version of Venn diagram.",
                "Publication‑ready figure.",
                "Same interpretation as PNG.")
        add_row(cat, "chromosome_barplot.png",
                "Bar chart of variant counts per chromosome.",
                "Visual comparison of variant distribution.",
                "RNA should roughly follow WGS pattern; large discrepancies may indicate technical issues.")
        add_row(cat, "chromosome_barplot.pdf",
                "Vector version of bar chart.",
                "Publication‑ready figure.",
                "Same interpretation as PNG.")
        add_row(cat, "variant_density_scatter.png",
                "Scatter plot of RNA vs WGS variant counts per chromosome.",
                "Each point is a chromosome; points should lie near the y=x line.",
                "Systematic deviation suggests bias (e.g., RNA missing many variants on large chromosomes).")
        add_row(cat, "variant_density_scatter.pdf",
                "Vector version of scatter plot.",
                "Publication‑ready figure.",
                "Same interpretation as PNG.")

        # ---------------------------------------------------------------------
        # study_snp_analysis outputs (for each study with run: true)
        # ---------------------------------------------------------------------
        study_configs = config.get('study_snps', [])
        for study in study_configs:
            if not study.get('run', False):
                continue
            study_name = study.get('name', 'unnamed_study')
            cat = f"study_snp_analysis/{study_name}"
            add_row(cat, "mapping_stats.txt",
                    f"Statistics for study '{study_name}': counts of SNPs present/missing in RNA VCF.",
                    "Quick overview of how well this specific panel is captured.",
                    "Low presence may indicate technical issues or sample mismatch.")
            add_row(cat, "snp_details.tsv",
                    f"Per‑SNP details for study '{study_name}': presence, allele match, missing rate.",
                    "Detailed view for troubleshooting or validation.",
                    "Check allele matches – mismatches could indicate sequencing errors.")
            add_row(cat, "study_sites.vcf.gz",
                    f"VCF subset of RNA data at the SNP positions of study '{study_name}'.",
                    "Useful for extracting genotypes for this panel.",
                    "Can be used for further analysis or validation.")

        # ---------------------------------------------------------------------
        # exploratory_plotting outputs (both PDF and PNG)
        # ---------------------------------------------------------------------
        cat = "exploratory_plotting"
        plot_basenames = [
            ("rna_metrics_dist", "Multi‑panel histogram of RNA‑specific metrics (VDB, RPB, MQB, BQB, MQ0F, DP, QD).",
             "Visual overview of quality metric distributions.",
             "Helps identify outliers and set filtering thresholds."),
            ("concordance_bars", "Bar plot of precision, recall, concordance per sample.",
             "Easy comparison of sample performance.",
             "Samples with low values should be investigated."),
            ("precision_recall", "Scatter plot of precision vs recall, each point a sample.",
             "Shows trade‑off between false positives and false negatives.",
             "Ideal samples cluster in top‑right corner."),
            ("titv_comparison", "Bar chart comparing RNA and WGS Ti/Tv ratios.",
             "Visual check of spectral consistency.",
             "Large discrepancies may indicate technical problems."),
            ("pca", "PCA plot of genotypes (red=RNA, blue=WGS).",
             "Confirms sample identity: RNA and WGS from same individual should cluster together.",
             "Mismatched clusters indicate sample swaps."),
            ("ibs_heatmap", "Heatmap of pairwise identity‑by‑state (IBS) similarity.",
             "Reveals sample relatedness and confirms expected relationships.",
             "Diagonal should be 1; off‑diagonal values reflect genetic similarity."),
            ("ase_scatter", "Scatter plot of WGS alt ratio vs RNA alt ratio at heterozygous sites.",
             "Visualises allelic balance.",
             "Points off the diagonal indicate ASE or bias."),
            ("ase_deviation", "Histogram of (RNA ratio – WGS ratio).",
             "Summarises overall allelic imbalance.",
             "Distribution centered near zero suggests little global bias.")
        ]

        for base, desc, imp, interp in plot_basenames:
            for ext in ['pdf', 'png']:
                filename = f"{base}.{ext}"
                add_row(cat, filename, desc, imp, interp)

    logger.info(f"Comprehensive output report written to: {report_path}")
    return report_path

def main():
    parser = argparse.ArgumentParser(description="Run variant QC pipeline")
    parser.add_argument("--config", required=True, help="Path to YAML config file")
    args = parser.parse_args()

    # Load configuration
    with open(args.config, 'r') as f:
        config = yaml.safe_load(f)

    # Create main output directory
    os.makedirs(config['output_dir'], exist_ok=True)

    # Set up logging
    log_path = config.get('log_file', 'pipeline.log')
    if not os.path.isabs(log_path):
        log_path = os.path.join(config['output_dir'], log_path)
    logger = setup_logging(log_path)
    logger.info("="*60)
    logger.info("RNA‑seq vs WGS Variant QC Pipeline started")
    logger.info("="*60)

    # Initialize pipeline statistics
    stats_path = init_stats(config['output_dir'])

    # Check dependencies
    check_dependencies()

    # Determine which modules to run
    modules_config = config.get('run_modules', {})
    run_initial = modules_config.get('initial_stats', False)
    run_variant = modules_config.get('variant_qc', False)
    run_compare = modules_config.get('compare_wgs_rnaseq', False)
    run_study = modules_config.get('study_snp_analysis', False)
    run_full = modules_config.get('full_vcf_comparison', False)
    run_plot = modules_config.get('exploratory_plotting', False)

    # -------------------------------------------------------------------------
    # Module 0: initial_stats – run first to diagnose input VCFs
    # -------------------------------------------------------------------------
    if run_initial:
        logger.info("Starting initial_stats module...")
        initial_stats.run(
            rna_vcf=config['rna_joint_vcf'],
            wgs_vcf=config['wgs_joint_vcf'],
            output_dir=config['output_dir'],
            params=config.get('initial_stats', {})
        )
        logger.info("initial_stats module completed.")
    else:
        logger.info("initial_stats module skipped (run: false)")

    # -------------------------------------------------------------------------
    # Module 1: variant_qc
    # -------------------------------------------------------------------------
    if run_variant:
        logger.info("Starting variant_qc module...")
        variant_qc.run(
            resource_bed=config['resource_bed'],
            rna_vcf=config['rna_joint_vcf'],
            wgs_vcf=config.get('wgs_joint_vcf'),   # pass WGS VCF if defined
            output_dir=config['output_dir'],
            stats_path=stats_path,
            params=config.get('variant_qc', {})
        )
        logger.info("variant_qc module completed.")
    else:
        logger.info("variant_qc module skipped (run: false)")

    # -------------------------------------------------------------------------
    # Module 2: compare_wgs_rnaseq
    # -------------------------------------------------------------------------
    if run_compare:
        logger.info("Starting compare_wgs_rnaseq module...")
        compare_wgs_rnaseq.run(
            rna_vcf=config['rna_joint_vcf'],
            wgs_vcf=config['wgs_joint_vcf'],
            sample_map=config['sample_map'],
            gene_expression=config.get('gene_expression', None),
            exome_bed=config.get('exome_bed', None),
            output_dir=config['output_dir'],
            stats_path=stats_path,
            params=config.get('compare_wgs_rnaseq', {})
        )
        logger.info("compare_wgs_rnaseq module completed.")
    else:
        logger.info("compare_wgs_rnaseq module skipped (run: false)")

    # -------------------------------------------------------------------------
    # Module 3: study_snp_analysis
    # -------------------------------------------------------------------------
    if run_study:
        logger.info("Starting study_snp_analysis module...")
        # Determine which RNA VCF to use: filtered if available and desired, otherwise original
        use_filtered = config.get('study_snp_analysis', {}).get('use_filtered_rna', False)
        if use_filtered:
            rna_vcf_for_study = os.path.join(config['output_dir'], 'variant_qc', 'rna_filtered.vcf.gz')
            if not os.path.exists(rna_vcf_for_study):
                logger.error("Filtered RNA VCF not found. Please run 'variant_qc' module first or set use_filtered_rna: false.")
                logger.error("Set 'variant_qc: true' in your config and rerun.")
                rna_vcf_for_study = None
        else:
            rna_vcf_for_study = config['rna_joint_vcf']

        if rna_vcf_for_study and os.path.exists(rna_vcf_for_study):
            study_snp_analysis.run(
                rna_filtered_vcf=rna_vcf_for_study,
                output_dir=config['output_dir'],
                study_configs=config.get('study_snps', []),
                stats_path=stats_path,
                params=config.get('variant_qc', {})  # reuse strip_chr parameter
            )
        else:
            logger.error("RNA VCF not available for study_snp_analysis. Check path and use_filtered_rna setting.")
        logger.info("study_snp_analysis module completed.")
    else:
        logger.info("study_snp_analysis module skipped (run: false)")

    # -------------------------------------------------------------------------
    # Module 4: full_vcf_comparison
    # -------------------------------------------------------------------------
    if run_full:
        logger.info("Starting full_vcf_comparison module...")
        # Determine which RNA VCF to use
        use_filtered = config.get('full_vcf_comparison', {}).get('use_filtered_rna', False)
        if use_filtered:
            rna_vcf_for_full = os.path.join(config['output_dir'], 'variant_qc', 'rna_filtered.vcf.gz')
            if not os.path.exists(rna_vcf_for_full):
                logger.error("Filtered RNA VCF not found. Please run 'variant_qc' module first or set use_filtered_rna: false.")
                logger.error("Using original RNA VCF instead.")
                rna_vcf_for_full = config['rna_joint_vcf']
        else:
            rna_vcf_for_full = config['rna_joint_vcf']

        full_vcf_comparison.run(
            rna_vcf=rna_vcf_for_full,
            wgs_vcf=config['wgs_joint_vcf'],
            sample_map=config.get('sample_map', None),
            output_dir=config['output_dir'],
            params=config.get('full_vcf_comparison', {})
        )
        logger.info("full_vcf_comparison module completed.")
    else:
        logger.info("full_vcf_comparison module skipped (run: false)")

    # -------------------------------------------------------------------------
    # Module 5: exploratory_plotting
    # -------------------------------------------------------------------------
    if run_plot:
        logger.info("Starting exploratory_plotting module...")
        exploratory_plotting.run(
            variant_qc_dir=os.path.join(config['output_dir'], 'variant_qc'),
            compare_dir=os.path.join(config['output_dir'], 'compare_wgs_rnaseq'),
            output_dir=config['output_dir'],
            params=config.get('exploratory_plotting', {})
        )
        logger.info("exploratory_plotting module completed.")
    else:
        logger.info("exploratory_plotting module skipped (run: false)")

    # Generate final report (always, regardless of which modules ran)
    logger.info("Generating comprehensive output report...")
    report_path = generate_report(config['output_dir'], config)
    logger.info(f"Report saved: {report_path}")

    # Write final stats summary
    write_stats_summary(stats_path, config['output_dir'])

    logger.info("Pipeline finished successfully.")
    logger.info("="*60)

if __name__ == "__main__":
    main()