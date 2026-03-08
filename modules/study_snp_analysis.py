"""
Module: study_snp_analysis
----------------------------
Analyzes user-defined SNP panels (from specific studies) in the RNA-seq VCF.
For each study with 'run: true', it:
  - Loads the SNP list (positions or rsIDs).
  - Converts rsIDs to positions if needed (using a provided dbSNP mapping file).
  - Intersects with the filtered RNA VCF.
  - Generates a summary table and statistics.
  - Saves results in a study-named subfolder.
All key metrics are recorded in the pipeline statistics.
"""

import os
import pandas as pd
import numpy as np
from . import utils

def run(rna_filtered_vcf, output_dir, study_configs, stats_path, params):
    """
    Main entry point for study_snp_analysis module.
    
    Args:
        rna_filtered_vcf (str): Path to the filtered RNA VCF (from variant_qc module).
        output_dir (str): Main pipeline output directory.
        study_configs (list): List of study dictionaries from config.
        stats_path (str): Path to pipeline statistics JSON file.
        params (dict): Module parameters (e.g., strip_chr).
    """
    logger = utils.get_logger()
    module_out = os.path.join(output_dir, "study_snp_analysis")
    os.makedirs(module_out, exist_ok=True)

    # Initialize statistics for this module
    stats = {}

    if not study_configs:
        logger.info("  No study SNP panels configured. Skipping module.")
        return

    # Check if filtered RNA VCF exists
    if not os.path.exists(rna_filtered_vcf):
        logger.error(f"Filtered RNA VCF not found: {rna_filtered_vcf}")
        logger.error("Please run the 'variant_qc' module first to generate this file.")
        logger.error("Set 'variant_qc: true' in your config and rerun the pipeline.")
        return

    # Optional: load dbSNP mapping once if any study uses rsID
    dbsnp_map = None
    for study in study_configs:
        if study.get('id_type') == 'rsid' and study.get('run', False):
            dbsnp_file = study.get('dbsnp_mapping')
            if dbsnp_file and os.path.exists(dbsnp_file):
                logger.info(f"  Loading dbSNP mapping from {dbsnp_file}")
                dbsnp_map = utils.load_dbsnp_mapping(dbsnp_file, study.get('build', 'GRCh38'))
                break
            else:
                logger.error(f"  rsID conversion requested but no valid dbSNP mapping file provided for study {study.get('name')}.")
                logger.error("  Please set 'dbsnp_mapping' to a file with columns: rsid, chrom, pos, ref, alt")
                study['run'] = False  # mark as not runnable

    for study in study_configs:
        if not study.get('run', False):
            logger.info(f"  Skipping study panel: {study.get('name', 'unnamed')} (run: false)")
            continue

        study_name = study.get('name', 'unnamed_study')
        logger.info(f"  Processing study panel: {study_name}")

        study_out = os.path.join(module_out, study_name)
        os.makedirs(study_out, exist_ok=True)

        # Initialize stats for this study
        study_stats = {
            'study_name': study_name,
            'total_snps': 0,
            'snps_present': 0,
            'snps_missing': 0,
            'fraction_present': 0.0
        }

        # Step 1: Load and validate SNP list
        snp_list_file = study.get('file')
        if not os.path.exists(snp_list_file):
            logger.error(f"    SNP list file not found: {snp_list_file}")
            continue

        id_type = study.get('id_type', 'position')
        has_header = study.get('has_header', False)

        if id_type == 'position':
            df_snps = _load_position_file(snp_list_file, has_header)
        elif id_type == 'rsid':
            if dbsnp_map is None:
                logger.error(f"    Cannot convert rsIDs for study {study_name} because no dbSNP mapping loaded.")
                continue
            df_snps = _convert_rsid_to_positions(snp_list_file, has_header, dbsnp_map, study_out)
        else:
            logger.error(f"    Unknown id_type '{id_type}' for study {study_name}. Skipping.")
            continue

        if df_snps is None or df_snps.empty:
            logger.error(f"    No valid SNPs loaded for study {study_name}. Skipping.")
            continue

        # Apply chromosome naming consistency (strip_chr if needed)
        if params.get('strip_chr', False):
            df_snps['chrom'] = df_snps['chrom'].str.replace('^chr', '', regex=True)

        study_stats['total_snps'] = len(df_snps)
        logger.info(f"    Loaded {len(df_snps)} SNPs for study {study_name}.")

        # Step 2: Create targets file for bcftools
        targets_file = os.path.join(study_out, "targets.txt")
        with open(targets_file, 'w') as f:
            for _, row in df_snps.iterrows():
                f.write(f"{row['chrom']}\t{row['pos']}\n")

        # Step 3: Intersect with filtered RNA VCF
        subset_vcf = os.path.join(study_out, "study_sites.vcf.gz")
        cmd = f"bcftools view -T {targets_file} {rna_filtered_vcf} -Oz -o {subset_vcf}"
        utils.run_cmd(cmd, f"Intersecting with study {study_name} SNPs")
        utils.run_cmd(f"bcftools index {subset_vcf}", f"Indexing subset VCF for {study_name}")

        # Step 4: Generate statistics
        stats_file = os.path.join(study_out, "mapping_stats.txt")
        present_count, missing_count = _generate_study_stats(df_snps, subset_vcf, stats_file, study_name)
        study_stats['snps_present'] = present_count
        study_stats['snps_missing'] = missing_count
        study_stats['fraction_present'] = present_count / study_stats['total_snps'] if study_stats['total_snps'] > 0 else 0.0

        # Step 5: Create detailed SNP table
        detail_file = os.path.join(study_out, "snp_details.tsv")
        _create_study_detail_table(df_snps, subset_vcf, rna_filtered_vcf, detail_file)

        # Store stats for this study
        stats[study_name] = study_stats
        logger.info(f"    Study {study_name} analysis completed. Results in {study_out}")

    # Update global statistics
    if stats:
        utils.update_stats(stats_path, 'study_snp_analysis', stats)

    logger.info(f"  Study SNP analysis module completed. Results in {module_out}")

def _load_position_file(file_path, has_header):
    """Load a position-based SNP file (chrom, pos, ref, alt)."""
    logger = utils.get_logger()
    try:
        if has_header:
            df = pd.read_csv(file_path, sep='\t')
            # Try to infer columns (common names: chrom, chr, pos, position, ref, alt)
            expected = ['chrom', 'pos', 'ref', 'alt']
            col_map = {}
            for col in expected:
                for variant in [col, col.upper(), col.capitalize()]:
                    if variant in df.columns:
                        col_map[col] = variant
                        break
                if col not in col_map:
                    raise ValueError(f"Required column '{col}' not found in file.")
            df = df[list(col_map.values())]
            df.columns = expected
        else:
            df = pd.read_csv(file_path, sep='\t', header=None,
                             names=['chrom', 'pos', 'ref', 'alt'])
        # Ensure pos is integer
        df['pos'] = df['pos'].astype(int)
        return df
    except Exception as e:
        logger.error(f"Error loading position file: {e}")
        return None

def _convert_rsid_to_positions(file_path, has_header, dbsnp_map, out_dir):
    """
    Convert rsIDs to chromosome positions using a pre-loaded dbSNP mapping DataFrame.
    The input file should contain rsIDs (one per line, optionally with a header).
    Returns a DataFrame with columns chrom, pos, ref, alt.
    """
    logger = utils.get_logger()
    try:
        if has_header:
            df_ids = pd.read_csv(file_path, sep='\t')
            # Assume first column contains rsIDs; allow column name 'rsid' or 'rs'
            if 'rsid' in df_ids.columns:
                rsids = df_ids['rsid'].astype(str)
            elif 'rs' in df_ids.columns:
                rsids = df_ids['rs'].astype(str)
            else:
                # take first column
                rsids = df_ids.iloc[:, 0].astype(str)
        else:
            rsids = pd.read_csv(file_path, sep='\t', header=None, names=['rsid'])['rsid'].astype(str)

        # Map to dbSNP
        mapped = dbsnp_map.reindex(rsids).dropna()
        if len(mapped) == 0:
            logger.error("    No rsIDs could be mapped to positions.")
            return None

        logger.info(f"    Mapped {len(mapped)} out of {len(rsids)} rsIDs to positions.")
        # Reset index to get rsid as column, then rename columns
        mapped = mapped.reset_index().rename(columns={'index': 'rsid'})
        # Keep only required columns
        result = mapped[['chrom', 'pos', 'ref', 'alt']].copy()
        return result
    except Exception as e:
        logger.error(f"Error converting rsIDs: {e}")
        return None

def _generate_study_stats(df_snps, subset_vcf, stats_file, study_name):
    """Generate presence/absence statistics."""
    logger = utils.get_logger()
    total_snps = len(df_snps)

    cmd = f"bcftools query -f '%CHROM\t%POS\n' {subset_vcf}"
    result = utils.run_cmd(cmd, f"Querying present SNPs for {study_name}")
    present_positions = set()
    for line in result.split('\n'):
        if line:
            chrom, pos = line.split()
            present_positions.add((chrom, int(pos)))

    present_count = len(present_positions)
    missing_count = total_snps - present_count

    with open(stats_file, 'w') as f:
        f.write(f"Total study SNPs\t{total_snps}\n")
        f.write(f"SNPs present in RNA VCF\t{present_count}\n")
        f.write(f"SNPs missing from RNA VCF\t{missing_count}\n")
        f.write(f"Fraction present\t{present_count/total_snps:.3f}\n")
    logger.info(f"    Statistics written to {stats_file}")
    return present_count, missing_count

def _create_study_detail_table(df_snps, subset_vcf, filtered_vcf, detail_file):
    """Create detailed per-SNP table (similar to variant_qc)."""
    logger = utils.get_logger()
    df_snps = df_snps.copy()
    df_snps['in_vcf'] = 'no'
    df_snps['vcf_filter'] = '.'
    df_snps['allele_match'] = 'NA'
    df_snps['vcf_ref'] = ''
    df_snps['vcf_alt'] = ''
    df_snps['f_missing'] = np.nan

    # Get info from subset VCF
    cmd = f"bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%FILTER\n' {subset_vcf}"
    result = utils.run_cmd(cmd, "Querying subset VCF")
    subset_data = {}
    for line in result.split('\n'):
        if line:
            parts = line.split('\t')
            if len(parts) >= 5:
                chrom, pos, ref, alt, filt = parts[:5]
                subset_data[(chrom, int(pos))] = (ref, alt, filt)

    # Get missing rates from filtered VCF (which may contain many samples)
    tmp_bed = os.path.join(os.path.dirname(detail_file), "tmp_targets.bed")
    with open(tmp_bed, 'w') as f:
        for _, row in df_snps.iterrows():
            f.write(f"{row['chrom']}\t{row['pos']-1}\t{row['pos']}\n")
    cmd = f"bcftools query -R {tmp_bed} -f '%CHROM\t%POS\t%F_MISSING\n' {filtered_vcf}"
    result = utils.run_cmd(cmd, "Querying missing rates")
    missing_data = {}
    for line in result.split('\n'):
        if line:
            chrom, pos, fmiss = line.split()
            missing_data[(chrom, int(pos))] = float(fmiss)
    if os.path.exists(tmp_bed):
        os.remove(tmp_bed)

    # Fill in details
    for idx, row in df_snps.iterrows():
        key = (row['chrom'], row['pos'])
        if key in subset_data:
            df_snps.at[idx, 'in_vcf'] = 'yes'
            vcf_ref, vcf_alt, vcf_filt = subset_data[key]
            df_snps.at[idx, 'vcf_ref'] = vcf_ref
            df_snps.at[idx, 'vcf_alt'] = vcf_alt
            df_snps.at[idx, 'vcf_filter'] = vcf_filt
            df_snps.at[idx, 'allele_match'] = 'yes' if (vcf_ref == row['ref'] and vcf_alt == row['alt']) else 'no'
        if key in missing_data:
            df_snps.at[idx, 'f_missing'] = missing_data[key]

    df_snps.to_csv(detail_file, sep='\t', index=False)
    logger.info(f"    Detailed SNP table saved to {detail_file}")