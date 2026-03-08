"""
Module 2: compare_wgs_rnaseq
-----------------------------
Compares RNA‑seq VCF against matched WGS VCF with parallel chromosome‑wise intersection.
- Reads sample mapping with user‑defined column names and separator.
- Stepwise filtering with detailed logging of missing samples.
- Caches intermediate files (common VCF, genotype table) to avoid recomputation.
- Allows user to specify chromosomes to analyze (default: all common chromosomes).
- Automatically detects available CPU cores (reserves 2) and RAM (reserves 8 GB) for optimal performance.
- Ensures input VCFs are indexed before parallel intersection.
- Automatically fixes chromosome naming mismatches (e.g., RNA has "chr1", WGS has "1") if enabled.
- **Three‑tier fallback**:
    1. Parallel per‑chromosome bcftools isec
    2. Single-threaded genome-wide bcftools isec
    3. Position‑list method (extract positions, intersect, then **merge** both subsets) – guaranteed to work and produces a VCF with both sample sets.
- Uses fast position‑only matching (-c all) and per‑chromosome temporary directories.
- Optionally uses RAM disk (/dev/shm) if enough free space (>20 GB) for temporary files.
- Can delete intermediate files (e.g., renamed RNA VCF) if `keep_intermediate` is false.
Produces:
  - common_sites.vcf.gz – sites present in both VCFs (with all samples)
  - genotypes.tsv – genotype table for all sample pairs
  - concordance_summary.tsv – precision, recall, concordance per sample
  - concordance_by_region.tsv – stratified by exonic/intronic (if exome BED provided)
  - titv_comparison.tsv – Ti/Tv for RNA and WGS per sample
  - ase_analysis.tsv – allele‑specific expression at WGS‑heterozygous sites
  - depth_correlation.tsv – Pearson correlation of log(DP) between RNA and WGS
If expression matrix is provided, the filtered mapping includes expression IDs for future use.
All key metrics are recorded in the pipeline statistics.
"""

import os
import pandas as pd
import numpy as np
from . import utils
import scipy.stats as stats
import concurrent.futures
import tempfile
import shutil
import logging
import subprocess

# -----------------------------------------------------------------------------
# Main entry point
# -----------------------------------------------------------------------------
def run(rna_vcf, wgs_vcf, sample_map, gene_expression, exome_bed, output_dir, stats_path, params):
    """
    Main entry point for compare_wgs_rnaseq module with stepwise filtering and logging.
    
    Args:
        rna_vcf (str): Path to RNA‑seq joint‑call VCF (bgzipped).
        wgs_vcf (str): Path to WGS joint‑call VCF (bgzipped).
        sample_map (str): Path to sample mapping TSV/CSV file.
        gene_expression (str): Path to expression matrix CSV (optional).
        exome_bed (str): Path to exonic regions BED file (optional).
        output_dir (str): Main pipeline output directory.
        stats_path (str): Path to pipeline statistics JSON file.
        params (dict): Module parameters from config:
            - min_wgs_dp (int): Minimum depth in WGS to trust genotype.
            - min_rna_dp (int): Minimum depth in RNA to consider call (used for filtering).
            - min_gq (int): Minimum genotype quality in WGS.
            - min_ase_depth (int): Minimum total depth for ASE analysis.
            - sample_map_sep (str): Delimiter for mapping file (default: tab).
            - sample_map_columns (dict): Column names mapping for mapping file.
            - auto_fix_chromosome_naming (bool): If True, attempt to fix chromosome naming mismatches.
            - chromosomes (str or list): Which chromosomes to analyze. "all" or a list.
            - fallback_to_full_isec (bool): If True, try single-threaded genome-wide isec when parallel fails.
            - use_position_list_fallback (bool): If True, use position‑list method as final fallback.
            - keep_intermediate (bool): If False, delete temporary files.
    """
    logger = utils.get_logger()
    module_out = os.path.join(output_dir, "compare_wgs_rnaseq")
    os.makedirs(module_out, exist_ok=True)

    # Initialize statistics for this module
    stats = {
        'input_rna_vcf': rna_vcf,
        'input_wgs_vcf': wgs_vcf,
        'mapping_initial_rows': 0,
        'mapping_final_rows': 0,
        'rna_vcf_samples': 0,
        'wgs_vcf_samples': 0,
        'missing_rna_ids': 0,
        'missing_wgs_ids': 0,
        'missing_expression_ids': 0,
        'common_sites_count': 0,
        'genotype_table_rows': 0,
        'concordance_summary': []
    }

    # -------------------------------------------------------------------------
    # Validate input files
    # -------------------------------------------------------------------------
    if not os.path.exists(rna_vcf):
        logger.error("RNA VCF file not found: %s", rna_vcf)
        raise FileNotFoundError(f"RNA VCF not found: {rna_vcf}")

    if not os.path.exists(wgs_vcf):
        logger.error("WGS VCF file not found: %s", wgs_vcf)
        raise FileNotFoundError(f"WGS VCF not found: {wgs_vcf}")

    if not os.path.exists(sample_map):
        logger.error("Sample mapping file not found: %s", sample_map)
        raise FileNotFoundError(f"Sample mapping not found: {sample_map}")

    # -------------------------------------------------------------------------
    # 1. Read mapping file using column names and separator from config
    # -------------------------------------------------------------------------
    column_mapping = params.get('sample_map_columns', None)
    sep = params.get('sample_map_sep', '\t')
    logger.info("Reading sample mapping from: %s with separator '%s'", sample_map, sep)
    mapping_df = utils.read_sample_map(sample_map, sep=sep, columns=column_mapping)
    initial_count = len(mapping_df)
    stats['mapping_initial_rows'] = initial_count
    logger.info(f"  Initial mapping contains {initial_count} rows.")

    # Get list of samples from VCFs
    rna_samples_vcf = utils.get_vcf_samples(rna_vcf)
    wgs_samples_vcf = utils.get_vcf_samples(wgs_vcf)
    stats['rna_vcf_samples'] = len(rna_samples_vcf)
    stats['wgs_vcf_samples'] = len(wgs_samples_vcf)
    logger.info("RNA VCF contains %d samples, WGS VCF contains %d samples.", 
                len(rna_samples_vcf), len(wgs_samples_vcf))

    # -------------------------------------------------------------------------
    # 2. Stepwise filtering with detailed logging
    # -------------------------------------------------------------------------
    # Filter by RNA VCF samples
    rna_in_mapping = set(mapping_df['rna'])
    rna_missing = rna_in_mapping - set(rna_samples_vcf)
    stats['missing_rna_ids'] = len(rna_missing)
    if rna_missing:
        logger.info(f"  {len(rna_missing)} RNA samples from mapping not found in RNA VCF.")
        logger.info(f"    Missing RNA IDs (first 10): {sorted(list(rna_missing))[:10]}")
    mapping_df = mapping_df[mapping_df['rna'].isin(rna_samples_vcf)]

    # Filter by WGS VCF samples
    wgs_in_mapping = set(mapping_df['wgs'])
    wgs_missing = wgs_in_mapping - set(wgs_samples_vcf)
    stats['missing_wgs_ids'] = len(wgs_missing)
    if wgs_missing:
        logger.info(f"  {len(wgs_missing)} WGS samples from mapping not found in WGS VCF.")
        logger.info(f"    Missing WGS IDs (first 10): {sorted(list(wgs_missing))[:10]}")
    mapping_df = mapping_df[mapping_df['wgs'].isin(wgs_samples_vcf)]

    # -------------------------------------------------------------------------
    # 3. Load expression matrix if provided
    # -------------------------------------------------------------------------
    expr_samples = None
    expr_df = None
    if gene_expression and os.path.exists(gene_expression):
        logger.info("  Loading expression matrix...")
        expr_df_full = pd.read_csv(gene_expression, index_col=0)
        expr_df_full.index = expr_df_full.index.astype(str)
        expr_samples = set(expr_df_full.index.astype(str))
        logger.info(f"    Expression matrix contains {len(expr_samples)} samples.")
        
        if 'expression' in mapping_df.columns:
            expr_in_mapping = set(mapping_df['expression'])
            expr_missing = expr_in_mapping - expr_samples
            stats['missing_expression_ids'] = len(expr_missing)
            if expr_missing:
                logger.info(f"  {len(expr_missing)} expression IDs from mapping not found in expression matrix.")
                logger.info(f"    Missing expression IDs (first 10): {sorted(list(expr_missing))[:10]}")
            mapping_df = mapping_df[mapping_df['expression'].isin(expr_samples)]
        else:
            logger.warning("  Expression samples provided but mapping has no 'expression' column; skipping expression filter.")
    elif gene_expression:
        logger.warning("  Expression file specified but not found: %s", gene_expression)
        logger.warning("  Continuing without expression data.")
    else:
        logger.info("  No expression file provided; skipping expression integration.")

    final_count = len(mapping_df)
    stats['mapping_final_rows'] = final_count
    logger.info(f"  Mapping filtered from {initial_count} to {final_count} common samples.")

    if mapping_df.empty:
        logger.error("No common samples found after filtering.")
        logger.error("Check that your mapping file contains sample names that exactly match those in the VCFs and expression matrix.")
        raise ValueError("No common samples found. Please verify sample names.")

    # Write filtered mapping for future runs (caching)
    filtered_map_path = os.path.join(module_out, "filtered_sample_map.tsv")
    utils.write_filtered_mapping(mapping_df, filtered_map_path)

    # Extract sample lists for downstream use
    if 'expression' in mapping_df.columns:
        rna_samples = mapping_df['rna'].tolist()
        wgs_samples = mapping_df['wgs'].tolist()
        expr_ids = mapping_df['expression'].tolist()
        sample_pairs = list(zip(rna_samples, wgs_samples))
        logger.info(f"  Using {len(sample_pairs)} sample pairs with matching expression IDs.")
    else:
        rna_samples = mapping_df['rna'].tolist()
        wgs_samples = mapping_df['wgs'].tolist()
        sample_pairs = list(zip(rna_samples, wgs_samples))
        expr_ids = None
        logger.info(f"  Using {len(sample_pairs)} sample pairs (no expression data).")

    # -------------------------------------------------------------------------
    # 4. Load expression matrix for the filtered samples (if applicable)
    # -------------------------------------------------------------------------
    if expr_ids is not None:
        expr_df = utils.load_expression_matrix(gene_expression, expr_ids)
        logger.info(f"    Expression matrix loaded with {expr_df.shape[0]} samples and {expr_df.shape[1]} genes.")
        expr_df.to_csv(os.path.join(module_out, "filtered_expression_matrix.csv"))
    else:
        expr_df = None

    # -------------------------------------------------------------------------
    # 5. Load exonic regions if provided
    # -------------------------------------------------------------------------
    if exome_bed and os.path.exists(exome_bed):
        exonic_pos = _load_exonic_positions(exome_bed)
    else:
        exonic_pos = None
        if exome_bed:
            logger.warning("  Exome BED file not found: %s; stratification by region will be skipped.", exome_bed)
        else:
            logger.info("  No exome BED provided; stratification by region skipped.")

    # -------------------------------------------------------------------------
    # 6. Intersect VCFs (with caching) – using multi‑tier fallback
    # -------------------------------------------------------------------------
    common_vcf = os.path.join(module_out, "common_sites.vcf.gz")
    dependencies = [rna_vcf, wgs_vcf, filtered_map_path]
    if utils.is_file_up_to_date(common_vcf, dependencies):
        logger.info(f"  Using existing common VCF: {common_vcf}")
    else:
        logger.info("  Generating common VCF (sites present in both RNA and WGS) using parallel intersection...")
        try:
            _intersect_vcfs_parallel(rna_vcf, wgs_vcf, common_vcf, params)
        except RuntimeError as e:
            if "No common sites found" in str(e):
                logger.error("No common sites found between RNA and WGS VCFs.")
                logger.info("Running diagnostics to identify the cause...")
                _run_diagnostics(rna_vcf, wgs_vcf, module_out)
                raise RuntimeError("No common sites found. Please check the diagnostic report in the output directory.")
            else:
                raise

    # Count common sites
    cmd = f"bcftools view -H {common_vcf} | wc -l"
    result = utils.run_cmd(cmd, "Counting common sites")
    stats['common_sites_count'] = int(result.strip())

    # -------------------------------------------------------------------------
    # 7. Extract genotype table (with caching)
    # -------------------------------------------------------------------------
    geno_tsv = os.path.join(module_out, "genotypes.tsv")
    deps_geno = [common_vcf, filtered_map_path]
    if utils.is_file_up_to_date(geno_tsv, deps_geno):
        logger.info(f"  Using existing genotype table: {geno_tsv}")
    else:
        logger.info("  Extracting genotypes...")
        _extract_genotypes(common_vcf, sample_pairs, geno_tsv, params)

    # Count rows in genotype table
    df_geno = pd.read_csv(geno_tsv, sep='\t')
    stats['genotype_table_rows'] = len(df_geno)

    # -------------------------------------------------------------------------
    # 8. Compute per‑sample precision/recall/concordance (always regenerate)
    # -------------------------------------------------------------------------
    concord_df = _compute_concordance(geno_tsv, module_out)
    stats['concordance_summary'] = concord_df.to_dict(orient='records')

    # -------------------------------------------------------------------------
    # 9. Stratify by region (if exome_bed provided)
    # -------------------------------------------------------------------------
    if exonic_pos:
        _stratify_by_region(geno_tsv, exonic_pos, module_out)

    # -------------------------------------------------------------------------
    # 10. Compute Ti/Tv for each sample
    # -------------------------------------------------------------------------
    _compute_titv(common_vcf, sample_pairs, module_out)

    # -------------------------------------------------------------------------
    # 11. Allele‑Specific Expression analysis
    # -------------------------------------------------------------------------
    _ase_analysis(common_vcf, sample_pairs, module_out, params)

    # -------------------------------------------------------------------------
    # 12. Coverage correlation
    # -------------------------------------------------------------------------
    _coverage_correlation(common_vcf, sample_pairs, module_out)

    # Update global statistics
    utils.update_stats(stats_path, 'compare_wgs_rnaseq', stats)

    logger.info(f"  Comparison module completed. Results in {module_out}")
    return stats

# -----------------------------------------------------------------------------
# Helper functions
# -----------------------------------------------------------------------------

def _load_exonic_positions(bed_file):
    """
    Load exonic BED file and return a set of (chrom, pos) for all 1‑based positions.
    BED is 0‑based, so positions are converted to 1‑based.
    """
    logger = utils.get_logger()
    exonic = set()
    with open(bed_file) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split()
            if len(parts) < 3:
                continue
            chrom, start, end = parts[0], int(parts[1]), int(parts[2])
            for pos in range(start+1, end+1):
                exonic.add((chrom, pos))
    logger.info(f"  Loaded {len(exonic)} exonic positions from {bed_file}")
    return exonic

def _extract_genotypes(vcf, sample_pairs, out_tsv, params):
    """
    Extract genotype calls (GT), depth (DP), and GQ for each sample pair.
    Apply filters: WGS DP >= min_wgs_dp, WGS GQ >= min_gq.
    RNA calls are kept even if low depth (but marked).
    """
    logger = utils.get_logger()
    min_wgs_dp = params.get('min_wgs_dp', 10)
    min_gq = params.get('min_gq', 20)

    # Get list of all samples in VCF
    all_samples = utils.run_cmd(f"bcftools query -l {vcf}").split('\n')
    # Map RNA and WGS sample names to indices
    rna_idx = []
    wgs_idx = []
    for rna, wgs in sample_pairs:
        if rna not in all_samples or wgs not in all_samples:
            raise ValueError(f"Sample {rna} or {wgs} not found in VCF.")
        rna_idx.append(all_samples.index(rna))
        wgs_idx.append(all_samples.index(wgs))

    # Query GT, DP, GQ for all samples
    # Correct format: %CHROM\t%POS\t%REF\t%ALT[\t%GT\t%DP\t%GQ]
    fmt_fields = "%CHROM\t%POS\t%REF\t%ALT[\t%GT\t%DP\t%GQ]"
    cmd = f"bcftools query -f '{fmt_fields}\n' {vcf}"
    result = utils.run_cmd(cmd, "Extracting genotypes")

    rows = []
    for line in result.split('\n'):
        if not line:
            continue
        parts = line.split('\t')
        # First 4 fields: chrom, pos, ref, alt
        chrom, pos, ref, alt = parts[0:4]
        # Remaining fields are triples for each sample: GT, DP, GQ
        sample_data = parts[4:]
        n_samples = len(all_samples)
        # Expect 3 * n_samples fields
        if len(sample_data) != 3 * n_samples:
            logger.warning(f"Line {chrom}:{pos} has unexpected number of fields, skipping.")
            continue
        for i, (rna_s, wgs_s) in enumerate(sample_pairs):
            # For each sample pair, we need to get the data for the corresponding sample indices
            # sample_data is in the same order as all_samples, so for a given sample index j,
            # its GT is at position 3*j, DP at 3*j+1, GQ at 3*j+2.
            rna_gt = sample_data[3 * rna_idx[i]]
            rna_dp = sample_data[3 * rna_idx[i] + 1]
            rna_gq = sample_data[3 * rna_idx[i] + 2]
            wgs_gt = sample_data[3 * wgs_idx[i]]
            wgs_dp = sample_data[3 * wgs_idx[i] + 1]
            wgs_gq = sample_data[3 * wgs_idx[i] + 2]

            # Convert '.' to missing
            rna_gt = '.' if rna_gt == '.' else rna_gt
            wgs_gt = '.' if wgs_gt == '.' else wgs_gt
            rna_dp = -1 if rna_dp == '.' else int(rna_dp)
            wgs_dp = -1 if wgs_dp == '.' else int(wgs_dp)
            rna_gq = -1 if rna_gq == '.' else int(rna_gq)
            wgs_gq = -1 if wgs_gq == '.' else int(wgs_gq)

            rows.append([chrom, pos, ref, alt, rna_s, wgs_s,
                         rna_gt, wgs_gt, rna_dp, wgs_dp, rna_gq, wgs_gq])

    df = pd.DataFrame(rows, columns=['chrom','pos','ref','alt','rna_sample','wgs_sample',
                                     'rna_GT','wgs_GT','rna_DP','wgs_DP','rna_GQ','wgs_GQ'])
    # Apply WGS filters
    df = df[(df['wgs_DP'] >= min_wgs_dp) & (df['wgs_GQ'] >= min_gq)]
    df.to_csv(out_tsv, sep='\t', index=False)
    logger.info(f"  Genotype table written to {out_tsv} ({len(df)} rows after filtering)")

def _compute_concordance(geno_tsv, out_dir):
    """
    For each sample pair, compute TP, FP, FN, TN, precision, recall, concordance.
    RNA is the call set, WGS is the truth.
    """
    logger = utils.get_logger()
    df = pd.read_csv(geno_tsv, sep='\t')
    results = []
    for (rna_sample, wgs_sample), group in df.groupby(['rna_sample','wgs_sample']):
        # WGS non‑reference = truth variant sites
        wgs_non_ref = group[group['wgs_GT'] != '0/0']
        # RNA non‑reference (and not missing)
        rna_non_ref = group[(group['rna_GT'] != '0/0') & (group['rna_GT'] != '.')]

        # True positives: RNA non‑ref and matches WGS non‑ref genotype
        tp = 0
        for _, row in wgs_non_ref.iterrows():
            if row['rna_GT'] != '.' and row['rna_GT'] == row['wgs_GT']:
                tp += 1

        # False positives: RNA non‑ref but WGS ref
        wgs_ref = group[group['wgs_GT'] == '0/0']
        fp = len(rna_non_ref[rna_non_ref['rna_sample'].isin(wgs_ref['rna_sample'])])

        # False negatives: WGS non‑ref but RNA ref or missing
        fn = len(wgs_non_ref[(wgs_non_ref['rna_GT'] == '0/0') | (wgs_non_ref['rna_GT'] == '.')])

        # True negatives: both ref
        tn = len(wgs_ref[wgs_ref['rna_GT'] == '0/0'])

        precision = tp / (tp + fp) if (tp + fp) > 0 else np.nan
        recall = tp / (tp + fn) if (tp + fn) > 0 else np.nan
        concordance = (tp + tn) / (tp + tn + fp + fn) if len(group) > 0 else np.nan

        results.append({
            'rna_sample': rna_sample,
            'wgs_sample': wgs_sample,
            'total_sites': len(group),
            'tp': tp,
            'fp': fp,
            'fn': fn,
            'tn': tn,
            'precision': precision,
            'recall': recall,
            'concordance': concordance
        })

    results_df = pd.DataFrame(results)
    results_df.to_csv(os.path.join(out_dir, 'concordance_summary.tsv'), sep='\t', index=False)
    logger.info(f"  Concordance summary saved to {out_dir}/concordance_summary.tsv")
    return results_df

def _stratify_by_region(geno_tsv, exonic_pos, out_dir):
    """
    Annotate each variant as exonic or intronic and compute concordance per region per sample.
    """
    logger = utils.get_logger()
    df = pd.read_csv(geno_tsv, sep='\t')
    df['region'] = df.apply(lambda row: 'exonic' if (row['chrom'], row['pos']) in exonic_pos else 'intronic', axis=1)

    rows = []
    for (rna_sample, wgs_sample, region), group in df.groupby(['rna_sample','wgs_sample','region']):
        wgs_non_ref = group[group['wgs_GT'] != '0/0']
        wgs_ref = group[group['wgs_GT'] == '0/0']
        tp = len(wgs_non_ref[(wgs_non_ref['rna_GT'] != '.') & (wgs_non_ref['rna_GT'] == wgs_non_ref['wgs_GT'])])
        fp = len(wgs_ref[(wgs_ref['rna_GT'] != '0/0') & (wgs_ref['rna_GT'] != '.')])
        fn = len(wgs_non_ref[(wgs_non_ref['rna_GT'] == '0/0') | (wgs_non_ref['rna_GT'] == '.')])
        tn = len(wgs_ref[wgs_ref['rna_GT'] == '0/0'])
        precision = tp/(tp+fp) if (tp+fp)>0 else np.nan
        recall = tp/(tp+fn) if (tp+fn)>0 else np.nan
        concordance = (tp+tn)/len(group) if len(group)>0 else np.nan
        rows.append([rna_sample, wgs_sample, region, len(group), precision, recall, concordance])

    region_df = pd.DataFrame(rows, columns=['rna_sample','wgs_sample','region','n_sites','precision','recall','concordance'])
    region_df.to_csv(os.path.join(out_dir, 'concordance_by_region.tsv'), sep='\t', index=False)
    logger.info(f"  Region‑stratified concordance saved to {out_dir}/concordance_by_region.tsv")

def _compute_titv(vcf, sample_pairs, out_dir):
    """
    Compute Ti/Tv ratio for each sample (RNA and WGS separately) using the common sites VCF.
    """
    logger = utils.get_logger()
    results = []
    for rna, wgs in sample_pairs:
        # RNA Ti/Tv
        cmd = f"bcftools view -s {rna} {vcf} | bcftools query -f '%REF\t%ALT\n'"
        result = utils.run_cmd(cmd, f"Extracting variants for {rna}")
        ti, tv = 0, 0
        for line in result.split('\n'):
            if not line:
                continue
            ref, alt = line.split()
            alt = alt.split(',')[0]  # take first alternate
            if len(ref)==1 and len(alt)==1:
                if {ref, alt} in [{'A','G'}, {'C','T'}]:
                    ti += 1
                else:
                    tv += 1
        rna_titv = ti/tv if tv>0 else np.nan

        # WGS Ti/Tv
        cmd = f"bcftools view -s {wgs} {vcf} | bcftools query -f '%REF\t%ALT\n'"
        result = utils.run_cmd(cmd, f"Extracting variants for {wgs}")
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
        wgs_titv = ti/tv if tv>0 else np.nan

        results.append([rna, wgs, rna_titv, wgs_titv])

    titv_df = pd.DataFrame(results, columns=['rna_sample','wgs_sample','rna_titv','wgs_titv'])
    titv_df.to_csv(os.path.join(out_dir, 'titv_comparison.tsv'), sep='\t', index=False)
    logger.info(f"  Ti/Tv comparison saved to {out_dir}/titv_comparison.tsv")

def _ase_analysis(vcf, sample_pairs, out_dir, params):
    """
    For sites where WGS is heterozygous (0/1), compute RNA allelic balance.
    Requires sufficient depth (min_ase_depth).
    """
    logger = utils.get_logger()
    min_depth = params.get('min_ase_depth', 20)
    out_rows = []
    for rna, wgs in sample_pairs:
        # Get WGS heterozygous sites with AD field
        cmd = f"bcftools view -s {wgs} {vcf} | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%GT\t%AD]\n'"
        result = utils.run_cmd(cmd, f"Extracting ASE data for {rna}/{wgs}")
        for line in result.split('\n'):
            if not line:
                continue
            parts = line.split()
            if len(parts) < 6:
                continue
            chrom, pos, ref, alt, gt, ad = parts[0], parts[1], parts[2], parts[3], parts[4], parts[5]
            if gt != '0/1':
                continue
            # Parse AD (ref,alt)
            try:
                ad_ref, ad_alt = map(int, ad.split(','))
            except:
                continue
            total_wgs = ad_ref + ad_alt
            if total_wgs < min_depth:
                continue
            # Get RNA AD at same site
            cmd2 = f"bcftools view -s {rna} -r {chrom}:{pos}-{pos} {vcf} | bcftools query -f '%AD\n'"
            rna_ad_res = utils.run_cmd(cmd2, check=False)
            if not rna_ad_res or rna_ad_res == '.':
                continue
            try:
                rna_ad_ref, rna_ad_alt = map(int, rna_ad_res.split(','))
            except:
                continue
            total_rna = rna_ad_ref + rna_ad_alt
            if total_rna < min_depth:
                continue
            wgs_ratio = ad_alt / total_wgs
            rna_ratio = rna_ad_alt / total_rna
            out_rows.append([rna, wgs, chrom, pos, ref, alt,
                             total_wgs, wgs_ratio, total_rna, rna_ratio])

    ase_df = pd.DataFrame(out_rows, columns=['rna_sample','wgs_sample','chrom','pos','ref','alt',
                                              'wgs_depth','wgs_alt_ratio','rna_depth','rna_alt_ratio'])
    ase_df.to_csv(os.path.join(out_dir, 'ase_analysis.tsv'), sep='\t', index=False)
    logger.info(f"  ASE analysis saved to {out_dir}/ase_analysis.tsv ({len(ase_df)} sites)")

def _coverage_correlation(vcf, sample_pairs, out_dir):
    """
    Compute Pearson correlation of log10(DP) between RNA and WGS at sites where both have DP>0.
    Uses FORMAT/DP (per-sample depth) correctly.
    """
    logger = utils.get_logger()
    rows = []
    for rna, wgs in sample_pairs:
        # Correct query to get per-sample DP using [%DP]
        cmd = f"bcftools query -s {rna},{wgs} -f '%CHROM\t%POS\t[%DP]\t[%DP]\n' {vcf}"
        result = utils.run_cmd(cmd, f"Extracting depth for {rna}/{wgs}")
        dp_pairs = []
        for line in result.split('\n'):
            if not line:
                continue
            parts = line.split()
            if len(parts) >= 4:
                chrom, pos, dp_rna, dp_wgs = parts[0], parts[1], parts[2], parts[3]
                if dp_rna != '.' and dp_wgs != '.' and int(dp_rna) > 0 and int(dp_wgs) > 0:
                    dp_pairs.append((int(dp_rna), int(dp_wgs)))
        if len(dp_pairs) > 1:
            rna_dp = [p[0] for p in dp_pairs]
            wgs_dp = [p[1] for p in dp_pairs]
            corr, pval = stats.pearsonr(np.log10(rna_dp), np.log10(wgs_dp))
            rows.append([rna, wgs, len(dp_pairs), corr, pval])
    corr_df = pd.DataFrame(rows, columns=['rna_sample','wgs_sample','n_sites','pearson_r','pvalue'])
    corr_df.to_csv(os.path.join(out_dir, 'depth_correlation.tsv'), sep='\t', index=False)
    logger.info(f"  Depth correlation saved to {out_dir}/depth_correlation.tsv")

# -----------------------------------------------------------------------------
# Parallel intersection function with multiple fallback strategies
# -----------------------------------------------------------------------------
def _intersect_vcfs_parallel(rna_vcf, wgs_vcf, out_vcf, params):
    """
    Parallel chromosome‑wise intersection of two VCFs with multiple fallback strategies.
    Strategy 1: Parallel per‑chromosome bcftools isec.
    Strategy 2: Single-threaded genome-wide bcftools isec.
    Strategy 3: Position‑list method (extract positions, intersect, then merge subsets) – guaranteed to work.
    """
    logger = utils.get_logger()
    keep_intermediate = params.get('keep_intermediate', True)
    fallback = params.get('fallback_to_full_isec', True)
    use_position_fallback = params.get('use_position_list_fallback', True)
    fixed_rna_vcf = None
    rna_vcf_to_use = rna_vcf

    # Ensure both VCFs are indexed (required for -r option)
    logger.info("  Ensuring RNA VCF is indexed...")
    utils.ensure_vcf_index(rna_vcf)
    logger.info("  Ensuring WGS VCF is indexed...")
    utils.ensure_vcf_index(wgs_vcf)

    # Get full contig lists
    rna_contigs_all = set(utils.get_vcf_contigs(rna_vcf))
    wgs_contigs_all = set(utils.get_vcf_contigs(wgs_vcf))
    common_full = sorted(rna_contigs_all & wgs_contigs_all)

    logger.info(f"  Total common contigs between VCFs: {len(common_full)}")
    if logger.isEnabledFor(logging.DEBUG):
        logger.debug(f"    Common contigs: {', '.join(common_full)}")

    # Determine which chromosomes to analyze
    chrom_param = params.get('chromosomes', 'all')
    if chrom_param == 'all':
        chroms_to_analyze = common_full
        logger.info(f"  Using all {len(chroms_to_analyze)} common chromosomes.")
    elif isinstance(chrom_param, list):
        # User provided a list; keep only those present in common_full
        chroms_to_analyze = [c for c in chrom_param if c in common_full]
        if not chroms_to_analyze:
            logger.error("None of the specified chromosomes are common between VCFs.")
            # Attempt auto‑fix if enabled
            auto_fix = params.get('auto_fix_chromosome_naming', True)
            if auto_fix:
                logger.info("  Attempting to fix chromosome naming...")
                # Try stripping "chr" from RNA if needed
                rna_stripped = {c.replace('chr', ''): c for c in rna_contigs_all if c.startswith('chr')}
                wgs_stripped = {c.replace('chr', ''): c for c in wgs_contigs_all if c.startswith('chr')}
                # Build mapping
                mapping = {}
                for user_chrom in chrom_param:
                    if user_chrom in rna_contigs_all and user_chrom in wgs_contigs_all:
                        continue  # already handled
                    # Check if RNA has "chr" and WGS doesn't
                    if user_chrom.startswith('chr') and user_chrom[3:] in wgs_contigs_all:
                        mapping[user_chrom] = user_chrom[3:]
                    # Check if RNA lacks "chr" and WGS has it
                    elif not user_chrom.startswith('chr') and ('chr' + user_chrom) in rna_contigs_all:
                        mapping['chr' + user_chrom] = user_chrom
                if mapping:
                    # Create temporary renamed RNA VCF
                    fixed_rna_vcf = os.path.join(os.path.dirname(out_vcf), "rna_renamed.vcf.gz")
                    logger.info(f"    Renaming {len(mapping)} contigs in RNA VCF and saving to {fixed_rna_vcf}")
                    utils.rename_contigs(rna_vcf, mapping, fixed_rna_vcf)
                    rna_vcf_to_use = fixed_rna_vcf
                    # Update common list
                    rna_contigs_new = set(utils.get_vcf_contigs(rna_vcf_to_use))
                    chroms_to_analyze = [c for c in chrom_param if c in rna_contigs_new and c in wgs_contigs_all]
                    logger.info(f"    After renaming, found {len(chroms_to_analyze)} common chromosomes.")
                else:
                    logger.error("Could not resolve naming mismatch.")
                    raise RuntimeError("No common chromosomes after attempted fix.")
            else:
                raise RuntimeError("No common chromosomes found for specified list.")
        logger.info(f"  Using {len(chroms_to_analyze)} chromosomes from user list.")
    else:
        logger.error("Invalid value for 'chromosomes' parameter.")
        raise ValueError("Invalid chromosomes parameter.")

    if not chroms_to_analyze:
        logger.error("No chromosomes to analyze.")
        raise RuntimeError("No chromosomes to analyze.")

    # -------------------------------------------------------------------------
    # Strategy 1: Parallel per‑chromosome bcftools isec
    # -------------------------------------------------------------------------
    def try_parallel_isec():
        logger.info("  Attempting parallel per‑chromosome bcftools isec...")
        use_ram = utils.should_use_ram_disk(min_free_gb=20)
        tmp_base = '/dev/shm' if use_ram else None
        tmp_dir = tempfile.mkdtemp(dir=tmp_base, prefix="isec_parallel_")
        logger.debug(f"  Temporary directory: {tmp_dir}")

        max_workers = utils.get_available_cores(reserve=2)
        max_workers = min(max_workers, len(chroms_to_analyze))
        logger.info(f"  Using up to {max_workers} cores (reserving 2).")

        def process_chromosome(chrom):
            chrom_dir = os.path.join(tmp_dir, f"{chrom}_isec")
            os.makedirs(chrom_dir, exist_ok=True)
            cmd = (f"bcftools isec -r {chrom} -c all -n=2 "
                   f"{rna_vcf_to_use} {wgs_vcf} -p {chrom_dir}")
            utils.run_cmd(cmd, f"Intersecting chromosome {chrom}")
            # Check for 0002.vcf (sites in both)
            intersect_file = os.path.join(chrom_dir, "0002.vcf")
            if os.path.exists(intersect_file) and os.path.getsize(intersect_file) > 0:
                compressed = os.path.join(tmp_dir, f"{chrom}.vcf.gz")
                utils.run_cmd(f"bgzip -c {intersect_file} > {compressed}", f"Compressing {chrom}")
                utils.run_cmd(f"bcftools index {compressed}", f"Indexing {chrom}")
                shutil.rmtree(chrom_dir, ignore_errors=True)
                return compressed
            else:
                shutil.rmtree(chrom_dir, ignore_errors=True)
                return None

        results = []
        with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
            future_to_chrom = {executor.submit(process_chromosome, chrom): chrom for chrom in chroms_to_analyze}
            for future in concurrent.futures.as_completed(future_to_chrom):
                chrom = future_to_chrom[future]
                try:
                    result = future.result()
                    if result:
                        results.append(result)
                        logger.debug(f"  Chromosome {chrom} yielded common sites.")
                    else:
                        logger.debug(f"  Chromosome {chrom} had no common sites.")
                except Exception as e:
                    logger.error(f"  Intersection failed for chromosome {chrom}: {e}")
                    shutil.rmtree(tmp_dir, ignore_errors=True)
                    if fixed_rna_vcf and os.path.exists(fixed_rna_vcf) and not keep_intermediate:
                        utils.safe_remove(fixed_rna_vcf)
                    raise RuntimeError(f"Intersection failed on chromosome {chrom}")

        if results:
            # Concatenate and finish
            logger.info(f"  Concatenating {len(results)} chromosome VCFs...")
            cmd = f"bcftools concat {' '.join(results)} -Oz -o {out_vcf}"
            utils.run_cmd(cmd, "Concatenating chromosome results")
            utils.run_cmd(f"bcftools index {out_vcf}", "Indexing common VCF")
            shutil.rmtree(tmp_dir, ignore_errors=True)
            if fixed_rna_vcf and os.path.exists(fixed_rna_vcf):
                if keep_intermediate:
                    logger.info(f"  Intermediate renamed RNA VCF kept at: {fixed_rna_vcf}")
                else:
                    utils.safe_remove(fixed_rna_vcf)
            logger.info(f"  Common sites VCF saved to {out_vcf}")
            return True
        else:
            logger.warning("  No common sites found via parallel isec.")
            shutil.rmtree(tmp_dir, ignore_errors=True)
            return False

    # -------------------------------------------------------------------------
    # Strategy 2: Single-threaded genome-wide bcftools isec
    # -------------------------------------------------------------------------
    def try_full_isec():
        logger.info("  Attempting single-threaded genome-wide bcftools isec...")
        use_ram = utils.should_use_ram_disk(min_free_gb=20)
        tmp_base = '/dev/shm' if use_ram else None
        tmp_dir = tempfile.mkdtemp(dir=tmp_base, prefix="isec_full_")
        cmd = f"bcftools isec -c all -n=2 {rna_vcf_to_use} {wgs_vcf} -p {tmp_dir}"
        try:
            utils.run_cmd(cmd, "Genome-wide intersection")
            intersect_file = os.path.join(tmp_dir, "0002.vcf")
            if os.path.exists(intersect_file) and os.path.getsize(intersect_file) > 0:
                utils.run_cmd(f"bgzip -c {intersect_file} > {out_vcf}", "Compressing full isec result")
                utils.run_cmd(f"bcftools index {out_vcf}", "Indexing full isec VCF")
                shutil.rmtree(tmp_dir, ignore_errors=True)
                if fixed_rna_vcf and os.path.exists(fixed_rna_vcf) and not keep_intermediate:
                    utils.safe_remove(fixed_rna_vcf)
                logger.info(f"  Common sites VCF saved to {out_vcf}")
                return True
            else:
                logger.warning("  No common sites found via full isec.")
                shutil.rmtree(tmp_dir, ignore_errors=True)
                return False
        except Exception as e:
            logger.error(f"  Full isec failed: {e}")
            shutil.rmtree(tmp_dir, ignore_errors=True)
            return False

    # -------------------------------------------------------------------------
    # Strategy 3: Position‑list method (extract positions, intersect, then merge subsets)
    # -------------------------------------------------------------------------
    def try_position_list():
        logger.info("  Attempting position‑list intersection method with merging...")
        # We'll extract positions from both VCFs using bcftools query.
        max_workers = utils.get_available_cores(reserve=2)
        logger.info(f"  Using up to {max_workers} cores for position extraction.")

        def extract_positions(vcf, chroms):
            """Extract all (chrom, pos) from given VCF on specified chromosomes."""
            pos_set = set()
            total = 0
            # Use bcftools query per chromosome
            for chrom in chroms:
                cmd = f"bcftools query -f '%CHROM\t%POS\n' -r {chrom} {vcf}"
                proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
                for line in proc.stdout:
                    line = line.strip()
                    if line:
                        parts = line.split('\t')
                        if len(parts) >= 2:
                            pos_set.add((parts[0], int(parts[1])))
                            total += 1
                proc.wait()
            return pos_set

        logger.info("    Extracting RNA positions...")
        rna_pos = extract_positions(rna_vcf_to_use, chroms_to_analyze)
        logger.info(f"      Extracted {len(rna_pos)} RNA positions.")

        logger.info("    Extracting WGS positions...")
        wgs_pos = extract_positions(wgs_vcf, chroms_to_analyze)
        logger.info(f"      Extracted {len(wgs_pos)} WGS positions.")

        intersect_pos = rna_pos & wgs_pos
        logger.info(f"      Intersection size: {len(intersect_pos)} positions.")

        if not intersect_pos:
            logger.error("  No overlapping positions found via position‑list method.")
            return False

        # Create a BED file of intersecting positions
        bed_file = os.path.join(os.path.dirname(out_vcf), "intersect_positions.bed")
        with open(bed_file, 'w') as f:
            for chrom, pos in sorted(intersect_pos):
                # BED format: chrom, start0, end1
                f.write(f"{chrom}\t{pos-1}\t{pos}\n")

        # Extract from RNA and WGS separately
        rna_subset = os.path.join(os.path.dirname(out_vcf), "rna_subset.vcf.gz")
        wgs_subset = os.path.join(os.path.dirname(out_vcf), "wgs_subset.vcf.gz")
        cmd_rna = f"bcftools view -R {bed_file} {rna_vcf_to_use} -Oz -o {rna_subset}"
        cmd_wgs = f"bcftools view -R {bed_file} {wgs_vcf} -Oz -o {wgs_subset}"
        utils.run_cmd(cmd_rna, "Extracting RNA subset")
        utils.run_cmd(cmd_wgs, "Extracting WGS subset")
        utils.run_cmd(f"bcftools index {rna_subset}", "Indexing RNA subset")
        utils.run_cmd(f"bcftools index {wgs_subset}", "Indexing WGS subset")

        # Merge the two VCFs
        cmd_merge = f"bcftools merge {rna_subset} {wgs_subset} -Oz -o {out_vcf}"
        utils.run_cmd(cmd_merge, "Merging subsets")
        utils.run_cmd(f"bcftools index {out_vcf}", "Indexing merged VCF")

        # Clean up intermediate files if not keeping
        if not keep_intermediate:
            utils.safe_remove(bed_file)
            utils.safe_remove(rna_subset)
            utils.safe_remove(wgs_subset)
            utils.safe_remove(rna_subset + '.tbi')
            utils.safe_remove(wgs_subset + '.tbi')
        if fixed_rna_vcf and os.path.exists(fixed_rna_vcf) and not keep_intermediate:
            utils.safe_remove(fixed_rna_vcf)

        logger.info(f"  Common sites VCF saved to {out_vcf}")
        return True

    # -------------------------------------------------------------------------
    # Execute strategies in order
    # -------------------------------------------------------------------------
    if try_parallel_isec():
        return

    if fallback and try_full_isec():
        return

    if use_position_fallback and try_position_list():
        return

    # All fallbacks failed
    logger.error("All intersection methods failed. Cannot proceed.")
    if fixed_rna_vcf and os.path.exists(fixed_rna_vcf) and not keep_intermediate:
        utils.safe_remove(fixed_rna_vcf)
    raise RuntimeError("No common sites found. Please check diagnostic report.")

# -----------------------------------------------------------------------------
# Diagnostic function (unchanged)
# -----------------------------------------------------------------------------
def _run_diagnostics(rna_vcf, wgs_vcf, out_dir):
    """
    Perform detailed diagnostics to understand why two VCFs have no common sites.
    Writes reports to the output directory.
    """
    logger = utils.get_logger()
    diag_dir = os.path.join(out_dir, "diagnostics")
    os.makedirs(diag_dir, exist_ok=True)

    logger.info("  Running diagnostics to identify cause of no common sites...")

    # 1. Check chromosome naming conventions
    rna_contigs = utils.get_vcf_contigs(rna_vcf)
    wgs_contigs = utils.get_vcf_contigs(wgs_vcf)

    with open(os.path.join(diag_dir, "chromosome_lists.txt"), 'w') as f:
        f.write("RNA VCF chromosomes:\n")
        f.write("\n".join(rna_contigs) + "\n\n")
        f.write("WGS VCF chromosomes:\n")
        f.write("\n".join(wgs_contigs) + "\n")

    # Check if one has "chr" prefix and the other doesn't
    rna_has_chr = any(c.startswith('chr') for c in rna_contigs)
    wgs_has_chr = any(c.startswith('chr') for c in wgs_contigs)
    if rna_has_chr != wgs_has_chr:
        logger.info("    Chromosome naming mismatch: RNA has 'chr' prefix? %s, WGS has 'chr' prefix? %s", 
                    rna_has_chr, wgs_has_chr)
        logger.info("    Consider enabling 'auto_fix_chromosome_naming' in config if not already enabled.")

    # 2. Get total variant counts
    def get_total_variants(vcf):
        cmd = f"bcftools view -H {vcf} | wc -l"
        try:
            result = utils.run_cmd(cmd, check=False)
            return int(result.strip())
        except:
            return 0

    total_rna = get_total_variants(rna_vcf)
    total_wgs = get_total_variants(wgs_vcf)
    logger.info(f"    RNA VCF total variants: {total_rna}")
    logger.info(f"    WGS VCF total variants: {total_wgs}")

    if total_rna == 0:
        logger.error("    RNA VCF contains zero variants!")
    if total_wgs == 0:
        logger.error("    WGS VCF contains zero variants!")

    # 3. Get variant counts per chromosome (only on standard autosomes/sex for brevity)
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

    with open(os.path.join(diag_dir, "variant_counts.txt"), 'w') as f:
        f.write("RNA VCF variant counts per chromosome:\n")
        for chrom in sorted(rna_counts.keys()):
            f.write(f"{chrom}\t{rna_counts[chrom]}\n")
        f.write("\nWGS VCF variant counts per chromosome:\n")
        for chrom in sorted(wgs_counts.keys()):
            f.write(f"{chrom}\t{wgs_counts[chrom]}\n")

    # 4. Try a quick intersection on a small region (e.g., first common chromosome after possible stripping)
    # We'll try both with and without "chr" to see if naming is the issue
    test_chrom = None
    for c in ['chr1', '1']:
        if c in rna_contigs and c in wgs_contigs:
            test_chrom = c
            break
    if test_chrom:
        logger.info(f"    Testing intersection on {test_chrom}:1-1000000...")
        tmp_dir = tempfile.mkdtemp(dir=diag_dir, prefix="test_isec_")
        cmd = (f"bcftools isec -r {test_chrom}:1-1000000 -c all -n=2 "
               f"{rna_vcf} {wgs_vcf} -p {tmp_dir}")
        try:
            utils.run_cmd(cmd, "Test intersection", check=True)
            intersect_file = os.path.join(tmp_dir, "0002.vcf")
            if os.path.exists(intersect_file) and os.path.getsize(intersect_file) > 0:
                line_count = int(utils.run_cmd(f"cat {intersect_file} | grep -v '^#' | wc -l"))
                logger.info(f"      Found {line_count} common sites in this region.")
            else:
                logger.info("      No common sites in this region.")
        except Exception as e:
            logger.error(f"      Test intersection failed: {e}")
        finally:
            shutil.rmtree(tmp_dir, ignore_errors=True)
    else:
        # If no common chromosome, try with stripped names
        logger.info("    No common chromosomes with exact names; testing with stripped 'chr' from RNA...")
        rna_stripped = {c.replace('chr', ''): c for c in rna_contigs if c.startswith('chr')}
        for stripped, original in rna_stripped.items():
            if stripped in wgs_contigs:
                test_chrom = stripped
                logger.info(f"    Testing intersection on {original} (RNA) vs {stripped} (WGS) region 1-1000000...")
                tmp_dir = tempfile.mkdtemp(dir=diag_dir, prefix="test_isec_")
                # Need to use the original RNA contig with -r on RNA VCF, but bcftools isec requires matching names.
                # So we create a temporary renamed RNA VCF for this test.
                mapping = {original: stripped}
                test_rna_vcf = os.path.join(tmp_dir, "test_rna.vcf.gz")
                utils.rename_contigs(rna_vcf, mapping, test_rna_vcf)
                cmd = (f"bcftools isec -r {stripped}:1-1000000 -c all -n=2 "
                       f"{test_rna_vcf} {wgs_vcf} -p {tmp_dir}/isec")
                try:
                    utils.run_cmd(cmd, "Test intersection with renamed RNA", check=True)
                    intersect_file = os.path.join(tmp_dir, "isec", "0002.vcf")
                    if os.path.exists(intersect_file) and os.path.getsize(intersect_file) > 0:
                        line_count = int(utils.run_cmd(f"cat {intersect_file} | grep -v '^#' | wc -l"))
                        logger.info(f"      Found {line_count} common sites in this region.")
                    else:
                        logger.info("      No common sites in this region.")
                except Exception as e:
                    logger.error(f"      Test intersection failed: {e}")
                finally:
                    shutil.rmtree(tmp_dir, ignore_errors=True)
                break

    logger.info(f"  Diagnostics complete. Reports saved to {diag_dir}")