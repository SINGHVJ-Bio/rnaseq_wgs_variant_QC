"""
Shared utility functions for the pipeline.
Includes logging setup, sample filtering, file caching, publication styling,
expression matrix handling, dbSNP mapping, VCF header parsing, contig retrieval,
system resource detection, and VCF contig renaming.
"""

import subprocess
import pandas as pd
import numpy as np
import os
import time
import logging
import json
import shutil
import tempfile
import psutil  # for memory and CPU detection

# -----------------------------------------------------------------------------
# Logging setup
# -----------------------------------------------------------------------------
def setup_logging(log_path, level=logging.INFO):
    """
    Configure logging to write to both file and console.
    Returns a logger object.
    """
    logger = logging.getLogger('rnaseq_wgs_pipeline')
    logger.setLevel(level)

    # Remove any existing handlers to avoid duplication
    if logger.hasHandlers():
        logger.handlers.clear()

    # File handler
    fh = logging.FileHandler(log_path)
    fh.setLevel(level)
    # Console handler
    ch = logging.StreamHandler()
    ch.setLevel(level)

    # Formatter
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    fh.setFormatter(formatter)
    ch.setFormatter(formatter)

    logger.addHandler(fh)
    logger.addHandler(ch)

    return logger

def get_logger():
    """Return the existing logger, or create a dummy one if not configured."""
    logger = logging.getLogger('rnaseq_wgs_pipeline')
    if not logger.hasHandlers():
        # Add a null handler to avoid "No handlers found" warning
        logger.addHandler(logging.NullHandler())
    return logger

# -----------------------------------------------------------------------------
# System resource detection
# -----------------------------------------------------------------------------
def get_available_cores(reserve=2):
    """
    Return the number of CPU cores to use, reserving `reserve` cores for system.
    Ensures at least 1 core is returned.
    """
    total = os.cpu_count()
    if total is None:
        total = 4  # conservative fallback
    available = max(1, total - reserve)
    logger = get_logger()
    logger.debug(f"Total cores: {total}, using {available} (reserving {reserve})")
    return available

def get_available_memory(reserve_gb=8):
    """
    Return available memory in bytes, considering a reserve.
    Uses psutil to get available RAM.
    Returns None if psutil not available or detection fails.
    """
    try:
        mem = psutil.virtual_memory()
        available_bytes = mem.available
        reserve_bytes = reserve_gb * 1024**3
        usable = max(0, available_bytes - reserve_bytes)
        logger = get_logger()
        logger.debug(f"Total available RAM: {available_bytes/(1024**3):.1f} GB, using up to {usable/(1024**3):.1f} GB (reserving {reserve_gb} GB)")
        return usable
    except ImportError:
        logger = get_logger()
        logger.warning("psutil not installed; cannot detect available memory. Install with: pip install psutil")
        return None
    except Exception as e:
        logger = get_logger()
        logger.warning(f"Could not detect available memory: {e}")
        return None

def should_use_ram_disk(min_free_gb=20):
    """
    Check if /dev/shm exists and has at least `min_free_gb` GB free.
    Returns True if suitable, False otherwise.
    """
    if not os.path.exists('/dev/shm'):
        return False
    try:
        usage = shutil.disk_usage('/dev/shm')
        free_gb = usage.free / (1024**3)
        if free_gb >= min_free_gb:
            logger = get_logger()
            logger.debug(f"/dev/shm has {free_gb:.1f} GB free; using as RAM disk.")
            return True
        else:
            return False
    except Exception:
        return False

# -----------------------------------------------------------------------------
# Shell command execution
# -----------------------------------------------------------------------------
def run_cmd(cmd, description="", check=True):
    """
    Run a shell command and return stdout.
    If check=True, raise exception on non-zero exit.
    """
    logger = get_logger()
    logger.debug(f"Running command: {cmd}")
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if check and result.returncode != 0:
        logger.error(f"{description} failed:\n{result.stderr}")
        raise RuntimeError(f"{description} failed:\n{result.stderr}")
    return result.stdout.strip()

# -----------------------------------------------------------------------------
# Statistics tracking
# -----------------------------------------------------------------------------
def init_stats(output_dir):
    """
    Initialize the pipeline statistics file.
    Returns the path to the stats file.
    """
    stats_path = os.path.join(output_dir, "pipeline_stats.json")
    if not os.path.exists(stats_path):
        with open(stats_path, 'w') as f:
            json.dump({}, f)
    return stats_path

def update_stats(stats_path, module, data):
    """
    Update the statistics file with data from a module.
    """
    logger = get_logger()
    try:
        with open(stats_path, 'r') as f:
            stats = json.load(f)
    except (FileNotFoundError, json.JSONDecodeError):
        stats = {}
    
    if module not in stats:
        stats[module] = {}
    stats[module].update(data)
    
    with open(stats_path, 'w') as f:
        json.dump(stats, f, indent=2)
    logger.debug(f"Updated stats for module {module}")

def write_stats_summary(stats_path, output_dir):
    """
    Generate a human-readable TSV summary of all statistics.
    """
    logger = get_logger()
    try:
        with open(stats_path, 'r') as f:
            stats = json.load(f)
    except (FileNotFoundError, json.JSONDecodeError):
        logger.warning("No statistics file found.")
        return
    
    summary_file = os.path.join(output_dir, "pipeline_stats_summary.tsv")
    with open(summary_file, 'w') as f:
        f.write("Module\tMetric\tValue\n")
        for module, metrics in stats.items():
            for key, value in metrics.items():
                # Convert dict/list to string for display
                if isinstance(value, (dict, list)):
                    value = json.dumps(value)
                f.write(f"{module}\t{key}\t{value}\n")
    logger.info(f"Statistics summary written to {summary_file}")

# -----------------------------------------------------------------------------
# VCF header parsing and contig handling
# -----------------------------------------------------------------------------
def get_vcf_info_fields(vcf):
    """
    Return a set of INFO field names present in the VCF header.
    """
    cmd = f"bcftools view -h {vcf} | grep '^##INFO=<ID=' | sed 's/^##INFO=<ID=\\([^,]*\\),.*$/\\1/'"
    result = run_cmd(cmd, "Extracting INFO fields", check=False)
    fields = set()
    for line in result.split('\n'):
        if line.strip():
            fields.add(line.strip())
    return fields

def get_vcf_samples(vcf):
    """Return list of sample names from a VCF."""
    return run_cmd(f"bcftools query -l {vcf}").split('\n')

def get_vcf_contigs(vcf):
    """
    Return a sorted list of chromosome names (contigs) present in the VCF.
    Uses bcftools index to list contigs (requires .csi or .tbi).
    """
    logger = get_logger()
    # Try to get contigs from index
    cmd = f"bcftools index -s {vcf} | cut -f1 | sort -V"
    result = run_cmd(cmd, "Extracting contigs", check=False)
    if result:
        contigs = result.split('\n')
        logger.debug(f"Found {len(contigs)} contigs in {vcf}")
        return contigs
    else:
        # Fallback: extract from header
        cmd = f"bcftools view -h {vcf} | grep '^##contig=<ID=' | sed 's/^##contig=<ID=\\([^,]*\\),.*$/\\1/'"
        result = run_cmd(cmd, "Extracting contigs from header", check=True)
        contigs = result.split('\n')
        logger.debug(f"Found {len(contigs)} contigs in header of {vcf}")
        return contigs

def rename_contigs(vcf, mapping, out_vcf):
    """
    Rename contigs in a VCF using a mapping dictionary.
    
    Args:
        vcf (str): Input VCF file (bgzipped).
        mapping (dict): Dictionary mapping old contig names to new names.
        out_vcf (str): Output VCF file (bgzipped).
    """
    logger = get_logger()
    # Create a temporary file with the mapping
    with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as map_file:
        map_path = map_file.name
        for old, new in mapping.items():
            map_file.write(f"{old} {new}\n")
    try:
        cmd = f"bcftools annotate --rename-chrs {map_path} {vcf} -Oz -o {out_vcf}"
        run_cmd(cmd, f"Renaming contigs in {vcf}")
        run_cmd(f"bcftools index {out_vcf}", f"Indexing renamed VCF {out_vcf}")
    finally:
        os.unlink(map_path)
    logger.info(f"Renamed contigs in {vcf} -> {out_vcf} using mapping of {len(mapping)} entries")

def safe_remove(path):
    """
    Safely remove a file or directory. Logs a warning if removal fails.
    """
    logger = get_logger()
    if path is None:
        return
    try:
        if os.path.isfile(path):
            os.unlink(path)
            logger.debug(f"Removed file: {path}")
        elif os.path.isdir(path):
            shutil.rmtree(path)
            logger.debug(f"Removed directory: {path}")
    except Exception as e:
        logger.warning(f"Could not remove {path}: {e}")
        
# -----------------------------------------------------------------------------
# Sample and VCF handling (existing functions)
# -----------------------------------------------------------------------------
def read_sample_map(path, sep=None, columns=None):
    """
    Read sample mapping TSV/CSV.
    
    Args:
        path (str): Path to mapping file.
        sep (str, optional): Delimiter. If None, attempts to detect (default: tab).
        columns (dict, optional): Dictionary with keys 'wgs', 'rna', and optionally 'expression'
                                   mapping to column names in the file. If not provided, the
                                   function falls back to positional interpretation (first 2 or 3 columns).
    
    Returns:
        DataFrame with columns 'wgs', 'rna' and optionally 'expression'.
    """
    logger = get_logger()
    
    # Use provided separator or default to tab
    if sep is None:
        sep = '\t'
    
    # Read with header (we'll detect later)
    df = pd.read_csv(path, sep=sep)
    
    # If columns dict is provided, use it to map
    if columns is not None:
        # Check required columns
        required = ['wgs', 'rna']
        missing = [col for col in required if col not in columns]
        if missing:
            raise ValueError(f"Missing required column keys in sample_map_columns: {missing}")
        
        # Map column names
        col_map = {}
        for key, col_name in columns.items():
            if col_name not in df.columns:
                raise ValueError(f"Column '{col_name}' specified for '{key}' not found in mapping file.")
            col_map[key] = col_name
        
        # Keep only mapped columns
        keep_cols = [col_map['wgs'], col_map['rna']]
        if 'expression' in col_map:
            keep_cols.append(col_map['expression'])
        df = df[keep_cols]
        
        # Rename to standard names
        rename_dict = {col_map['wgs']: 'wgs', col_map['rna']: 'rna'}
        if 'expression' in col_map:
            rename_dict[col_map['expression']] = 'expression'
        df = df.rename(columns=rename_dict)
        
        logger.info(f"  Sample mapping loaded using column mapping: {columns}")
        return df
    
    # Fallback: positional interpretation (old behaviour)
    logger.warning("  No sample_map_columns provided. Falling back to positional interpretation (first 2 or 3 columns).")
    # Re-read without header to avoid issues
    df = pd.read_csv(path, sep=sep, header=None)
    
    # Check if first row looks like header
    first_row = df.iloc[0].astype(str)
    if first_row.str.contains('WGS|RNA|ID', case=False, na=False).any():
        logger.info("  Detected header line in sample mapping; removing it.")
        df = df.iloc[1:].reset_index(drop=True)
    
    n_cols = df.shape[1]
    if n_cols < 2:
        raise ValueError(f"Sample mapping file must have at least 2 columns. Found {n_cols}.")
    
    if n_cols > 3:
        logger.warning(f"  Sample mapping file has {n_cols} columns. Using only the first 3 columns: WGS, RNA, Expression (if present).")
        df = df.iloc[:, :3]
        n_cols = 3
    
    if n_cols == 2:
        df.columns = ['rna', 'wgs']
        logger.info("  Sample mapping has 2 columns: RNA, WGS.")
    elif n_cols == 3:
        df.columns = ['wgs', 'rna', 'expression']
        logger.info("  Sample mapping has 3 columns: WGS, RNA, Expression.")
    
    return df

def filter_mapping_by_samples(mapping_df, rna_samples, wgs_samples, expr_samples=None):
    """
    Keep only rows where:
      - rna_sample is in rna_samples
      - wgs_sample is in wgs_samples
      - if expr_samples provided, expression ID is in expr_samples
    Returns filtered DataFrame.
    """
    mask = mapping_df['rna'].isin(rna_samples) & mapping_df['wgs'].isin(wgs_samples)
    if expr_samples is not None and 'expression' in mapping_df.columns:
        mask &= mapping_df['expression'].isin(expr_samples)
    return mapping_df[mask]

def write_filtered_mapping(filtered_df, output_path):
    """Write filtered mapping to TSV (no header)."""
    filtered_df.to_csv(output_path, sep='\t', header=False, index=False)
    logger = get_logger()
    logger.info(f"  Filtered sample mapping written to {output_path}")

def ensure_vcf_index(vcf_path):
    """
    Ensure a VCF file is indexed (create .tbi if missing).
    """
    logger = get_logger()
    if not os.path.exists(vcf_path + '.tbi'):
        logger.info(f"  Indexing {vcf_path}...")
        run_cmd(f"bcftools index {vcf_path}", "Indexing VCF")
    else:
        logger.debug(f"  Index already exists for {vcf_path}")

def is_file_up_to_date(file_path, dependency_paths):
    """
    Check if file exists and is newer than all dependencies.
    Returns True if file is up‑to‑date, False otherwise.
    """
    if not os.path.exists(file_path):
        return False
    file_mtime = os.path.getmtime(file_path)
    for dep in dependency_paths:
        if not os.path.exists(dep):
            return False
        if os.path.getmtime(dep) > file_mtime:
            return False
    return True

# -----------------------------------------------------------------------------
# Expression matrix handling
# -----------------------------------------------------------------------------
def load_expression_matrix(path, sample_list):
    """
    Load gene expression matrix from a CSV file.
    Expected format: rows = samples, columns = genes, first column is sample index (no header name).
    Returns DataFrame with index=sample_id, columns=gene_id, values = expression.
    Only keeps samples present in sample_list.
    """
    df = pd.read_csv(path, index_col=0)  # first column becomes index (sample IDs)
    # Ensure sample IDs are strings (they might be read as objects)
    df.index = df.index.astype(str)
    # Find common samples
    common = [s for s in sample_list if s in df.index]
    if not common:
        raise ValueError("No samples from mapping file found in expression matrix.")
    return df.loc[common]

# -----------------------------------------------------------------------------
# dbSNP mapping for rsID conversion (optional)
# -----------------------------------------------------------------------------
def load_dbsnp_mapping(file_path, build='GRCh38'):
    """
    Load a dbSNP mapping file (rsID to chrom, pos, ref, alt).
    Expected format: tab-separated with header: rsid, chrom, pos, ref, alt
    Returns a DataFrame indexed by rsid.
    """
    df = pd.read_csv(file_path, sep='\t')
    required = ['rsid', 'chrom', 'pos', 'ref', 'alt']
    for col in required:
        if col not in df.columns:
            raise ValueError(f"dbSNP mapping file must contain column: {col}")
    df['pos'] = df['pos'].astype(int)
    return df.set_index('rsid')

# -----------------------------------------------------------------------------
# Publication styling for matplotlib
# -----------------------------------------------------------------------------
def set_publication_style():
    """Apply matplotlib style for publication-quality figures."""
    import matplotlib.pyplot as plt
    plt.rcParams['font.size'] = 12
    plt.rcParams['axes.labelsize'] = 14
    plt.rcParams['axes.titlesize'] = 16
    plt.rcParams['legend.fontsize'] = 10
    plt.rcParams['figure.dpi'] = 300
    plt.rcParams['savefig.dpi'] = 300
    plt.rcParams['savefig.bbox'] = 'tight'