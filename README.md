# RNA‑seq vs WGS Variant QC Pipeline

A modular pipeline for comprehensive quality control of RNA‑seq variant calls using matched whole‑genome sequencing (WGS) data. It validates RNA variants against a trusted DNA reference, detects sample swaps, assesses allele‑specific expression, and evaluates custom SNP panels.

## Features
- **Variant QC** – Compare RNA‑seq VCF with a reference SNP panel (e.g., VerifyBAMID2 10k SNPs)
- **RNA‑WGS Comparison** – Precision, recall, concordance; stratified by region and expression
- **Allele‑Specific Expression (ASE)** – Detect imbalance at heterozygous sites
- **Custom SNP Panels** – Analyse user‑defined lists (e.g., forensic identification SNPs)
- **Publication‑Ready Plots** – PCA, IBS heatmap, Ti/Tv, concordance, and ASE visualisations
- **Comprehensive Report** – Self‑documenting `output_report.tsv` with file descriptions and interpretation
- **Modular & Cached** – Run only the modules you need; intermediate files are reused

## Requirements
- **Linux** or **macOS**
- **Conda** (for environment management)
- At least **8 GB RAM** (16+ GB recommended)
- **50+ GB free disk space** (depends on VCF sizes)

## Installation

1. **Clone or copy** the pipeline to your working directory:
   ```bash
   git clone https://github.com/your-repo/rnaseq-wgs-variant-pipeline.git
   cd rnaseq-wgs-variant-pipeline

2. **Run the setup script** the pipeline to your working directory:
    ```bash
    chmod +x setup.sh
    ./setup.sh

2. **Activate the environment** 
    ```bash
    conda activate rnaseq_wgs_pipeline

3. **Configuration**
    ```yaml
    resource_bed: "/path/to/1000g.phase3.10k.b38.vcf.gz.dat.bed"
    rna_joint_vcf: "/path/to/RNA_jointcall.vcf.gz"
    wgs_joint_vcf: "/path/to/WGS_jointcall.snp.recalibrated.vcf.gz"
    sample_map: "/path/to/sample_mapping.tsv"
    output_dir: "/path/to/results"

    run_modules:
    variant_qc: true
    compare_wgs_rnaseq: true
    study_snp_analysis: true
    exploratory_plotting: true
    .
    .

4. **Running the Pipeline** 
    ```bash
    python pipeline.py --config config.yaml

## Citation
If you use this pipeline in your work, please cite the original publications for the methods used (e.g., GATK, bcftools) and reference this repository.